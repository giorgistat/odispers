"""Extract downloaded study tables without harmonising their scientific meaning.

Run with a Python environment containing openpyxl and pdfplumber. Originals
remain unchanged. CSVs retain source header rows, categories and missing-value
markers; captions, footnotes and row counts are recorded in tables.json.
"""

import csv
import io
import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path
from zipfile import ZipFile

import openpyxl
import pdfplumber


ROOT = Path(__file__).resolve().parent
WORD = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"


def text(element):
    if element is None:
        return ""
    return re.sub(r"\s+", " ", "".join(element.itertext())).strip()


def save_table(source, name, rows, **notes):
    if not rows:
        raise ValueError(f"Empty table: {source}, {name}")
    output = source.parent.parent / "tables" / f"{name}.csv"
    output.parent.mkdir(exist_ok=True)
    width = max(map(len, rows))
    with output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerows(list(row) + [""] * (width - len(row)) for row in rows)
    return dict(source=str(source.relative_to(ROOT)),
                file=str(output.relative_to(ROOT)), rows_including_headers=len(rows),
                columns=width, **notes)


def jats_tables(source):
    records = []
    root = ET.parse(source).getroot()
    for number, wrapper in enumerate(root.findall(".//table-wrap"), 1):
        table = wrapper.find(".//table")
        if table is None:
            raise ValueError(f"Table has no XML cells: {source}, {number}")
        rows, spanning = [], {}
        for row_number, row in enumerate(table.findall(".//tr")):
            cells = {col: value for (r, col), value in spanning.items()
                     if r == row_number}
            column = 0
            for cell in row:
                if cell.tag not in ("td", "th"):
                    continue
                while column in cells:
                    column += 1
                value = text(cell)
                colspan = int(cell.get("colspan", 1))
                rowspan = int(cell.get("rowspan", 1))
                for offset in range(colspan):
                    cells[column + offset] = value if offset == 0 else ""
                    for future in range(1, rowspan):
                        spanning[row_number + future, column + offset] = (
                            value if offset == 0 else "")
                column += colspan
            rows.append([cells.get(col, "") for col in range(max(cells) + 1)])
        records.append(save_table(
            source, f"article_table_{number:02d}", rows,
            label=text(wrapper.find("label")), caption=text(wrapper.find("caption")),
            footnotes=text(wrapper.find("table-wrap-foot"))))
    return records


def word_tables(source):
    with ZipFile(source) as archive:
        root = ET.fromstring(archive.read("word/document.xml"))
    paragraphs = [" ".join(t.text or "" for t in p.iter(WORD + "t"))
                  for p in root.findall(f"{WORD}body/{WORD}p")]
    records = []
    for number, table in enumerate(root.iter(WORD + "tbl"), 1):
        rows, merged = [], {}
        for row in table.findall(WORD + "tr"):
            values = []
            for cell in row.findall(WORD + "tc"):
                column = len(values)
                value = " ".join(t.text or "" for t in cell.iter(WORD + "t"))
                merge = cell.find(f"{WORD}tcPr/{WORD}vMerge")
                span = cell.find(f"{WORD}tcPr/{WORD}gridSpan")
                width = int(span.get(WORD + "val")) if span is not None else 1
                if merge is not None:
                    if merge.get(WORD + "val") == "restart":
                        merged[column] = value
                    else:
                        value = merged.get(column, value)
                else:
                    merged.pop(column, None)
                values.extend([value] + [""] * (width - 1))
            rows.append(values)
        records.append(save_table(source, f"{source.stem}_table_{number:02d}",
                                  rows, document_notes=paragraphs))
    return records


def spreadsheet_tables(source):
    records = []
    # Read cached values only; the workbook is never edited or recalculated.
    workbook = openpyxl.load_workbook(source, read_only=True, data_only=True)
    for sheet in workbook:
        name = re.sub(r"[^A-Za-z0-9_-]+", "_", sheet.title)
        rows = list(sheet.iter_rows(values_only=True))
        records.append(save_table(source, f"{source.stem}_{name}", rows,
                                  sheet=sheet.title, values="cached workbook values"))
    workbook.close()
    return records


def pdf_tables(source):
    records, pages = [], []
    output = source.parent.parent / "tables"
    output.mkdir(exist_ok=True)
    for previous in output.glob(f"{source.stem}_page_*_table_*.csv"):
        previous.unlink()
    with pdfplumber.open(source) as document:
        for page_number, page in enumerate(document.pages, 1):
            pages.append(f"PAGE {page_number}\n{page.extract_text(layout=True)}")
            # These PLOS supplements draw cell borders as thin rectangles.
            # Default detection also treats text-background boxes as borders.
            vertical = [r["x0"] for r in page.rects
                        if r["width"] < 1 and r["height"] > 10]
            horizontal = [r["top"] for r in page.rects
                          if r["height"] < 1 and r["width"] > 35]
            def distinct_edges(values):
                edges = []
                for value in sorted(values):
                    if not edges or value - edges[-1] > 1:
                        edges.append(value)
                return edges
            columns, boundaries = distinct_edges(vertical), distinct_edges(horizontal)
            if len(columns) >= 2 and len(boundaries) >= 2:
                tables = page.extract_tables({
                    "vertical_strategy": "explicit", "horizontal_strategy": "explicit",
                    "explicit_vertical_lines": columns,
                    "explicit_horizontal_lines": boundaries})
            else:
                tables = page.extract_tables()
            for number, rows in enumerate(tables, 1):
                records.append(save_table(
                    source, f"{source.stem}_page_{page_number:02d}_table_{number:02d}",
                    rows, page=page_number,
                    note="Automatic PDF extraction; line wrapping can split records."))
    (output / f"{source.stem}.txt").write_text("\n\n".join(pages), encoding="utf-8")
    if not records:
        print(f"PDF retained as text only: {source.name}")
    return records


def native_csv(source):
    raw = source.read_bytes()
    try:
        content, encoding = raw.decode("utf-8-sig"), "utf-8-sig"
    except UnicodeDecodeError:
        content, encoding = raw.decode("cp1252"), "cp1252"
    delimiter = ";" if source.name == "S1_epidemiological_data.csv" else ","
    rows = list(csv.reader(io.StringIO(content), delimiter=delimiter))
    return [save_table(source, source.stem + "_utf8", rows,
                       original_encoding=encoding, original_delimiter=delimiter,
                       note="Converted to UTF-8 comma-delimited CSV; cells unchanged.")]


def main():
    records = []
    readers = {".xml": jats_tables, ".docx": word_tables,
               ".xlsx": spreadsheet_tables, ".pdf": pdf_tables, ".csv": native_csv}
    for source in sorted(ROOT.glob("*/originals/*")):
        if source.suffix not in readers:
            continue
        if source.name == "article.pdf":
            continue  # The corresponding XML provides the article tables.
        extracted = readers[source.suffix](source)
        records.extend(extracted)
        print(source.parent.parent.name, source.name, len(extracted), "tables")
    (ROOT / "tables.json").write_text(json.dumps(records, indent=2, ensure_ascii=False)
                                      + "\n", encoding="utf-8")
    print(f"Saved {len(records)} CSV tables; see tables.json for source notes.")


if __name__ == "__main__":
    main()
