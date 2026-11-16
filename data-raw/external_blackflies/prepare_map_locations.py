"""Build site-level map inputs from downloaded tables, preserving source coordinates.

Run from any directory with Python 3. Source tables remain unchanged. Repeated
observations at the same coordinates are combined within each study and data
type. Missing/invalid locations are recorded rather than geocoded or guessed.
"""

import csv
import json
import math
import re
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "maps"
FIELDS = ["study", "country", "site", "latitude", "longitude", "data_type",
          "period", "river", "source_rows", "coordinate_note", "needs_review",
          "source_file", "source_doi"]
records, exclusions = [], []


def rows(relative):
    with (ROOT / relative).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.reader(stream))


def clean(value):
    return re.sub(r"\s+", " ", str(value)).strip()


def angle(value):
    """Read hemisphere-tagged degrees/minutes or decimal degrees as printed."""
    value = value.strip().replace(":", "")
    sign = -1 if value.startswith(("S", "W")) else 1
    value = re.sub(r"^[NSEW]\s*", "", value)
    if "°" in value:
        degrees, rest = value.split("°", 1)
        degrees = float(degrees.replace(" ", ""))
        parts = re.findall(r"\d+(?:\.\d+)?", rest)
        minutes = float(parts[0]) if parts else 0
        seconds = float(parts[1]) if len(parts) > 1 else 0
        if minutes >= 60 or seconds >= 60:
            raise ValueError("minutes or seconds outside valid range")
        return sign * (degrees + minutes / 60 + seconds / 3600)
    return sign * float(value.strip())


def coordinate_pair(value):
    match = re.search(r"([NS].*?)\s*([EW].*)", value, re.S)
    if not match:
        raise ValueError("cannot read latitude/longitude pair")
    return angle(match[1]), angle(match[2])


def add(study, country, site, latitude, longitude, data_type, period, river,
        source_file, doi, note="Reported coordinates", needs_review=False):
    try:
        latitude, longitude = float(latitude), float(longitude)
        if not (-90 <= latitude <= 90 and -180 <= longitude <= 180):
            raise ValueError("coordinates outside geographic range")
        if not (math.isfinite(latitude) and math.isfinite(longitude)):
            raise ValueError("coordinates are not finite")
    except ValueError as error:
        exclusions.append(dict(study=study, site=clean(site), reason=str(error),
                               source_file=source_file))
        return
    records.append(dict(zip(FIELDS, [
        study, country, clean(site), latitude, longitude, data_type, clean(period),
        clean(river), 1, note, needs_review, source_file, doi])))


# Sudan and Ethiopia: the source's spaced degree digits are typography only.
source = "sudan_atbara_2019/tables/breeding_sites_table_01.csv"
for row in rows(source)[2:]:
    latitude, longitude = coordinate_pair(row[1])
    add("Sudan: Atbara breeding survey", "Sudan", row[0], latitude, longitude,
        "Fly / vector survey", "2009; 2013; 2015", "Atbara", source,
        "10.1186/s12879-019-4113-1", "Converted from reported degrees and minutes")

source = "ethiopia_northwest_2024/tables/S2_sample_sites_table_01.csv"
for row in rows(source)[1:]:
    add("Ethiopia: Metema / Metekel", "Ethiopia", row[2], row[10], row[11],
        "Fly / vector survey", row[9], row[6], source, "10.1371/journal.pntd.0011868")

source = "south_sudan_naam_yei_2025/tables/article_table_01.csv"
for row in rows(source)[2:]:
    if not row[0] or not row[0][0].isdigit():
        continue
    add("South Sudan: Naam / Yei", "South Sudan", row[2], row[3], row[4],
        "Fly / vector survey", row[0], row[1], source, "10.1371/journal.pntd.0013309")

source = "uganda_tanzania_traps_2017/tables/article_table_01.csv"
for row in rows(source)[1:]:
    add("Uganda / Tanzania: trap study", row[0], row[2], angle(row[3]), angle(row[4]),
        "Fly / vector survey", row[6], row[7], source, "10.1371/journal.pntd.0005688",
        "Converted from reported degrees and minutes")

source = "sudan_seasonality_2016/originals/article.xml"
article_text = " ".join(ET.parse(ROOT / source).getroot().itertext())
for site in ["Nady", "Mograt", "Gumaiza", "Hilat Khatir"]:
    match = re.search(re.escape(site) + r"(?: village)?\s*\((N.*?)\s*,\s*(E.*?)\s*,",
                      article_text)
    if match is None:
        raise ValueError(f"Cannot locate published coordinates for {site}")
    add("Sudan: seasonal biting study", "Sudan", site, angle(match[1]), angle(match[2]),
        "Fly / vector survey", "2007–2010", "Nile" if site in ["Nady", "Mograt"]
        else "Atbara", source, "10.1371/journal.pone.0150309",
        "Converted from the study-area paragraph; source values not corrected")

# DRC PDF supplements: section/year rows are not site observations.
directory = ROOT / "drc_uganda_kakoi_koda_2022/tables"
for path in sorted(directory.glob("S2_species_by_date_locality_page_*_table_01.csv")):
    source = str(path.relative_to(ROOT))
    for row in rows(source):
        if len(row) < 3 or not re.search(r"N\s*\d", row[2]):
            continue
        latitude, longitude = coordinate_pair(row[2])
        add("DRC: species collections", "DRC", row[1], latitude, longitude,
            "Fly / vector survey", row[0], row[1].split("\n")[0], source,
            "10.1371/journal.pntd.0010684")

for path in sorted(directory.glob("S3_crab_trapping_page_*_table_01.csv")):
    source = str(path.relative_to(ROOT))
    for row in rows(source):
        if len(row) < 4 or not re.search(r"N\s*\d", row[2]):
            continue
        latitude, longitude = coordinate_pair(row[2])
        suspect = longitude < 30.7
        note = ("Reported longitude is an outlier among the Kakoi-Koda crab sites; "
                "retained for review, not corrected") if suspect else "Reported coordinates"
        add("DRC: crab / vector survey", "DRC", row[1], latitude, longitude,
            "Fly / vector survey", row[3], row[0], source,
            "10.1371/journal.pntd.0010684", note, suspect)

country, period = "", ""
for path in sorted(directory.glob("S4_crossborder_vector_collections_page_*_table_01.csv")):
    source = str(path.relative_to(ROOT))
    for row in rows(source):
        if len(row) < 5 or not re.match(r"\d+\.\d+/\d+\.\d+", row[4]):
            continue
        country = row[0] or country
        period = row[1] or period
        latitude, longitude = map(float, row[4].split("/"))
        add("DRC / Uganda: cross-border vector survey", country, row[3], latitude,
            longitude, "Fly / vector survey", period, row[2], source,
            "10.1371/journal.pntd.0010684")

# Only geographical and survey metadata are exported from the human records.
source = "drc_kakoi_koda_2026/tables/S1_epidemiological_data_utf8.csv"
with (ROOT / source).open(newline="", encoding="utf-8") as stream:
    for row in csv.DictReader(stream):
        add("DRC: human infection surveys", "DRC", row["Locality"], row["Latitude"],
            row["Longitude"], "Human infection survey", row["Year of study implementation"],
            "", source, "10.1371/journal.pntd.0014164")

source = "gbif_regional_occurrences_2026/tables/georeferenced_simuliidae.csv"
country_names = {"CF": "CAR", "KE": "Kenya", "UG": "Uganda", "SD": "Sudan"}
with (ROOT / source).open(newline="", encoding="utf-8") as stream:
    for row in csv.DictReader(stream):
        name = row["locality"] or row["scientificName"]
        note = "Presence record; source issues: " + (row["issues"] or "none reported")
        if row["coordinateUncertaintyInMeters"]:
            note += "; coordinate uncertainty (m): " + row["coordinateUncertaintyInMeters"]
        add("GBIF: occurrence records", country_names[row["countryCode"]], name,
            row["decimalLatitude"], row["decimalLongitude"], "GBIF presence record",
            row["year"], "", source, "", note)

# Studies without published coordinates in the downloaded tables are explicit.
for study, count, source in [
    ("Chad: Oulibangala / Kaitia", 2, "chad_oulibangala_2026/tables/article_table_01.csv"),
    ("Sudan: Merowe", 11, "sudan_merowe_2014/tables/article_table_01.csv"),
    ("South Sudan: Maridi", 5, "south_sudan_maridi_2023/originals/Seroprevalence_dataset.xlsx"),
    ("Uganda: western foci", 16, "uganda_western_2020/tables/article_table_01.csv"),
]:
    exclusions.append(dict(study=study, site=f"{count} named sites", source_file=source,
                           reason="No site coordinates in the downloaded table; not geocoded"))

# Collapse repeated source observations at one position within the same layer.
groups = defaultdict(list)
for record in records:
    key = (record["study"], record["country"], record["data_type"],
           round(record["latitude"], 6), round(record["longitude"], 6), record["needs_review"])
    groups[key].append(record)
locations = []
for group in groups.values():
    record = group[0].copy()
    for field in ["site", "period", "river", "source_file", "coordinate_note"]:
        record[field] = "; ".join(dict.fromkeys(r[field] for r in group if r[field]))
    record["source_rows"] = len(group)
    locations.append(record)

OUTPUT.mkdir(exist_ok=True)
with (OUTPUT / "locations.csv").open("w", newline="", encoding="utf-8") as stream:
    writer = csv.DictWriter(stream, fieldnames=FIELDS)
    writer.writeheader()
    writer.writerows(locations)
with (OUTPUT / "unmapped_locations.csv").open("w", newline="", encoding="utf-8") as stream:
    writer = csv.DictWriter(stream, fieldnames=["study", "site", "reason", "source_file"])
    writer.writeheader()
    writer.writerows(exclusions)

summary = {"raw_georeferenced_rows": len(records), "mapped_site_entries": len(locations),
           "source_rows_without_valid_coordinates": len(exclusions) - 4,
           "coordinate_review_entries": sum(r["needs_review"] for r in locations),
           "by_study": {study: sum(r["study"] == study for r in locations)
                        for study in sorted({r["study"] for r in locations})}}
(OUTPUT / "location_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2))
