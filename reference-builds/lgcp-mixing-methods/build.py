"""Build the review from its TeX source and the BibTeX embedded in reference.md.

Requires Python 3, pdflatex and bibtex. Intermediate files stay in a temporary
directory; the only generated deliverable is lgcp_methods_review.pdf.
"""

from pathlib import Path
import re
import shutil
import subprocess
import tempfile


def main():
    root = Path(__file__).resolve().parent
    reference = (root / "reference.md").read_text(encoding="utf-8")
    blocks = re.findall(r"```bibtex\n(.*?)\n```", reference, flags=re.DOTALL)
    if len(blocks) != 1:
        raise ValueError("Expected exactly one embedded BibTeX block.")

    stem = "lgcp_methods_review"
    with tempfile.TemporaryDirectory(prefix="lgcp-literature-") as temporary:
        build = Path(temporary)
        (build / "reference.bib").write_text(blocks[0] + "\n", encoding="utf-8")
        shutil.copy2(root / f"{stem}.tex", build)
        latex = ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", f"{stem}.tex"]
        for command in (latex, ["bibtex", stem], latex, latex):
            result = subprocess.run(
                command, cwd=build, text=True, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            )
            if result.returncode:
                raise RuntimeError(result.stdout)

        log = (build / f"{stem}.log").read_text(errors="replace")
        if "undefined" in log.lower() or "Overfull" in log:
            raise RuntimeError("Unresolved references or overflowing text; inspect the TeX source.")
        shutil.copy2(build / f"{stem}.pdf", root / f"{stem}.pdf")
    print(f"Built {root / (stem + '.pdf')}")


if __name__ == "__main__":
    main()
