#!/usr/bin/env python3
"""Build a complete PDF appendix and ZIP archive for Silicon Virovore P4 CSVs.

Run from the Silicon-Virovore repository root:
    python /path/to/build_p4_csv_appendix.py

Outputs:
    prototype4/results/p4_csv_data_appendix.pdf
    prototype4/results/p4_csv_source_archive.zip

The script reads CSV files only. It does not modify source CSVs, SQLite databases,
model settings, or optimizer configuration. Every CSV row and field is included
in the ZIP archive; the PDF renders each CSV record for human review.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT, TA_CENTER
from reportlab.lib.pagesizes import landscape, letter
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak,
    KeepTogether, LongTable
)
from reportlab.pdfbase.pdfmetrics import stringWidth

SKIP_DIRS = {".git", ".venv", "venv", "__pycache__", "node_modules"}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def detect_dialect(path: Path) -> csv.Dialect:
    with path.open("r", encoding="utf-8-sig", newline="", errors="replace") as f:
        sample = f.read(65536)
    try:
        return csv.Sniffer().sniff(sample, delimiters=",;\t|")
    except csv.Error:
        class CommaDialect(csv.excel):
            delimiter = ","
        return CommaDialect()


def inspect_csv(path: Path) -> dict[str, Any]:
    dialect = detect_dialect(path)
    n = 0
    columns: list[str] = []
    with path.open("r", encoding="utf-8-sig", newline="", errors="replace") as f:
        reader = csv.reader(f, dialect)
        try:
            columns = next(reader)
        except StopIteration:
            return {"path": path, "rows": 0, "columns": [], "delimiter": dialect.delimiter}
        for _ in reader:
            n += 1
    return {"path": path, "rows": n, "columns": columns, "delimiter": dialect.delimiter}


def iter_csv_records(path: Path, dialect: csv.Dialect):
    with path.open("r", encoding="utf-8-sig", newline="", errors="replace") as f:
        reader = csv.DictReader(f, dialect=dialect)
        for row in reader:
            yield row


def safe_text(value: Any) -> str:
    if value is None:
        return ""
    if not isinstance(value, str):
        value = str(value)
    # Render control characters visibly without changing source archive bytes.
    return value.replace("\r\n", "\\n").replace("\n", "\\n").replace("\r", "\\r").replace("\t", "\\t")


def chunked(iterable, size: int):
    batch = []
    for item in iterable:
        batch.append(item)
        if len(batch) >= size:
            yield batch
            batch = []
    if batch:
        yield batch


def main() -> int:
    parser = argparse.ArgumentParser(description="Export every CSV under prototype4 to PDF + ZIP")
    parser.add_argument("--repo", default=".", help="Silicon-Virovore repository root (default: current directory)")
    parser.add_argument("--rows-per-page", type=int, default=14,
                        help="Maximum CSV records per planned PDF page (default: 14)")
    parser.add_argument("--min-pages", type=int, default=21,
                        help="Require at least this many PDF pages; default 21")
    args = parser.parse_args()

    repo = Path(args.repo).expanduser().resolve()
    p4root = repo / "prototype4"
    if not p4root.is_dir():
        print(f"ERROR: expected folder not found: {p4root}", file=sys.stderr)
        return 2
    if args.rows_per_page < 1 or args.min_pages < 1:
        print("ERROR: --rows-per-page and --min-pages must be positive.", file=sys.stderr)
        return 2

    csv_paths = []
    for path in p4root.rglob("*"):
        if any(part in SKIP_DIRS for part in path.parts):
            continue
        if path.is_file() and path.suffix.lower() == ".csv":
            csv_paths.append(path)
    csv_paths = sorted(set(csv_paths), key=lambda p: str(p).lower())
    if not csv_paths:
        print(f"ERROR: no CSV files found beneath {p4root}", file=sys.stderr)
        return 2

    output_dir = p4root / "results"
    output_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = output_dir / "p4_csv_data_appendix.pdf"
    zip_path = output_dir / "p4_csv_source_archive.zip"

    print(f"Found {len(csv_paths)} CSV files beneath {p4root}")
    metadata = []
    for path in csv_paths:
        info = inspect_csv(path)
        info["relative"] = path.relative_to(repo).as_posix()
        info["bytes"] = path.stat().st_size
        info["sha256"] = sha256_file(path)
        metadata.append(info)
        print(f"  {info['relative']}: {info['rows']} data rows, {len(info['columns'])} columns")

    total_rows = sum(x["rows"] for x in metadata)
    if total_rows < args.min_pages:
        print(f"ERROR: only {total_rows} CSV records found; cannot honestly generate {args.min_pages} data pages without padding.", file=sys.stderr)
        return 3

    # Dynamically choose a page batch size so datasets with at least min_pages rows
    # naturally produce the requested page count without blank-page padding.
    rows_per_page = min(args.rows_per_page, max(1, total_rows // args.min_pages))
    planned_pages = sum(max(1, (x["rows"] + rows_per_page - 1) // rows_per_page) for x in metadata)
    if planned_pages < args.min_pages:
        print("ERROR: planned pagination is below required minimum; refusing to pad with blank pages.", file=sys.stderr)
        return 3

    # ZIP preserves all original CSV bytes and includes a SHA-256 inventory.
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
        manifest_rows = []
        for info in metadata:
            rel = info["relative"]
            zf.write(info["path"], arcname=rel)
            manifest_rows.append({
                "relative_path": rel,
                "rows": info["rows"],
                "columns": len(info["columns"]),
                "delimiter": info["delimiter"],
                "bytes": info["bytes"],
                "sha256": info["sha256"],
            })
        import io
        buf = io.StringIO(newline="")
        writer = csv.DictWriter(buf, fieldnames=list(manifest_rows[0].keys()))
        writer.writeheader()
        writer.writerows(manifest_rows)
        zf.writestr("CSV_MANIFEST.csv", buf.getvalue().encode("utf-8"))
        zf.writestr("README.txt", (
            "Silicon Virovore P4 CSV source archive\n"
            "All source CSV files beneath prototype4 were copied byte-for-byte.\n"
            "CSV_MANIFEST.csv records relative paths, row/column counts, delimiters, sizes, and SHA-256 hashes.\n"
            "No source CSV was edited by this export.\n"
        ))

    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="CoverTitle2", parent=styles["Title"], fontName="Helvetica-Bold",
                              fontSize=22, leading=27, alignment=TA_CENTER, spaceAfter=14))
    styles.add(ParagraphStyle(name="SectionHead2", parent=styles["Heading1"], fontSize=14, leading=18,
                              textColor=colors.HexColor("#17365D"), spaceBefore=8, spaceAfter=8))
    styles.add(ParagraphStyle(name="SubHead2", parent=styles["Heading2"], fontSize=10, leading=13,
                              textColor=colors.HexColor("#285078"), spaceBefore=5, spaceAfter=4))
    styles.add(ParagraphStyle(name="BodySmall2", parent=styles["BodyText"], fontSize=8, leading=10, spaceAfter=4))
    styles.add(ParagraphStyle(name="Cell2", parent=styles["BodyText"], fontName="Courier", fontSize=5.5,
                              leading=7.0, wordWrap="CJK", splitLongWords=1, spaceAfter=0))
    styles.add(ParagraphStyle(name="CellHead2", parent=styles["BodyText"], fontName="Helvetica-Bold",
                              fontSize=7, leading=8, textColor=colors.white))
    styles.add(ParagraphStyle(name="Note2", parent=styles["BodyText"], fontSize=7, leading=9,
                              textColor=colors.HexColor("#444444")))

    page_size = landscape(letter)
    doc = SimpleDocTemplate(str(pdf_path), pagesize=page_size,
                            rightMargin=0.38*inch, leftMargin=0.38*inch,
                            topMargin=0.48*inch, bottomMargin=0.48*inch,
                            title="Silicon Virovore P4 Complete CSV Data Appendix",
                            author="Silicon Virovore research project")
    story = []
    story.append(Spacer(1, 0.55*inch))
    story.append(Paragraph("Silicon Virovore P4", styles["CoverTitle2"]))
    story.append(Paragraph("Complete CSV Data Appendix", styles["Heading1"]))
    story.append(Spacer(1, 0.15*inch))
    story.append(Paragraph(
        "This appendix renders every CSV file found beneath the prototype4 directory. "
        "The companion ZIP contains byte-for-byte copies of the original CSV files plus a SHA-256 inventory. "
        "The PDF is a human-readable view; use the ZIP as the exact machine-readable source.", styles["BodyText"]))
    story.append(Spacer(1, 0.12*inch))
    story.append(Paragraph(f"Repository: {repo}", styles["BodySmall2"]))
    story.append(Paragraph(f"CSV files: {len(metadata)} | Total data records: {total_rows:,} | "
                           f"Records per planned page: {rows_per_page} | Planned data sections: {planned_pages:,}",
                           styles["BodySmall2"]))
    story.append(Paragraph(f"Generated UTC: {datetime.now(timezone.utc).isoformat(timespec='seconds')}", styles["BodySmall2"]))
    story.append(Spacer(1, 0.15*inch))
    story.append(Paragraph("Source inventory", styles["SectionHead2"]))
    inventory_data = [["CSV file", "Rows", "Columns", "Bytes", "SHA-256 (first 16)"]]
    for info in metadata:
        inventory_data.append([
            Paragraph(info["relative"].replace("&", "&amp;"), styles["BodySmall2"]),
            f"{info['rows']:,}", str(len(info["columns"])), f"{info['bytes']:,}", info["sha256"][:16]
        ])
    inv = Table(inventory_data, colWidths=[4.65*inch, 0.65*inch, 0.65*inch, 0.8*inch, 1.35*inch], repeatRows=1)
    inv.setStyle(TableStyle([
        ("BACKGROUND", (0,0), (-1,0), colors.HexColor("#17365D")),
        ("TEXTCOLOR", (0,0), (-1,0), colors.white),
        ("FONTNAME", (0,0), (-1,0), "Helvetica-Bold"),
        ("FONTSIZE", (0,0), (-1,-1), 7),
        ("GRID", (0,0), (-1,-1), 0.35, colors.HexColor("#B7C9D6")),
        ("VALIGN", (0,0), (-1,-1), "TOP"),
        ("ROWBACKGROUNDS", (0,1), (-1,-1), [colors.white, colors.HexColor("#F2F6F9")]),
        ("LEFTPADDING", (0,0), (-1,-1), 4),
        ("RIGHTPADDING", (0,0), (-1,-1), 4),
        ("TOPPADDING", (0,0), (-1,-1), 3),
        ("BOTTOMPADDING", (0,0), (-1,-1), 3),
    ]))
    story.append(inv)
    story.append(Spacer(1, 0.1*inch))
    story.append(Paragraph(
        "Interpretation note: CSVs may contain source data, derived analysis tables, fold predictions, "
        "or figure-ready graph values. The appendix does not infer which columns were plotted; compare the "
        "listed filenames and columns against each graph's generating script. No values are recalculated.",
        styles["Note2"]))

    for file_index, info in enumerate(metadata, start=1):
        story.append(PageBreak())
        story.append(Paragraph(f"Dataset {file_index}/{len(metadata)}: {info['path'].name}", styles["SectionHead2"]))
        story.append(Paragraph(f"Path: {info['relative']}", styles["Note2"]))
        story.append(Paragraph(f"Rows: {info['rows']:,} | Columns: {len(info['columns'])} | Delimiter: {repr(info['delimiter'])} | "
                               f"SHA-256: {info['sha256']}", styles["Note2"]))
        story.append(Paragraph("Fields: " + ", ".join(str(c) for c in info["columns"]), styles["Note2"]))
        story.append(Spacer(1, 5))

        dialect = detect_dialect(info["path"])
        records = iter_csv_records(info["path"], dialect)
        row_start = 1
        section_index = 0
        for batch in chunked(records, rows_per_page):
            section_index += 1
            row_end = row_start + len(batch) - 1
            story.append(Paragraph(f"Records {row_start:,}-{row_end:,} | section {section_index}", styles["SubHead2"]))
            table_data = [[Paragraph("CSV row", styles["CellHead2"]),
                           Paragraph("Complete record (all columns)", styles["CellHead2"])]]
            for offset, record in enumerate(batch):
                # JSON retains every column/value. Insert spaces after commas for legible wrapping.
                payload = json.dumps(record, ensure_ascii=True, separators=(", ", ": "))
                payload = payload.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
                para = Paragraph(payload, styles["Cell2"])
                table_data.append([str(row_start + offset), para])
            tbl = LongTable(table_data, colWidths=[0.48*inch, 9.75*inch], repeatRows=1, splitByRow=1,
                            hAlign="LEFT")
            tbl.setStyle(TableStyle([
                ("BACKGROUND", (0,0), (-1,0), colors.HexColor("#17365D")),
                ("TEXTCOLOR", (0,0), (-1,0), colors.white),
                ("GRID", (0,0), (-1,-1), 0.25, colors.HexColor("#C4D0D8")),
                ("VALIGN", (0,0), (-1,-1), "TOP"),
                ("ROWBACKGROUNDS", (0,1), (-1,-1), [colors.white, colors.HexColor("#F7F9FA")]),
                ("LEFTPADDING", (0,0), (-1,-1), 3),
                ("RIGHTPADDING", (0,0), (-1,-1), 3),
                ("TOPPADDING", (0,0), (-1,-1), 2),
                ("BOTTOMPADDING", (0,0), (-1,-1), 2),
            ]))
            story.append(tbl)
            row_start = row_end + 1
            if row_start <= info["rows"]:
                story.append(PageBreak())

        if info["rows"] == 0:
            story.append(Paragraph("This CSV contains a header but no data rows.", styles["BodySmall2"]))

    def footer(canvas, document):
        canvas.saveState()
        w, h = page_size
        canvas.setStrokeColor(colors.HexColor("#D4DEE5"))
        canvas.line(0.38*inch, 0.34*inch, w-0.38*inch, 0.34*inch)
        canvas.setFont("Helvetica", 7)
        canvas.setFillColor(colors.HexColor("#58636B"))
        canvas.drawString(0.4*inch, 0.20*inch, "Silicon Virovore P4 | CSV data appendix")
        canvas.drawRightString(w-0.4*inch, 0.20*inch, f"Page {document.page}")
        canvas.restoreState()

    doc.build(story, onFirstPage=footer, onLaterPages=footer)

    # Verify PDF page count using pypdf when available, otherwise use pdfinfo.
    page_count = None
    try:
        from pypdf import PdfReader
        page_count = len(PdfReader(str(pdf_path)).pages)
    except Exception:
        try:
            import subprocess
            out = subprocess.check_output(["pdfinfo", str(pdf_path)], text=True, stderr=subprocess.DEVNULL)
            for line in out.splitlines():
                if line.startswith("Pages:"):
                    page_count = int(line.split(":", 1)[1].strip())
                    break
        except Exception:
            page_count = None

    print("\nEXPORT COMPLETE")
    print("PDF:", pdf_path)
    print("ZIP:", zip_path)
    print("CSV files:", len(metadata))
    print("CSV data rows:", total_rows)
    print("Records per planned section:", rows_per_page)
    print("Planned data sections:", planned_pages)
    print("Verified PDF pages:", page_count if page_count is not None else "page-count utility unavailable")
    if page_count is not None and page_count < args.min_pages:
        print(f"ERROR: PDF has {page_count} pages; required at least {args.min_pages}.", file=sys.stderr)
        return 4
    print("Source data files and SQLite databases were not modified.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
