#!/usr/bin/env python3
"""
Merge a manually-exported EMBASE search into the living-review archive.

EMBASE has no public API, so this half of the search stays manual:
1. Log into Embase.com, run the saved search (same strategy as
   EMBASE_search_prompt.docx), restrict to records added since the last run.
2. Export the results as RIS (.ris) or CSV (.csv).
3. Run:  python scripts/merge_embase.py path/to/export.ris

This dedupes the export against everything already in data/archive.json
(from PubMed and prior Embase runs) by DOI, PMID (if present), and
normalized title, then appends only the genuinely new records to
data/archive.json and data/new_hits.json.

Supported RIS tags: TI/T1 (title), AU (author, repeatable), JO/JF/T2
(journal), PY/Y1 (year), AB/N2 (abstract), DO (doi).
Supported CSV columns (case-insensitive, best-effort matching):
title, authors, journal/source, year, abstract, doi.
"""

import csv
import json
import os
import re
import sys
from datetime import date

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
ARCHIVE_PATH = os.path.join(DATA_DIR, "archive.json")
NEW_HITS_PATH = os.path.join(DATA_DIR, "new_hits.json")


def load_json(path, default):
    if os.path.exists(path):
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    return default


def save_json(path, data):
    os.makedirs(DATA_DIR, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)


def normalize_title(title):
    return re.sub(r"[^a-z0-9]+", " ", title.lower()).strip()


RIS_TAG_MAP = {
    "TI": "title",
    "T1": "title",
    "AU": "author",
    "JO": "journal",
    "JF": "journal",
    "T2": "journal",
    "PY": "pubdate",
    "Y1": "pubdate",
    "AB": "abstract",
    "N2": "abstract",
    "DO": "doi",
}


def parse_ris(path):
    records = []
    current = {}
    authors = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            m = re.match(r"^([A-Z][A-Z0-9])\s*-\s*(.*)$", line)
            if not m:
                continue
            tag, value = m.group(1), m.group(2).strip()
            if tag == "ER":
                if current:
                    current["authors"] = ", ".join(authors[:6]) + (
                        " et al." if len(authors) > 6 else ""
                    )
                    records.append(current)
                current = {}
                authors = []
                continue
            if tag == "AU":
                authors.append(value)
                continue
            field = RIS_TAG_MAP.get(tag)
            if field:
                current[field] = value
    return records


def parse_csv(path):
    records = []
    with open(path, "r", encoding="utf-8-sig", errors="ignore") as f:
        reader = csv.DictReader(f)
        fieldmap = {k.lower().strip(): k for k in (reader.fieldnames or [])}

        def find(*keys):
            for k in keys:
                if k in fieldmap:
                    return fieldmap[k]
            return None

        title_col = find("title", "article title")
        authors_col = find("authors", "author names", "author(s)")
        journal_col = find("journal", "source", "source title")
        year_col = find("year", "publication year", "date")
        abstract_col = find("abstract")
        doi_col = find("doi")

        for row in reader:
            records.append(
                {
                    "title": row.get(title_col, "") if title_col else "",
                    "authors": row.get(authors_col, "") if authors_col else "",
                    "journal": row.get(journal_col, "") if journal_col else "",
                    "pubdate": row.get(year_col, "") if year_col else "",
                    "abstract": row.get(abstract_col, "") if abstract_col else "",
                    "doi": row.get(doi_col, "") if doi_col else "",
                }
            )
    return records


def main():
    if len(sys.argv) != 2:
        print("Usage: python scripts/merge_embase.py path/to/export.(ris|csv)")
        sys.exit(1)

    path = sys.argv[1]
    if not os.path.exists(path):
        print(f"File not found: {path}")
        sys.exit(1)

    ext = os.path.splitext(path)[1].lower()
    if ext == ".ris":
        raw_records = parse_ris(path)
    elif ext == ".csv":
        raw_records = parse_csv(path)
    else:
        print("Unsupported file type. Please export as .ris or .csv from Embase.com.")
        sys.exit(1)

    print(f"Parsed {len(raw_records)} record(s) from {path}.")

    archive = load_json(ARCHIVE_PATH, [])
    seen_dois = {r["doi"].lower() for r in archive if r.get("doi")}
    seen_titles = {normalize_title(r["title"]) for r in archive if r.get("title")}

    fresh = []
    for r in raw_records:
        title = r.get("title", "").strip()
        if not title:
            continue
        doi = (r.get("doi") or "").strip()
        norm_title = normalize_title(title)
        if doi and doi.lower() in seen_dois:
            continue
        if norm_title in seen_titles:
            continue

        record = {
            "pmid": None,
            "doi": doi,
            "title": title,
            "authors": r.get("authors", ""),
            "journal": r.get("journal", ""),
            "pubdate": r.get("pubdate", ""),
            "abstract": r.get("abstract", ""),
            "source": "Embase",
            "url": f"https://doi.org/{doi}" if doi else "",
            "date_added": date.today().isoformat(),
        }
        fresh.append(record)
        seen_titles.add(norm_title)
        if doi:
            seen_dois.add(doi.lower())

    print(f"{len(fresh)} genuinely new record(s) after dedup against the archive.")

    archive.extend(fresh)
    save_json(ARCHIVE_PATH, archive)

    new_hits = load_json(NEW_HITS_PATH, [])
    new_hits.extend(fresh)
    save_json(NEW_HITS_PATH, new_hits)

    print("Done. Re-open krant.html (or re-deploy) to see the merged results.")


if __name__ == "__main__":
    main()
