#!/usr/bin/env python3
"""
Import one or more historical search export CSVs (PubMed-style or
Embase-style, auto-detected) as the baseline record set for the living
review, and update the cumulative "records identified" counter used by
the automated PRISMA flow diagram.

Unlike a weekly/incremental run, this is meant for bulk historical
imports — e.g. the full 2010-2025 and 2010-2026 vintage exports used to
write the manuscript. Records are merged into data/archive.json (deduped
by PMID / DOI / normalized title) and tagged "prior_search": true, so
they never appear as new candidates in "To Screen". They also do NOT
get added to data/new_hits.json.

Every raw row across every file passed in this run is added to
state.json's cumulative "total_records_identified" counter — this
represents "records identified" for the PRISMA diagram, deliberately
including records that turn out to be duplicates, since that's the
whole point of tracking "duplicates removed" as a separate, derived
number (total_records_identified - len(archive.json)).

Usage:
    python scripts/import_historical.py file1.csv file2.csv ...

Auto-detected formats:
- PubMed CSV export: has a "PMID" column.
- Embase CSV export (row-per-record variant): has "Author Names" and/or
  "AiP/IP Entry Date" columns, no PMID column, no abstract.
"""

import csv
import json
import os
import re
import sys
from datetime import date

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
ARCHIVE_PATH = os.path.join(DATA_DIR, "archive.json")
STATE_PATH = os.path.join(DATA_DIR, "state.json")


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


def detect_format(fieldnames):
    lower = [f.lower().strip() for f in (fieldnames or [])]
    if "pmid" in lower:
        return "pubmed"
    if "author names" in lower or "aip/ip entry date" in lower:
        return "embase"
    return "unknown"


def parse_pubmed_csv(path):
    records = []
    with open(path, "r", encoding="utf-8-sig", errors="ignore") as f:
        reader = csv.DictReader(f)
        for row in reader:
            title = (row.get("Title") or "").strip()
            if not title:
                continue
            pmid = (row.get("PMID") or "").strip() or None
            records.append(
                {
                    "pmid": pmid,
                    "doi": (row.get("DOI") or "").strip(),
                    "title": title,
                    "authors": (row.get("Authors") or "").strip(),
                    "journal": (row.get("Journal/Book") or "").strip(),
                    "pubdate": (row.get("Publication Year") or "").strip(),
                    "abstract": "",
                    "source": "PubMed",
                    "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else "",
                }
            )
    return records


def parse_embase_csv(path):
    records = []
    with open(path, "r", encoding="utf-8-sig", errors="ignore") as f:
        reader = csv.DictReader(f)
        for row in reader:
            title = (row.get("Title") or "").strip()
            if not title:
                continue
            records.append(
                {
                    "pmid": None,
                    "doi": (row.get("DOI") or "").strip(),
                    "title": title,
                    "authors": (row.get("Author Names") or "").strip(),
                    "journal": (row.get("Source title") or "").strip(),
                    "pubdate": (row.get("Publication Year") or "").strip(),
                    "abstract": "",
                    "source": "Embase",
                    "url": "",
                }
            )
    return records


def main():
    if len(sys.argv) < 2:
        print("Usage: python scripts/import_historical.py file1.csv [file2.csv ...]")
        sys.exit(1)

    paths = sys.argv[1:]
    for p in paths:
        if not os.path.exists(p):
            print(f"File not found: {p}")
            sys.exit(1)

    archive = load_json(ARCHIVE_PATH, [])
    state = load_json(STATE_PATH, {"last_run": date.today().isoformat(), "seen_pmids": []})

    total_raw_this_run = 0
    all_new = []

    for path in paths:
        with open(path, "r", encoding="utf-8-sig", errors="ignore") as f:
            fieldnames = csv.DictReader(f).fieldnames
        fmt = detect_format(fieldnames)

        if fmt == "pubmed":
            raw_records = parse_pubmed_csv(path)
        elif fmt == "embase":
            raw_records = parse_embase_csv(path)
        else:
            print(f"Could not detect format for {path} (columns: {fieldnames}). Skipping.")
            continue

        total_raw_this_run += len(raw_records)
        print(f"{os.path.basename(path)}: detected {fmt}, {len(raw_records)} raw record(s).")

        seen_pmids = {r["pmid"] for r in archive if r.get("pmid")}
        seen_dois = {r["doi"].lower() for r in archive if r.get("doi")}
        seen_titles = {normalize_title(r["title"]) for r in archive if r.get("title")}

        for r in raw_records:
            norm_title = normalize_title(r["title"])
            if r.get("pmid") and r["pmid"] in seen_pmids:
                continue
            if r.get("doi") and r["doi"].lower() in seen_dois:
                continue
            if norm_title in seen_titles:
                continue

            record = dict(r)
            record["date_added"] = date.today().isoformat()
            record["prior_search"] = True
            archive.append(record)
            all_new.append(record)
            seen_titles.add(norm_title)
            if r.get("doi"):
                seen_dois.add(r["doi"].lower())
            if r.get("pmid"):
                seen_pmids.add(r["pmid"])

    print(f"\nTotal raw records across all files this run: {total_raw_this_run}")
    print(f"Genuinely new unique records added to the archive: {len(all_new)}")
    print(f"Duplicates in this run (within files + against existing archive): "
          f"{total_raw_this_run - len(all_new)}")

    save_json(ARCHIVE_PATH, archive)

    state["total_records_identified"] = state.get("total_records_identified", 0) + total_raw_this_run
    save_json(STATE_PATH, state)

    print(f"\nCumulative total_records_identified is now: {state['total_records_identified']}")
    print(f"Cumulative archive size (unique records): {len(archive)}")
    print(f"=> Cumulative duplicates removed: {state['total_records_identified'] - len(archive)}")


if __name__ == "__main__":
    main()
