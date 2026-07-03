#!/usr/bin/env python3
"""
Import a PubMed CSV export (the standard "Send to > Citation manager" /
Rayyan-style export, with columns PMID, Title, Authors, Journal/Book,
Publication Year, DOI, ...) as *already-known* records into the archive.

Use this when you have results from a previous manual search that should
count as "already seen" — so the automated search never re-surfaces them
as new, even though they were never run through pubmed_search.py itself.

Usage:
    python scripts/import_prior_csv.py path/to/export.csv
"""

import csv
import json
import os
import re
import sys
from datetime import date

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
ARCHIVE_PATH = os.path.join(DATA_DIR, "archive.json")


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


def main():
    if len(sys.argv) != 2:
        print("Usage: python scripts/import_prior_csv.py path/to/export.csv")
        sys.exit(1)

    path = sys.argv[1]
    if not os.path.exists(path):
        print(f"File not found: {path}")
        sys.exit(1)

    archive = load_json(ARCHIVE_PATH, [])
    seen_pmids = {r["pmid"] for r in archive if r.get("pmid")}
    seen_dois = {r["doi"].lower() for r in archive if r.get("doi")}
    seen_titles = {normalize_title(r["title"]) for r in archive if r.get("title")}

    added = []
    with open(path, "r", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            pmid = (row.get("PMID") or "").strip() or None
            doi = (row.get("DOI") or "").strip()
            title = (row.get("Title") or "").strip()
            if not title:
                continue
            norm_title = normalize_title(title)
            if pmid and pmid in seen_pmids:
                continue
            if doi and doi.lower() in seen_dois:
                continue
            if norm_title in seen_titles:
                continue

            record = {
                "pmid": pmid,
                "doi": doi,
                "title": title,
                "authors": (row.get("Authors") or "").strip(),
                "journal": (row.get("Journal/Book") or "").strip(),
                "pubdate": (row.get("Publication Year") or "").strip(),
                "abstract": "",
                "source": "PubMed",
                "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else "",
                "date_added": date.today().isoformat(),
                "prior_search": True,
            }
            added.append(record)
            seen_titles.add(norm_title)
            if pmid:
                seen_pmids.add(pmid)
            if doi:
                seen_dois.add(doi.lower())

    archive.extend(added)
    save_json(ARCHIVE_PATH, archive)
    print(f"Imported {len(added)} new record(s) into the archive as already-known (prior_search).")


if __name__ == "__main__":
    main()
