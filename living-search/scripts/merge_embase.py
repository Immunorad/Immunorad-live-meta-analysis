#!/usr/bin/env python3
"""
Merge a manually-exported EMBASE search into the living-review archive.

EMBASE has no public API, so this half of the search stays manual:
1. Log into Embase.com, run the saved search (same strategy as
   EMBASE_search_prompt.docx), restrict to records added since the last run.
2. Export the results — RIS, a plain CSV, or Embase.com's native
   "complete reference" CSV export (one field per row, e.g. "TITLE",
   "AUTHOR NAMES", "ABSTRACT", "DOI", "MEDLINE PMID", ...).
3. Run:  python scripts/merge_embase.py path/to/export.(ris|csv)

This dedupes the export against everything already in data/archive.json
(from PubMed and prior Embase runs) by PMID (if the Embase record is also
MEDLINE-indexed), DOI, and normalized title, then appends only the
genuinely new records to data/archive.json and data/new_hits.json.

Supported formats (auto-detected):
- RIS: TI/T1 (title), AU (author, repeatable), JO/JF/T2 (journal),
  PY/Y1 (year), AB/N2 (abstract), DO (doi).
- Plain CSV: columns title, authors, journal/source, year, abstract, doi
  (case-insensitive, best-effort matching).
- Embase.com native CSV export: one field per row, first column is the
  field name (TITLE, AUTHOR NAMES, SOURCE TITLE, PUBLICATION YEAR,
  ABSTRACT, DOI, MEDLINE PMID, EMBASE LINK, ...), records separated by
  a blank line. Detected automatically if the file starts with a
  "SEARCH QUERY" row.
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


# ---------------------------------------------------------------------------
# Format 1: RIS
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Format 2: plain CSV (one row per record, normal header row)
# ---------------------------------------------------------------------------

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
        pmid_col = find("pmid", "medline pmid")

        for row in reader:
            pmid_val = (row.get(pmid_col, "") if pmid_col else "").strip()
            records.append(
                {
                    "title": row.get(title_col, "") if title_col else "",
                    "authors": row.get(authors_col, "") if authors_col else "",
                    "journal": row.get(journal_col, "") if journal_col else "",
                    "pubdate": row.get(year_col, "") if year_col else "",
                    "abstract": row.get(abstract_col, "") if abstract_col else "",
                    "doi": row.get(doi_col, "") if doi_col else "",
                    "pmid": pmid_val if pmid_val.isdigit() else None,
                }
            )
    return records


# ---------------------------------------------------------------------------
# Format 3: Embase.com native "field-per-row" export
# First column of each row is the field name (TITLE, AUTHOR NAMES, ...),
# remaining columns are that field's value(s). Records are separated by a
# blank line. Detected by a leading "SEARCH QUERY" row.
# ---------------------------------------------------------------------------

def is_embase_native_format(path):
    with open(path, "r", encoding="utf-8-sig", errors="ignore", newline="") as f:
        reader = csv.reader(f)
        for row in reader:
            if row and row[0].strip():
                return row[0].strip() == "SEARCH QUERY"
    return False


def parse_embase_native(path):
    raw_records = []
    current = {}
    with open(path, "r", encoding="utf-8-sig", errors="ignore", newline="") as f:
        reader = csv.reader(f)
        for row in reader:
            if not row or all(not c.strip() for c in row):
                if current.get("TITLE"):
                    raw_records.append(current)
                current = {}
                continue
            field = row[0].strip()
            if field in ("SEARCH QUERY",) or field.startswith("---"):
                continue
            values = [v for v in row[1:] if v.strip()]
            if not values:
                continue
            if field in current:
                current[field].extend(values)
            else:
                current[field] = values
        if current.get("TITLE"):
            raw_records.append(current)

    records = []
    for r in raw_records:
        title = r.get("TITLE", [""])[0]
        authors_list = r.get("AUTHOR NAMES", [])
        authors = ", ".join(authors_list[:6]) + (" et al." if len(authors_list) > 6 else "")
        journal = r.get("SOURCE TITLE", [""])[0].strip()
        pubdate = r.get("PUBLICATION YEAR", [""])[0]
        abstract = r.get("ABSTRACT", [""])[0]
        doi = r.get("DOI", [""])[0].strip()
        pmid_raw = r.get("MEDLINE PMID", [""])
        pmid = pmid_raw[0].strip() if pmid_raw and pmid_raw[0].strip().isdigit() else None

        records.append(
            {
                "title": title,
                "authors": authors,
                "journal": journal,
                "pubdate": pubdate,
                "abstract": abstract,
                "doi": doi,
                "pmid": pmid,
            }
        )
    return records


# ---------------------------------------------------------------------------
# Merge logic
# ---------------------------------------------------------------------------

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
        fmt = "RIS"
    elif ext == ".csv":
        if is_embase_native_format(path):
            raw_records = parse_embase_native(path)
            fmt = "Embase native CSV"
        else:
            raw_records = parse_csv(path)
            fmt = "plain CSV"
    else:
        print("Unsupported file type. Please export as .ris or .csv from Embase.com.")
        sys.exit(1)

    print(f"Parsed {len(raw_records)} record(s) from {path} (detected format: {fmt}).")

    state = load_json(STATE_PATH, {"last_run": date.today().isoformat(), "seen_pmids": []})
    state["total_records_identified"] = state.get("total_records_identified", 0) + len(raw_records)
    save_json(STATE_PATH, state)

    archive = load_json(ARCHIVE_PATH, [])
    seen_pmids = {r["pmid"] for r in archive if r.get("pmid")}
    seen_dois = {r["doi"].lower() for r in archive if r.get("doi")}
    seen_titles = {normalize_title(r["title"]) for r in archive if r.get("title")}

    fresh = []
    for r in raw_records:
        title = (r.get("title") or "").strip()
        if not title:
            continue
        doi = (r.get("doi") or "").strip()
        pmid = r.get("pmid")
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
            "authors": r.get("authors", ""),
            "journal": r.get("journal", ""),
            "pubdate": r.get("pubdate", ""),
            "abstract": r.get("abstract", ""),
            "source": "Embase",
            "url": (
                f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
                if pmid
                else (f"https://doi.org/{doi}" if doi else "")
            ),
            "date_added": date.today().isoformat(),
        }
        fresh.append(record)
        seen_titles.add(norm_title)
        if doi:
            seen_dois.add(doi.lower())
        if pmid:
            seen_pmids.add(pmid)

    print(f"{len(fresh)} genuinely new record(s) after dedup against the archive "
          f"(by PMID, DOI, and normalized title).")

    archive.extend(fresh)
    save_json(ARCHIVE_PATH, archive)

    new_hits = load_json(NEW_HITS_PATH, [])
    new_hits.extend(fresh)
    save_json(NEW_HITS_PATH, new_hits)

    print("Done. Re-open immunorad_times.html (or re-deploy) to see the merged results.")


if __name__ == "__main__":
    main()
