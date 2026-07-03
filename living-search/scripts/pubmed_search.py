#!/usr/bin/env python3
"""
Living meta-analysis: automated PubMed search.

Runs the EXACT literal PubMed search strategy from the manuscript
(Manuscript_EJC_immunorad_clear_13052026.docx / PubMed_Search_prompt.docx),
retrieves any records newly indexed in PubMed since the last run, and
writes them to data/new_hits.json (this run's new candidates) while
appending everything to data/archive.json (full history, for dedup).

No API key is required for light use (< 3 requests/second), but you can
set an NCBI_API_KEY repo secret / environment variable to raise the rate
limit and reduce the chance of throttling.

Usage:
    python scripts/pubmed_search.py
"""

import json
import os
import re
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from datetime import date, datetime, timedelta, timezone

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
STATE_PATH = os.path.join(DATA_DIR, "state.json")
ARCHIVE_PATH = os.path.join(DATA_DIR, "archive.json")
NEW_HITS_PATH = os.path.join(DATA_DIR, "new_hits.json")

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
API_KEY = os.environ.get("NCBI_API_KEY", "").strip()
TOOL_EMAIL = os.environ.get("NCBI_EMAIL", "")  # optional, NCBI etiquette

# ---------------------------------------------------------------------------
# EXACT literal search strategy, transcribed field-for-field from
# PubMed_Search_prompt.docx. Do not "clean up" or paraphrase this — it must
# match the original systematic review search 1:1.
# ---------------------------------------------------------------------------

DOMAIN = (
    '(Cancer*[Title/Abstract] OR tumor*[Title/Abstract] OR tumour*[Title/Abstract] '
    'OR neoplasm*[Title/Abstract] OR malignan*[Title/Abstract] OR carcinom*[Title/Abstract] '
    'OR adenocarcinom*[Title/Abstract] OR sarcom*[Title/Abstract] OR squamous[Title/Abstract] '
    'OR SCC[Title/Abstract] OR Glio*[Title/Abstract] OR Meningiom*[Title/Abstract] '
    'OR Neuroblastom*[Title/Abstract] OR Medulloblastom*[Title/Abstract] OR Schwannom*[Title/Abstract] '
    'OR Melanom*[Title/Abstract] OR Retinoblastom*[Title/Abstract] OR Mesotheliom*[Title/Abstract] '
    'OR "Solid tumor*"[Title/Abstract] OR "Solid tumour*"[Title/Abstract] OR NSCLC[Title/Abstract])'
)

DETERMINANT_1_RADIOTHERAPY = (
    '(Radiotherap*[Title/Abstract] OR Radiat*[Title/Abstract] OR Irradiat*[Title/Abstract] '
    'OR "X-ray-therap*"[Title/Abstract] OR "X-ray therap*"[Title/Abstract] OR chemoradiat*[Title/Abstract] '
    'OR Chemoradio*[Title/Abstract] OR "Chemo-radio*"[Title/Abstract] OR SBRT[Title/Abstract] '
    'OR SABR[Title/Abstract] OR SRT[Title/Abstract] OR SRS[Title/Abstract] OR "Gamma knife"[Title/Abstract] '
    'OR GKRS[Title/Abstract] OR LINAC[Title/Abstract] OR CyberKnife[Title/Abstract] '
    'OR radiosurgery[Title/Abstract] OR stereostatic[Title/Abstract])'
)

DETERMINANT_2_IMMUNOTHERAPY = (
    '(Immunotherap*[Title/Abstract] OR "Immuno-therap*"[Title/Abstract] OR Chemoimmunotherap*[Title/Abstract] '
    'OR Immunochemotherap*[Title/Abstract] OR "Immuno-chemotherap*"[Title/Abstract] '
    'OR "Chemo-immunotherap*"[Title/Abstract] OR "Checkpoint inhibitor*"[Title/Abstract] '
    'OR "Checkpoint-Inhibitor*"[Title/Abstract] OR "Immune checkpoint inhibitor*"[Title/Abstract] '
    'OR "Immune-checkpoint inhibitor*"[Title/Abstract] OR "Immune-checkpoint-inhibitor*"[Title/Abstract] '
    'OR "CTLA-4"[Title/Abstract] OR "PD-1"[Title/Abstract] OR "PD-L1"[Title/Abstract] '
    'OR "Immuno checkpoint inhibitor*"[Title/Abstract] OR "Immuno-checkpoint inhibitor*"[Title/Abstract] '
    'OR "Immuno-checkpoint-inhibitor*"[Title/Abstract] OR "Immuno checkpoint-inhibitor*"[Title/Abstract] '
    'OR bevacizumab[Title/Abstract] OR pembrolizumab[Title/Abstract] OR durvalumab[Title/Abstract] '
    'OR ipilimumab[Title/Abstract])'
)

DETERMINANT_3_RADIOIMMUNOTHERAPY = (
    '(Radioimmunotherap*[Title/Abstract] OR "Radio-immunotherap*"[Title/Abstract] '
    'OR Immunoradiotherap*[Title/Abstract])'
)

OUTCOME = (
    '("Overall survival"[Title/Abstract] OR Mortality[Title/Abstract] OR Surviv*[Title/Abstract])'
)

DETERMINANT = f'(({DETERMINANT_1_RADIOTHERAPY} AND {DETERMINANT_2_IMMUNOTHERAPY}) OR {DETERMINANT_3_RADIOIMMUNOTHERAPY})'

FILTERS = (
    'AND (Humans[Filter]) '
    'AND (English[Language] OR Dutch[Language]) '
    'AND (Randomized Controlled Trial[Publication Type] OR Controlled Clinical Trial[Publication Type])'
)

# Lower date bound fixed at the original review's start date (2010-01-01).
# The upper bound / "what's new" is handled separately via datetype=edat
# (date a record was ADDED to PubMed), so we catch newly-indexed older
# articles too, not just newly-published ones.
BASE_QUERY = f'{DOMAIN} AND {DETERMINANT} AND {OUTCOME} {FILTERS} AND ("2010/01/01"[Date - Publication] : "3000"[Date - Publication])'


def _eutils_params(extra):
    params = dict(extra)
    if API_KEY:
        params["api_key"] = API_KEY
    if TOOL_EMAIL:
        params["email"] = TOOL_EMAIL
    params["tool"] = "immunorad-living-review"
    return params


def _get(url, params, retries=3):
    query = urllib.parse.urlencode(_eutils_params(params))
    full_url = f"{url}?{query}"
    last_err = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(full_url, timeout=30) as resp:
                return resp.read()
        except Exception as e:  # noqa: BLE001
            last_err = e
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"Failed to fetch {url}: {last_err}")


def load_state():
    if os.path.exists(STATE_PATH):
        with open(STATE_PATH, "r", encoding="utf-8") as f:
            return json.load(f)
    return {"last_run": "2026-01-15", "seen_pmids": []}


def save_state(state):
    os.makedirs(DATA_DIR, exist_ok=True)
    with open(STATE_PATH, "w", encoding="utf-8") as f:
        json.dump(state, f, indent=2, ensure_ascii=False)


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


def esearch_new_pmids(mindate, maxdate, datetype="edat"):
    """Run the literal query, restricted to records added to PubMed (edat)
    or published (pdat) in the given window, depending on datetype."""
    all_ids = []
    retstart = 0
    retmax = 200
    while True:
        raw = _get(
            f"{EUTILS}/esearch.fcgi",
            {
                "db": "pubmed",
                "term": BASE_QUERY,
                "datetype": datetype,
                "mindate": mindate,
                "maxdate": maxdate,
                "retmode": "json",
                "retstart": retstart,
                "retmax": retmax,
                "sort": "most recent",
            },
        )
        payload = json.loads(raw)
        idlist = payload.get("esearchresult", {}).get("idlist", [])
        all_ids.extend(idlist)
        count = int(payload.get("esearchresult", {}).get("count", 0))
        retstart += retmax
        if retstart >= count or not idlist:
            break
        time.sleep(0.4)
    return all_ids


def efetch_details(pmids):
    """Fetch title/authors/journal/pubdate/abstract/doi for a batch of PMIDs."""
    if not pmids:
        return []
    records = []
    batch_size = 100
    for i in range(0, len(pmids), batch_size):
        batch = pmids[i : i + batch_size]
        raw = _get(
            f"{EUTILS}/efetch.fcgi",
            {
                "db": "pubmed",
                "id": ",".join(batch),
                "retmode": "xml",
                "rettype": "abstract",
            },
        )
        root = ET.fromstring(raw)
        for article in root.findall(".//PubmedArticle"):
            pmid_el = article.find(".//PMID")
            pmid = pmid_el.text if pmid_el is not None else None

            title_el = article.find(".//ArticleTitle")
            title = "".join(title_el.itertext()).strip() if title_el is not None else "(no title)"

            journal_el = article.find(".//Journal/Title")
            journal = journal_el.text if journal_el is not None else ""

            year_el = article.find(".//JournalIssue/PubDate/Year")
            medline_date_el = article.find(".//JournalIssue/PubDate/MedlineDate")
            pubdate = (
                year_el.text
                if year_el is not None
                else (medline_date_el.text if medline_date_el is not None else "")
            )

            authors = []
            for author in article.findall(".//AuthorList/Author"):
                last = author.find("LastName")
                initials = author.find("Initials")
                if last is not None:
                    name = last.text
                    if initials is not None:
                        name += f" {initials.text}"
                    authors.append(name)
            authors_str = ", ".join(authors[:6]) + (" et al." if len(authors) > 6 else "")

            abstract_parts = [
                "".join(el.itertext()) for el in article.findall(".//Abstract/AbstractText")
            ]
            abstract = " ".join(abstract_parts).strip()

            doi = ""
            for el_id in article.findall(".//ArticleIdList/ArticleId"):
                if el_id.attrib.get("IdType") == "doi":
                    doi = el_id.text

            records.append(
                {
                    "pmid": pmid,
                    "doi": doi,
                    "title": title,
                    "authors": authors_str,
                    "journal": journal,
                    "pubdate": pubdate,
                    "abstract": abstract,
                    "source": "PubMed",
                    "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                    "date_added": date.today().isoformat(),
                }
            )
        time.sleep(0.4)
    return records


def dedupe_against_archive(records, archive):
    """Drop anything already present by PMID, DOI, or normalized title."""
    seen_pmids = {r["pmid"] for r in archive if r.get("pmid")}
    seen_dois = {r["doi"].lower() for r in archive if r.get("doi")}
    seen_titles = {normalize_title(r["title"]) for r in archive if r.get("title")}

    fresh = []
    for r in records:
        if r.get("pmid") and r["pmid"] in seen_pmids:
            continue
        if r.get("doi") and r["doi"].lower() in seen_dois:
            continue
        if normalize_title(r["title"]) in seen_titles:
            continue
        fresh.append(r)
        seen_pmids.add(r.get("pmid"))
        if r.get("doi"):
            seen_dois.add(r["doi"].lower())
        seen_titles.add(normalize_title(r["title"]))
    return fresh


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Living meta-analysis PubMed search")
    parser.add_argument(
        "--since",
        default=None,
        help="Override the start date (YYYY-MM-DD). Use for a one-time backfill, "
        "e.g. --since 2025-01-01. If omitted, uses the last run date from state.json "
        "(normal weekly incremental behaviour).",
    )
    parser.add_argument(
        "--datetype",
        default="edat",
        choices=["edat", "pdat"],
        help="edat = date added to PubMed (default, best for catching new hits). "
        "pdat = publication date (use for a one-time backfill, to also catch "
        "records PubMed indexed late under an older publication date).",
    )
    args = parser.parse_args()

    state = load_state()
    archive = load_json(ARCHIVE_PATH, [])

    mindate = (args.since or state["last_run"]).replace("-", "/")
    maxdate = date.today().isoformat().replace("-", "/")

    print(f"Searching PubMed ({args.datetype}) for records between {mindate} and {maxdate} ...")
    pmids = esearch_new_pmids(mindate, maxdate, datetype=args.datetype)
    print(f"esearch returned {len(pmids)} candidate PMIDs.")

    records = efetch_details(pmids)
    fresh = dedupe_against_archive(records, archive)
    print(f"{len(fresh)} genuinely new record(s) after dedup against the archive.")

    # Merge into archive (full history) and write this run's new hits.
    archive.extend(fresh)
    save_json(ARCHIVE_PATH, archive)

    existing_new_hits = load_json(NEW_HITS_PATH, [])
    # new_hits.json accumulates *unreviewed* candidates across runs; the
    # "krant" front-end / your screening step is responsible for clearing
    # entries out of it once reviewed.
    combined_new = existing_new_hits + fresh
    save_json(NEW_HITS_PATH, combined_new)

    # Only advance last_run for the normal incremental (edat) mode. A manual
    # --since backfill on pdat shouldn't change what next Monday's regular
    # run considers "already covered".
    if args.since is None:
        state["last_run"] = date.today().isoformat()
        save_state(state)

    print("Done.")


if __name__ == "__main__":
    main()
