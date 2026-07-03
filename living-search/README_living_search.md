# The ImmunoRad Times — living search setup

This adds to your **existing** repository — `Immunorad/Immunorad-live-meta-analysis` — no new GitHub project needed.

```
.github/workflows/weekly_pubmed_search.yml   # runs the PubMed search every week
scripts/pubmed_search.py                      # the exact query from your manuscript
scripts/merge_embase.py                       # processes your manual EMBASE export
immunorad_times.html                          # the front page
data/                                         # starts empty, fills in automatically
```

## 1. One-time setup

1. Copy the files/folders above into your existing repo (root level, or adjust paths if your repo is structured differently).
2. Turn on GitHub Pages for the repo (Settings → Pages), so `immunorad_times.html` becomes a real, reachable page. This matters because the page fetches `data/*.json` with a normal `fetch()` call — that only works when served over `https://`, not when the file is opened locally by double-clicking it.
3. Go to the "Actions" tab of your repo. GitHub will pick up the workflow in `.github/workflows/weekly_pubmed_search.yml` automatically. You can click "Run workflow" to trigger it once by hand right away.
4. (Optional) Get a free NCBI API key (https://www.ncbi.nlm.nih.gov/account/settings/) and add it as a repo secret named `NCBI_API_KEY` (Settings → Secrets and variables → Actions). Not required, but it avoids any rate-limit throttling.

From then on, the PubMed search runs **automatically every Monday morning**, using the exact search strategy from your manuscript, and commits new hits to `data/new_hits.json` and `data/archive.json`. Open your `immunorad_times.html` page to see new candidates.

## 2. Adding EMBASE (manual, as agreed)

EMBASE has no public API, so this half stays a manual step on your end:

1. Log into Embase.com, re-run the saved search (same strategy as `EMBASE_search_prompt.docx`), limited to records added since the last time.
2. Export as **RIS** or **CSV**.
3. Send me that file (or paste its contents) in Claude — I'll run `scripts/merge_embase.py` on it, deduplicate against everything already in `data/archive.json` (PubMed + prior EMBASE runs, matched on DOI/PMID/title), and hand you back the updated `data/` files to commit to your repo.

If you'd rather run it yourself locally:
```bash
python scripts/merge_embase.py path/to/your-export.ris
```

## 3. Reading the paper

- **This Week**: unreviewed candidates only (`new_hits.json`). Clicking "Mark as read" keeps the trial visible in the archive but removes the "NEW" flag (remembered per browser via localStorage).
- **Archive**: everything ever found, grouped by the week it was added.
- The page does **not** auto-include anything — that stays your call. It's purely a screening aid: title, authors, journal, abstract, link to the source.

## 4. What this does and doesn't guarantee

- The query in `scripts/pubmed_search.py` is copied **verbatim** from `PubMed_Search_prompt.docx` (Title/Abstract fields, same filters: Humans, English/Dutch, RCT/CCT, from 2010-01-01). Not simplified or "cleaned up."
- New records are detected by **date added to PubMed** (`datetype=edat`), not just publication date — so you won't miss older articles that were only just indexed.
- Deduplication runs on PMID, DOI, and normalized title (fuzzy on spelling/punctuation only, not on content) — so no duplicate trials between PubMed and EMBASE.
- This system screens candidates; it doesn't judge them. Inclusion/exclusion decisions always stay with you.
