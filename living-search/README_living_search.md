# The ImmunoRad Times — living search

Lives inside `Immunorad/Immunorad-live-meta-analysis`, in its own `living-search/` folder, so it never touches your existing site files.

```
living-search/
├── immunorad_times.html          # the front page (This Week / Included / Excluded / Archive)
├── scripts/
│   ├── pubmed_search.py           # exact query from the manuscript; supports weekly + one-off backfill
│   ├── merge_embase.py            # processes a manual EMBASE export (RIS/CSV)
│   └── import_prior_csv.py        # imports a prior PubMed/Rayyan CSV export as "already known"
└── data/
    ├── state.json                 # last run date
    ├── archive.json               # every record ever found (full history)
    ├── new_hits.json              # candidates awaiting a decision
    └── decisions.json             # shared Include/Exclude decisions (see below)

.github/workflows/
├── immunorad_living_search.yml    # runs every Monday automatically
└── immunorad_backfill_2025.yml    # manual, one-off: re-search from 2025-01-01
```

Already deployed at: `https://immunorad.github.io/Immunorad-live-meta-analysis/living-search/immunorad_times.html`

## Screening: Include / Exclude, shared across the team

Every candidate has **Include** and **Exclude** buttons. Once clicked, the trial moves out of "This Week" into the **Included** or **Excluded** tab, with who decided and when.

These decisions are **not** stored per-browser — they're written directly into `living-search/data/decisions.json` in the GitHub repo, via the GitHub API, so every colleague sees the same shared state and there's a full audit trail (every decision is a git commit).

### For each colleague who wants to screen

1. They need **write access to the repo** (add them as a collaborator: repo → Settings → Collaborators, or they're already in the `Immunorad` organization).
2. They create a **fine-grained personal access token**: GitHub → their profile → Settings → Developer settings → Personal access tokens → Fine-grained tokens → **Generate new token** → Repository access: **only this repo** → Permissions: **Contents: Read and write**.
3. On the ImmunoRad Times page, they click **"Screening access"** (top right), paste their name and that token, and click Save. It's stored only in their own browser (localStorage) and only ever sent to `api.github.com`.
4. They can now click Include/Exclude and it saves for everyone.

No token = read-only: they can browse and read abstracts, but clicking Include/Exclude will prompt them to add access first.

## One-time backfill from 2025-01-01

You asked for this because PubMed sometimes indexes older records late, so a pure "what's new" search could miss things the original 15-01-2026 search didn't catch either.

**Step 1 — import your prior CSV results as a known baseline:**
Your `Pubmed_2026_NEW_vs_2025_RAYYAN.csv` (236 records) is committed at `living-search/data/prior_imports/`. Run:
**Actions → "ImmunoRad import prior CSV" → Run workflow**
This adds all 236 records to `data/archive.json` tagged `"prior_search": true`, on top of whatever the live archive currently contains — so it never conflicts with what the weekly Action has already found, and none of these 236 will ever reappear as "new".

**Step 2 — run the actual backfill:**
**Actions → "ImmunoRad one-time backfill (since 2025-01-01)" → Run workflow**
This runs the exact literal query with `mindate=2025-01-01`, `datetype=pdat` (publication date, to also catch late-indexed older records), dedupes against the full archive (including the 236 just imported), and adds only genuinely new hits to `new_hits.json`. It's a one-off — it does not affect the normal Monday schedule.

Run Step 1 before Step 2, so the dedup already knows about your prior results.

## Why I couldn't just run the search myself in this chat

My tools here can't reach the PubMed API directly (network restrictions in this environment) and can't call arbitrary URLs I haven't been given or found via search. GitHub Actions runs on GitHub's own servers, which have normal internet access — so that's where the actual API calls, and the CSV import, happen. This also avoids merge conflicts: the Action always works from whatever the live `archive.json` currently is, instead of a local copy that can drift out of sync.

## Adding EMBASE (manual, as before)

1. Export your saved Embase.com search as RIS or CSV.
2. Send it to me — I'll run `scripts/merge_embase.py` on it and hand back updated `data/` files, or run it yourself:
   ```bash
   python living-search/scripts/merge_embase.py path/to/export.ris
   ```

## Importing another prior CSV later

```bash
python living-search/scripts/import_prior_csv.py path/to/export.csv
```
Marks everything in it as already-known (same effect as what was done with your Rayyan export), without touching `new_hits.json`.

## Automated PRISMA flow diagram

A new **PRISMA** tab on the page shows a live flow diagram:

- **Records identified**, **Duplicates removed**, **Unique records screened**, **Excluded on title/abstract** — fully automatic, computed from `state.json`'s cumulative `total_records_identified` counter and `archive.json`/`decisions.json`.
- **"Eligible for full text"** (renamed from "Include") — the count of trials marked eligible, also automatic.
- **Eligibility (full-text) exclusions** and **Inclusion** (final trial count, patient count, subgroups) — these require human full-text judgment, so they're **not** automatic. Fill them in yourself in `data/manual_prisma_stats.json` once full-text review for a round is done; the page picks them up automatically.

### Importing the historical 2025/2026 search vintages

Your original 2010–2025 and 2010–2026 full search exports (both databases) are committed at `living-search/data/historical_imports/`. Run once:

**Actions → "ImmunoRad import historical vintages (2025 + 2026)" → Run workflow**

This merges all four files into `archive.json` (deduped by PMID/DOI/title, tagged `prior_search: true` — they won't show up in "To Screen"), and adds their **raw** record counts to `state.json`'s `total_records_identified` counter, which is what powers "Records identified" and "Duplicates removed" on the PRISMA tab.

**Order matters a bit:** run this historical import once, early — ideally before or alongside the other one-off imports (Embase CSV import, PubMed backfill), so the cumulative counter reflects the full picture. Every script that touches PubMed/Embase data now updates this same counter, so future weekly runs and any new Embase imports keep it accurate automatically.

**Known small gap:** the two PubMed workflow runs from before this counter existed (the very first manual run and the first scheduled Monday run) aren't retroactively included, since incremental runs advance their search window forward and can't be safely re-run for the same historical dates. If you want that exact-to-the-record precision, check those two runs' logs (Actions → click the run → "esearch returned N candidate PMIDs") and let me know the numbers — I'll fold them into `state.json` by hand once.

### `data/manual_prisma_stats.json` format

```json
{
  "full_text_excluded_total": 218,
  "full_text_exclusion_reasons": [
    { "label": "Not most recent report", "n": 79 }
  ],
  "included_trials": 41,
  "included_patients": 15049,
  "subgroups": [
    { "label": "ICI vs no ICI", "trials": 35, "outcomes": ["OS", "PFS", "EFS"] }
  ]
}
```
It's currently seeded with the **last completed round's** (published manuscript) numbers as a placeholder — update it once you've done full-text review on this round's new "Eligible for full text" candidates.

## What this does and doesn't guarantee

- The PubMed query is copied **verbatim** from `PubMed_Search_prompt.docx` — not simplified.
- Weekly runs detect new records by **date added to PubMed** (`edat`); the one-off backfill uses **publication date** (`pdat`) to also catch late-indexed older records.
- Deduplication runs on PMID, DOI, and normalized title.
- Include/Exclude is a screening aid with a shared audit trail — it does not itself alter the manuscript or the pooled analysis. That step is still yours.
