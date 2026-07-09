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

**What I already did on my end (no network access needed for this part):**
- Imported your uploaded `Pubmed_2026_NEW_vs_2025_RAYYAN.csv` (236 records) into `data/archive.json`, tagged `"prior_search": true`, so none of those will ever reappear as "new" — they're the known baseline.

**What still needs to run for real** (I can't call the PubMed API directly from this chat — see below):
1. Repo → **Actions** tab → **"ImmunoRad one-time backfill (since 2025-01-01)"** → **Run workflow**
2. This runs the exact literal query with `mindate=2025-01-01`, `datetype=pdat` (publication date, to also catch late-indexed older records), dedupes against the full archive (including your 236 CSV records), and adds only genuinely new hits to `new_hits.json`
3. It's a one-off — it does **not** affect the normal Monday schedule or `state.json`'s incremental tracking

## Why I couldn't just run the search myself in this chat

My tools here can't reach the PubMed API directly (network restrictions in this environment) and can't call arbitrary URLs I haven't been given or found via search. GitHub Actions runs on GitHub's own servers, which have normal internet access — so that's where the actual API calls happen. Everything I *could* do locally (parsing your CSV, writing the scripts, updating the archive with your prior results) is already done.

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

## What this does and doesn't guarantee

- The PubMed query is copied **verbatim** from `PubMed_Search_prompt.docx` — not simplified.
- Weekly runs detect new records by **date added to PubMed** (`edat`); the one-off backfill uses **publication date** (`pdat`) to also catch late-indexed older records.
- Deduplication runs on PMID, DOI, and normalized title.
- Include/Exclude is a screening aid with a shared audit trail — it does not itself alter the manuscript or the pooled analysis. That step is still yours.
<!-- redeploy trigger -->
