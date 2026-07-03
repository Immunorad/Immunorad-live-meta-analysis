# Eerst dit lezen: hoe dit zich verhoudt tot je bestaande website

Je bestaande living-meta-analyse website blijft **volledig intact**. Niets wordt overschreven.

## Waarom dit veilig is

Alles zit nu in een eigen map: `living-search/`. Geen van de nieuwe bestanden heet zoals iets wat een normale website al zou hebben (`index.html`, `style.css`, etc.), dus er is geen overlap mogelijk — behalve de GitHub Action, die moet wel op een vaste plek staan (`.github/workflows/`), zie hieronder.

```
jouw-repo/
├── index.html                        ← jouw bestaande site, blijft ongewijzigd
├── (jouw overige bestaande bestanden) ← blijven ongewijzigd
├── .github/
│   └── workflows/
│       └── immunorad_living_search.yml   ← NIEUW
└── living-search/                        ← NIEUW, alles hierin
    ├── immunorad_times.html
    ├── scripts/
    │   ├── pubmed_search.py
    │   └── merge_embase.py
    └── data/
        ├── state.json
        ├── archive.json
        └── new_hits.json
```

De krant wordt straks bereikbaar op een eigen sub-URL, bijvoorbeeld:
`https://immunorad.github.io/immunorad-live-meta-analysis/living-search/immunorad_times.html`

Je bestaande homepage blijft gewoon op de hoofd-URL staan.

## Twee dingen die je even moet checken vóór het uploaden

1. **Heb je al een map genaamd `living-search` of `.github/workflows/` met een bestaand bestand?**
   Zo goed als zeker niet, maar check even in je repo of die namen al in gebruik zijn. Zo ja, laat het weten en ik hernoem alles.

2. **Heb je GitHub Pages al aanstaan voor deze repo?**
   Ga naar **Settings → Pages** en kijk wat er bij "Source" staat:
   - Staat het al op **Deploy from a branch → main → / (root)**? Dan hoef je niets te veranderen — de krant verschijnt vanzelf op de sub-URL zodra je de bestanden hebt geüpload.
   - Staat het op iets anders (bv. een `/docs`-map, een aparte `gh-pages`-branch, of een custom domain)? Laat het me weten voordat je iets aanpast in de instellingen — dan pas ik de mapstructuur aan zodat 'ie op de juiste plek terechtkomt, in plaats van dat je de Pages-instelling van je bestaande site aanraakt.

## Stappen (als beide bovenstaande punten oké zijn)

1. Ga naar je repo → **Add file → Upload files**
2. Sleep de hele map `living-search/` erin, én de map `.github/` (of alleen het bestand `.github/workflows/immunorad_living_search.yml` als je al een `.github`-map hebt — dan moet je 'm daarin plaatsen, niet vervangen)
3. Commit changes
4. Tabblad **Actions** → "ImmunoRad living search (weekly PubMed)" → **Run workflow** om de eerste keer handmatig te starten
5. Open de krant op `.../living-search/immunorad_times.html`

Zie `README_living_search.md` (in de `living-search`-map) voor de volledige uitleg over de EMBASE-koppeling en hoe de krant werkt.
