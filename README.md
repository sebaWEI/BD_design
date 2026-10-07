# BD_design

Energy scoring, polymer–kinetics landscapes, and retrospective BD design for
**iGEM PekingHSC 2026** (HEPHA-RNA).

BLAST / variant **safety screening** lives in
[`BD_screening`](https://github.com/sebaWEI/BD_screening) (`bsst`). This
repo takes candidate sites (or wet-lab tiles) and scores them.

## Scope

- ViennaRNA **RNAup** site scoring (`scripts/score_rnaup.py`)
- Closed-loop / constrained-ring positional scores and ED capture–release
  (`scripts/mRNA_model*.pdf`, `plot_hepha_landscape.py`,
  `try_closedloop_occupancy.py`)
- Engineering-cycle figures and UTR secondary-structure plots
- LETM1 / NSD2 wet-lab retrospective tables under `data/retrospective_study/`
- BD primer / construct dossiers under `examples/`

## Layout

Meant to sit **inside** a local `BD_screening` clone as
`BD_screening/BD_design/`, so scripts can read sibling `runs/` and
`examples/*.fasta`:

```bash
cd BD_screening
git clone https://github.com/sebaWEI/BD_design.git
```

## RNAup after screening

```bash
# score sites that already passed BLAST (+ optional variants)
../.venv/bin/python scripts/score_rnaup.py \
  --utr ../examples/LETM1.fasta \
  --table ../runs/<id>/candidates.tsv \
  --out data/letm1_rnaup.tsv

# score every wet-lab tile (ignore BLAST for retrospective energy tables)
../.venv/bin/python scripts/score_rnaup.py \
  --utr ../examples/LETM1.fasta \
  --sites ../examples/LETM1.sites.fasta \
  --out data/retrospective_study/letm1_tiles_rnaup.tsv
```

Needs `RNAup` on `PATH` (ViennaRNA).

## Other figures

```bash
../.venv/bin/python scripts/plot_engineering_figures.py
../.venv/bin/python scripts/plot_hepha_landscape.py
../.venv/bin/python scripts/try_closedloop_occupancy.py
../.venv/bin/python scripts/plot_utr_secondary_structure.py
```

Outputs land in `docs/figures/engineering/`.
