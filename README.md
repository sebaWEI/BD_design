# BD_design

Energy scoring, polymer–kinetics landscapes, and retrospective BD design for
**iGEM PekingHSC 2026** (HEPHA-RNA).

BLAST / variant **safety screening** lives in
[`BD_screening`](https://github.com/sebaWEI/BD_screening) (`bsst`). This
repo takes candidate sites (or wet-lab tiles) and scores them.

## Scope

- ViennaRNA **RNAup** site scoring (`model/score_rnaup.py`)
- Closed-loop / constrained-ring positional scores and ED capture–release
  (`model/mRNA_model*.pdf`, `model/main.tex`,
  `simulation/plot_hepha_landscape.py`,
  `simulation/try_closedloop_occupancy.py`)
- Engineering-cycle figures and UTR secondary-structure plots
  (`simulation/`)
- LETM1 / NSD2 wet-lab retrospective tables under `data/retrospective_study/`
- BD primer / construct dossiers under `data/constructs/`

## Layout

```
data/            # constructs, retrospective tables
docs/            # engineering notes + figures
model/           # RNAup scorer + theory (main.tex, mRNA_model*.pdf)
simulation/      # closed-loop occupancy, landscapes, score builders, plots
test/            # pytest
```

Meant to sit **inside** a local `BD_screening` clone as
`BD_screening/BD_design/`, so modules can read sibling `runs/` and
`examples/*.fasta`:

```bash
cd BD_screening
git clone https://github.com/sebaWEI/BD_design.git
```

Run commands from the `BD_design` root (parent `.venv` shared with screening).

## RNAup after screening

```bash
# score sites that already passed BLAST (+ optional variants)
../.venv/bin/python -m model.score_rnaup \
  --utr ../examples/LETM1.fasta \
  --table ../runs/<id>/candidates.tsv \
  --out data/retrospective_study/letm1_rnaup.tsv

# score every wet-lab tile (ignore BLAST for retrospective energy tables)
../.venv/bin/python -m model.score_rnaup \
  --utr ../examples/LETM1.fasta \
  --sites ../examples/LETM1.sites.fasta \
  --out data/retrospective_study/letm1_tiles_rnaup.tsv
```

Needs `RNAup` on `PATH` (ViennaRNA).

## Simulation / figures

```bash
../.venv/bin/python -m simulation.build_retrospective_scores
../.venv/bin/python -m simulation.try_closedloop_occupancy
../.venv/bin/python -m simulation.plot_hepha_landscape
../.venv/bin/python -m simulation.plot_engineering_figures
../.venv/bin/python -m simulation.plot_utr_secondary_structure
```

Outputs land in `docs/figures/engineering/`.

## Tests

```bash
../.venv/bin/python -m pytest test/
```
