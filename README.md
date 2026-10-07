# BD_design

Energy, polymer–kinetics scoring, and retrospective BD design for
**iGEM PekingHSC 2026** (HEPHA-RNA).

This repository is **not** the BLAST / variant filter. Screening lives in
[`BD_screening`](https://github.com/sebaWEI/BD_screening) (`bsst`).

## Scope

- Closed-loop / constrained-ring positional scores and ED capture–release
  (`scripts/mRNA_model*.pdf`, `scripts/plot_hepha_landscape.py`,
  `scripts/try_closedloop_occupancy.py`)
- Engineering-cycle figures and UTR secondary-structure plots
- LETM1 / NSD2 wet-lab retrospective tables under `data/retrospective_study/`
- BD primer / construct dossiers under `examples/`

## Layout

This checkout is meant to sit **inside** a local `BD_screening` clone as
`BD_screening/BD_design/`, so scripts can read sibling `runs/` and
`examples/*.fasta` for site coordinates. Clone either way you prefer:

```bash
# nested (scripts that need screening runs work out of the box)
cd BD_screening
git clone https://github.com/sebaWEI/BD_design.git

# standalone
git clone https://github.com/sebaWEI/BD_design.git
# then pass paths / copy runs as needed
```

## Reproduce figures

Use the parent screening `.venv` (or any env with pandas / matplotlib / scipy):

```bash
../.venv/bin/python scripts/plot_engineering_figures.py
../.venv/bin/python scripts/plot_hepha_landscape.py
../.venv/bin/python scripts/try_closedloop_occupancy.py
../.venv/bin/python scripts/plot_utr_secondary_structure.py
```

Outputs land in `docs/figures/engineering/`.
