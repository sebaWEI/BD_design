# BD_design

Energy scoring, polymer–kinetics landscapes, and BD design for
**iGEM PekingHSC 2026** (HEPHA-RNA).

- [Team model repo](https://github.com/sebaWEI/PekingHSC-2026-Model)
- [iGEM team](https://teams.igem.org/6371)
- [Wiki](https://2026.igem.wiki/pekinghsc/)


## Scope

- ViennaRNA **RNAup** site scoring (`model/energy/`)
- Polymer positional models from `model/main.tex` (`model/main/`: FJC,
  Gaussian ring, constrained ring)
- Closed-loop / ED landscapes in `simulation/`
- Clean wet-lab tables in `data/` (cis luciferase + 1018 LETM1 BD scan)

## Layout

```
data/            # cleaned CSVs only (see data/README.md)
docs/figures/    # wiki figures (cis + BD-scan fits)
model/
  energy/        # RNAup scoring
  main/          # FJC → GaussianRing → ConstrainedRing (main.tex)
  main.tex       # theory write-up
simulation/      # fits + wiki plots
test/            # pytest
```

Meant to sit **inside** a local `BD_screening` clone as
`BD_screening/BD_design/`, so modules can read sibling `runs/` and
`examples/*.fasta`.

Run commands from the `BD_design` root (parent `.venv` shared with screening).

## Cis / BD-scan fits

```bash
../.venv/bin/python -m simulation.fit_cis_sites
../.venv/bin/python -m simulation.plot_cis_sites_wiki
../.venv/bin/python -m simulation.fit_letm1_bd_scan
```

Figures → `docs/figures/model/cis_fit_wiki.*` and `letm1_bd_scan_fit.*`.

## Polymer classes

```python
from model.main import FreelyJointedChain, GaussianRing, ConstrainedRing

fjc = FreelyJointedChain(n_segments=200)
c_loop = fjc.effective_concentration()

ring = GaussianRing(n_total=120)
c_ring = ring.effective_concentration(n=40)

constrained = ConstrainedRing.from_params(n_total=120, n_5a=20, k_match=0.02)
c_star = constrained.effective_concentration(n_b3=20)
```

## Tests

```bash
../.venv/bin/python -m pytest test/
```
