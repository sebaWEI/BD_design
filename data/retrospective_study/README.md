# Retrospective wet-lab + RNAup scores

| File | Contents |
|------|----------|
| `retrospective_scores.tsv` | Per-site mean fold + RNAup energies (Engineering Cycle 3 figures) |
| `retrospective_fluc_replicates_20ng.tsv` | Per-replicate luciferase folds |
| `*fluc.xlsx` / `分子动力学结果.xlsx` | Raw wet-lab / MD sheets |
| `*.cdna.fa` / `*.gbk` | Transcript / construct references |

Rebuild energy columns with `scripts/score_rnaup.py` (not `bsst filter`):

```bash
# From BD_design/, with RNAup on PATH and parent screening checkout
../.venv/bin/python scripts/score_rnaup.py \
  --utr ../examples/LETM1.fasta \
  --sites ../examples/LETM1.sites.fasta \
  --out data/retrospective_study/letm1_tiles_rnaup.tsv

../.venv/bin/python scripts/score_rnaup.py \
  --utr ../examples/NSD2.FASTA \
  --sites ../examples/NSD2.sites.fasta \
  --out data/retrospective_study/nsd2_tiles_rnaup.tsv

../.venv/bin/python scripts/build_retrospective_scores.py
```

Then:

```bash
../.venv/bin/python scripts/plot_engineering_figures.py
```
