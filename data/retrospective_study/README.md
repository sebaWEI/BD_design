# Retrospective wet-lab + RNAup scores

| File | Contents |
|------|----------|
| `260722-24、31 LETM1fluc.xlsx` | Raw LETM1 luciferase plates |
| `260814—NSD2fluc.xlsx` | Raw NSD2 luciferase plates |
| `LETM1-3l.gbk` / `NSD2-3s+3l.gbk` | SnapGene constructs |
| `retrospective_fluc_replicates_20ng.tsv` | Per-replicate fold (20 ng/well blocks) |
| `retrospective_scores.tsv` | Per-site mean fold + RNAup energies (input to Engineering Cycle 3 figures) |

Build / refresh::

```bash
# Score every wet-lab tile (raise BLAST length so RNAup sees all sites)
uv run bsst filter --gene LETM1 --utr examples/LETM1.fasta \
  --sites examples/LETM1.sites.fasta --offtarget-min-length 100 --rnaup
uv run bsst filter --gene NSD2 --utr examples/NSD2.fasta \
  --sites examples/NSD2.sites.fasta --offtarget-min-length 100 --rnaup

uv run python scripts/build_retrospective_scores.py \
  --letm1-run runs/<LETM1_run_id> --nsd2-run runs/<NSD2_run_id>
uv run python scripts/plot_engineering_figures.py
```

Figures land in `docs/figures/engineering/`.
