# Clean tables

| File | Content |
|------|---------|
| `cis_luciferase.csv` | Cis UTR luciferase replicates (Moderna / TMSB10) |
| `cis_fit_params.csv` | Constrained-ring fit parameters |
| `cis_fit_predictions.csv` | Per-site constrained-ring predictions |
| `1018_bd_scan.csv` | 1018 LETM1 BD scan (RLU / folds / sequences / MMPBSA when available) |
| `fit_letm1_bd_scan_params.csv` | Constrained-ring params for the BD scan |
| `fit_letm1_bd_scan_predictions.csv` | Per-site BD-scan predictions |

`mmpbsa_kcal_mol` in `1018_bd_scan.csv` is filled for plate samples `1`–`9` (LETM1-3l-1…9). New tiles (`-1`…`-10`) are left blank until those runs finish.
