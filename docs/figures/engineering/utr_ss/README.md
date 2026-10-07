# UTR secondary structures (RNAfold / RNAplot)

MFE structures, 37 °C, ViennaRNA 2.7. Coordinates on 3′UTR match `examples/*.fasta`.

| File | What |
|------|------|
| `LETM1_5utr_ss` / `NSD2_5utr_ss` | Full 5′UTR |
| `LETM1_3utr_ss` / `NSD2_3utr_ss` | Full 3′UTR (crowded; use zooms) |
| `*_3l_zoom_ss` / `*_3s_zoom_ss` | ±80 nt around BD ladders |

Marks (EPS `omark`):

- **blue** — 5′ end / cap (first 8 nt of 5′UTR)
- **green** — AUG-adjacent (last 8 nt of 5′UTR) and stop-adjacent (first 8 nt of 3′UTR)
- **orange** — polyA-proximal (last 15 nt of 3′UTR)
- **red** — 3l BD tiles (name at tile midpoint)
- **purple** — 3s BD tiles (name at tile midpoint)
