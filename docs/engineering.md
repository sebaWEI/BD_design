# Engineering cycles (wiki)

This page is the engineering narrative for the iGEM wiki. It matches the
default behaviour of the current `bsst` CLI.

## Cycle 1 — Where can it bind?

**Design.** The 3′UTR of a target transcript can span thousands of
nucleotides. A binding site should support target recognition and minimise
off-target interactions. Sequence similarity via BLAST was the first filter.

**Build.** The first Binding Site Selection Tool (BSST) searched candidates
against reference transcripts and removed sites with potentially problematic
matches.

**Test.** Alignments were produced, but not every match is a relevant
off-target. A CDS hit on another transcript is not the interaction we set out
to avoid for a 3′UTR-acting element. Removing every BLAST hit was too
conservative.

**Learn.** Similarity alone is not enough. Later versions needed transcript
annotation so that 5′UTR and 3′UTR matches could be separated from other
regions.

## Cycle 2 — Will the site remain stable?

**Design.** Genetic variation inside a selected site may change binding across
individuals, so the tool should also ask whether a site is genetically stable.

**Build.** Population variation from dbSNP `common_all_20180418` was added as
an optional module. When that filter is enabled, sites overlapping annotated
variable positions can be removed. BLAST was revised to use transcript
annotation and to treat only 5′UTR and 3′UTR hits as relevant. Same-strand
alignments alone are kept: a reverse-strand hit is the antisense element
sequence on another transcript, not a site the element can bind.

**Test.** The tool then combined off-target context with optional genetic
robustness. Wet-lab users still found the raw BLAST table hard to read.

**Learn.** Correctness is not enough without interpretability. Later versions
summarise gene, hit length, identity, and transcript region (5′UTR / 3′UTR).
The UTR-hit length threshold is user-adjustable (`--offtarget-min-length`).

## Cycle 3 — How well will it work?

**Design.** After the safety filter, we asked whether BSST could rank sites by
expected translational upregulation using RNA thermodynamics.

**Build.** RNAup estimated opening and interaction free energy; candidates were
ranked by total ΔG.

**Test.** Against wet-lab upregulation for LETM1 and NSD2, the ranking did not
match experiment. Alternative reweightings of the energy terms also failed to
give a shared, reproducible rule.

**Learn.** Equilibrium binding free energy alone does not describe functional
activity. A favourable interaction is not a favourable translational outcome.

## Cycle 3.5 — What can we reliably claim?

**Design.** A more dynamical model would add parameters we cannot constrain.
For a therapeutic design pipeline, an unstable prediction is worse than no
prediction. We therefore limited the tool to safety-oriented filtering.

**Build.** RNAup remains an optional exploratory module (`--rnaup`), off by
default. The default core is same-strand BLAST against the pinned GENCODE 45
transcriptome, with only 5′UTR / 3′UTR hits above an adjustable length
threshold dropping a site. Variation filtering (`--variants`) stays optional
and off by default. Dropped sites name gene, transcript, UTR, identity, and
both intervals.

**Test.** The workflow is a computational pre-screen before wet-lab validation.
The thermodynamic score is not presented as a validated activity predictor.

**Learn.** A useful design tool should not predict what its model cannot
reliably explain. BSST is a filter and decision-support tool, not a black-box
predictor of translational upregulation.
