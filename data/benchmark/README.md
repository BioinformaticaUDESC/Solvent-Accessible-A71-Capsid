# Benchmark records

Every number in Table 4, Figure 5 and Supplementary Tables S5, S6 and S7 of
the manuscript is produced by a file in this directory. Nothing here was
transcribed by hand.

## Where each published number comes from

| Published object | File |
|---|---|
| Table 4, all columns | `results/V94_tabela4.tsv` |
| Table 4, inclusion accounting | `results/V94_inclusao.tsv` |
| Figure 5, points and bars | `results/V94_tabela4.tsv` (columns `dif_pareada`, `dif_lo`, `dif_hi`) |
| Supplementary Table S5 | `results/V98_tabela4_26.tsv` |
| Supplementary Table S6 | `results/V94_grade3x3.tsv` |
| Supplementary Table S7 | `results/V101_grade_operacao.tsv` |
| Every descriptor, per residue | `per_residue/<PDB>_descritores.tsv` |
| Residue depth, per residue | `per_residue/<PDB>_profundidade.tsv` |

## Columns of the per-residue tables

`residuo` is the author-numbered residue of the reference chain.
`blocking_pub` is the occupancy as deposited in the earlier per-residue
files; `cone_45_30` is the same quantity recomputed from the reconstructed
assembly for this comparison. `cone_<half-angle>_<radius>` are the nine
grid cells. `cone_45_30_sem_vizinhos` excludes the vertex residue and its
two sequential neighbours. `sasa_1.4`, `sasa_5.0` and `sasa_10.0` are
solvent-accessible areas on the assembled particle at those probe radii.
`hse_fora` is half-sphere exposure along the Ca-to-Cb vector. `shell_pct`
is the radial percentile as published. `contato` is the observed contact
label. `prof_media` and `prof_CA` are residue depth, mean over the atoms of
the residue and at the Ca.

## How the numbers were produced

Assemblies were rebuilt from the deposited symmetry operators and verified
against the published per-residue counts before any descriptor was
computed: a complex entered the comparison only when the recomputed
occupancy differed by at most one atom at every residue and no residue
crossed the 106-atom threshold. `results/V94_inclusao.tsv` records that
check for every complex, including the four that did not pass.

Solvent-accessible area on the assembled particle is exact rather than
truncated: an atom can influence another atom's area only when the
probe-expanded spheres intersect, at most r_i + r_j + 2*probe = 23.8 A for
the 10 A probe, so the calculation was restricted to the reference chain
plus every particle atom within 25 A of it.

Residue depth is the mean distance from the atoms of a residue to the
nearest point of bulk solvent connected to the exterior of the particle,
on a 1 A grid, with the flood fill seeded only outside the particle so
that the 25 A truncation cannot create solvent inside the capsid.

## Scripts

The scripts in `code/benchmark/` run in the order of their numbers and
read only the files named above. Seeds are stated in each script: 20261116
for the calibration/held-out permutation (Python's `random.Random`, applied
to the validation-table order with the two assembly-requirement exclusions
dropped first, second half held out) and 20261125 for the 5000-resample
intervals.
