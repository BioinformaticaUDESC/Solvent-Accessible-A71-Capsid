# Per-residue measurement tables

One file per complex: `<PDB>_<chain>_pegada_V49.tsv`, 45 files for 45
complexes, 13499 data rows in total.

The chain in the file name is the reference chain the analysis measured.
Where a deposit labels the same protein under more than one chain
identifier, only the chain used by the analysis is included here.

## Columns

| Column | Meaning |
|---|---|
| `residue` | author-numbered residue of the reference chain |
| `resname` | three-letter residue name |
| `radius_A` | radius of the Ca from the particle centroid, in angstrom |
| `shell_pct` | percentile rank of that radius among the Ca radii of all capsid chains |
| `blocking` | outward-cone atomic occupancy: atoms of the assembled particle inside a 45 degree cone of 30 angstrom radius along the outward radial axis |
| `bfactor` | mean B-factor of the residue as deposited |
| `rank` | rank of the residue by occupancy within the complex |
| `in_footprint` | whether the receptor of that complex contacts the residue |
| `passes` | decision under the adopted parameterization |
| `passes_V49` | decision under the alternative parameterization the manuscript reconciles |

`passes` and `passes_V49` differ for 2640 rows, 19.56% of the total. That
fraction is a property of the two parameterizations, not of any single
complex.

## Counts to expect

Anyone recomputing from this deposit should find 45 files, 45 complexes
and 13499 rows. A different total means a different set of files was summed:
working copies of this project have circulated with duplicate tables under
`_V49_02` and `_V49_03` suffixes and with a second chain label for one
deposit, and none of those are part of this deposit.

The descriptor tables used for the comparative benchmark, with every
alternative descriptor computed on the assembled particle, are in
`../benchmark/`, which also maps each published table and figure to the
file that produces it.
