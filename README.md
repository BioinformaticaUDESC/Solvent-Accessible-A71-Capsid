# Solvent-accessible yet externally occluded

Existing data and analysis scripts for the Enterovirus A71 capsid study.
This archive preserves selected existing files; its preparation does not run
simulations, recompute measurements, or validate scientific outputs.

## Descriptor

Outward-cone atomic occupancy counts particle atoms in an external cone with
a 45-degree half-angle, its vertex at the residue C-alpha atom, its axis along
the outward radial vector from the particle centroid, and a search up to 30 Å.

## Layout

| Path | Contents |
|---|---|
| data/per_residue/ | Existing per-residue footprint tables |
| data/radial_profiles/ | Existing radial profiles |
| data/ | Existing validation and 3VBS SASA-versus-radial tables |
| scripts/ | Original project and analysis scripts |
| figures/ | Requested individual panels and composite figures |
| manuscript/ | Requested manuscript, cover letter and methods supplement |
| MANIFEST.sha256 | SHA-256 of every repository file except the manifest itself |
| LICENSE | License scope and applicable license files |

Copied files by top-level folder: {"data": 53, "figures": 23, "manuscript": 1, "scripts": 18}.
Unresolved missing items or conflicting selections: 6.
See the separate preparation report before publication; it is outside the ZIP.

## Reproducing figures

The supplied figure scripts are V58_gera_figuras.py, V59_conserta_figuras.py,
V60_fecha_figuras.py, V62_cone_e_fundo.py, V64_acabamento.py and
V65_monta_paineis.py, when present. Inspect the actual scripts, required inputs,
output locations and dependencies before executing them. Local paths may need
adaptation. The preparation utility has not executed or independently verified
these scripts. V63_sasa_3vbs.py is the named source of the SASA table; a missing
table is reported rather than regenerated.

Before reproduction, obtain the necessary coordinates directly from the PDB,
adapt the source paths, and run the figure scripts in a separate working copy.
Successful end-to-end reproduction is not asserted by this packaging step.

## Excluded exploratory docking

The exploratory docking layer is **not reproducible** and is not included.
Its original execution tree is unavailable, and surviving summaries are
insufficient to reconstruct it. No docking summaries or reconstructed runs
are included by this preparation utility.

## Structures

Coordinates originate from the Protein Data Bank and are not redistributed
in this repository. The requested deposits are 3VBF, 3VBH and 3VBS
(Wang et al., Nature Structural & Molecular Biology, 2012,
https://doi.org/10.1038/nsmb.2255) and 6I2K
(Zhou et al., Nature Microbiology, 2019,
https://doi.org/10.1038/s41564-018-0319-z).

## Licenses and citation

Code in scripts/: MIT (LICENSE-CODE.txt).
Data in data/: CC BY 4.0 (LICENSE-DATA.txt).
No additional license is assigned here to figures/ or manuscript/.
Existing third-party notices retain their effect. See LICENSE and CITATION.cff.
No release or DOI is claimed by this preparation step.
