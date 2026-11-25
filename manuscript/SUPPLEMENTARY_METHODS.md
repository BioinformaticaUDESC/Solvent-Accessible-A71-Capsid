# Supplementary Methods — directional occupancy of capsid surfaces

*Companion to "Solvent-Accessible Yet Externally Occluded: Directional Analysis of
a Docking-Selected Site in the Mature Enterovirus A71 Capsid". Two kinds of
material are collected here and they carry different weight. Sections S2 to S7 and
S9 document the exploratory layer that generated the hypothesis and procedures not
used in the central inference: that campaign is provenance for the hypothesis and
is not offered as evidence for it. Section S8 is different in kind — it is a
control on the implementation of the main article's method, and it supports the
assembly step that the central analysis depends on. Values marked `[CHECK]`,
`[REQUIRED]` or `[CONFIRM]` are not verified and must be resolved or removed
before submission.*

---

## S1. Scope of this supplement

The main article reports one contribution: the definition, calibration and
evaluation of outward-cone atomic occupancy as a criterion for excluding
docking-selected sites that are solvent-accessible but externally occluded in the
assembled particle. The exploratory campaign that produced the candidate site is
reported here as provenance.

Three statements bound everything below. The intermediate run tree of the
exploratory docking layer is no longer accessible: its preserved output tables can
be **re-analysed**, but the original docking execution is **not reproducible** from
a table of results, and the two must not be conflated. The docking benchmark has
complete result files for six complexes and no aggregate table is presented here;
moving it to this supplement does not make it resolved. And no unverified value is
presented as confirmed.

---

## S2. Protein panels

The host panel comprises thirteen human proteins drawn from three functional axes
of nigrostriatal vulnerability: tyrosine hydroxylase, DOPA decarboxylase,
monoamine oxidase B, catechol-O-methyltransferase, the dopamine transporter
DAT/SLC6A3 and the vesicular monoamine transporter VMAT2/SLC18A2 for dopaminergic
homeostasis; beclin-1, MAP1LC3B and glucocerebrosidase/GBA1 for proteostasis and
the autophagy–lysosome pathway; and PINK1, parkin/PRKN, PARK7/DJ-1 and VDAC1 for
mitochondrial quality control. The viral panel comprises thirteen proteins from
six viruses: NS1, NP, M1 and the PB2 cap-binding domain of influenza A; 3Cpro,
3Dpol and VP1 of enterovirus A71; UL30 and UL42 of herpes simplex virus 1; the 3C
protease of coxsackievirus B3 and of B4; and the NS3 protease and NS5 polymerase
of Japanese encephalitis virus.

Two proteins planned for the panel were segregated on structural grounds and are
absent from the matrix: the full-length LRRK2 model, whose audit found 40,136
readable atoms together with 2,615 lines carrying coordinates incompatible with
fixed-field format and spatial extents inconsistent with robust parsing, and the
PB1 domain of p62/SQSTM1. Both are recorded as technical limitations rather than
negative biological results. The exclusion of LRRK2 removes the genetic-risk axis
from the matrix.

## S3. Structural preparation of the exploratory panel

The matrix reported here was computed on a preparation procedure reconstructed to
be re-executable, because the original preparation could not be reproduced from
the surviving records. The reconstructed procedure starts from the deposited
mmCIF file; preserves the original residue numbering, on which the identity of
site 22 and of the 214/216/217 footprint depend; converts selenomethionine to
methionine before heteroatom removal; strips alternative conformations, hydrogens,
ligands and waters; and fails explicitly when a requested chain is absent instead
of proceeding with whatever it finds. Each prepared structure is registered in a
manifest carrying the deposit identifier, the chains retained, the residue count
and checksum prefixes of both the source and the output file. Twenty-eight
structures of the protein panel were prepared, of which twenty-six enter the
matrix. These counts refer to the panel and are unrelated to the twenty-eight
capsid–receptor complexes and twenty-six-complex split described below, which
have their own preparation and manifest.

The two preparations are not comparable with each other. In several cases they act
on different molecular objects, with different chain composition and residue
content, so a score difference between them reflects a difference of input and not
instability of the procedure. The coordinate files of the earlier preparation are
not available; no claim is made here about what those coordinates contained, and
no value from the earlier preparation is combined with any value from the
reconstructed one.

## S4. Host × viral matrix and decomposition of effects

The matrix was recomputed in full on the reconstructed preparation, one chain per
structure and with frozen parameters, so that every score is traceable from the
deposit code to the final value. **`[REQUIRED: name the docking engine, its`**
`version and the scoring function that produced these 169 scores, with the`
`parameter file. The Software section lists the tools used across the project and`
`does not establish which engine generated this matrix.]` Its dimensions are thirteen host proteins by
thirteen viral proteins, giving 169 pairs. Score variation was decomposed
additively into an effect attributable to the host partner, an effect attributable
to the viral partner and a combination term, with effects obtained as deviations
from the marginal means. The design contains no within-cell replication, so the
combination term and the run-to-run error are not separable from the matrix alone
and the combination term is reported as a ceiling on the pair effect rather than
as an estimate of it.

## S5. Seed replicates of the matrix

To estimate the score variation of a single cell across runs, a subset of pairs
was re-executed with additional seeds and all other parameters held. Six pairs
were re-run with two additional seeds each, with a median deviation of
3.832 units across the three runs of a cell. `[CHECK: the project log recording
this value notes that an earlier text reported 3.83; these are the same
measurement at different precision and one of them should be carried
consistently.]` These pairs were chosen
after the matrix had been computed, so as to span the observed range of
residuals; the choice therefore depended on the outcome the replicates assess and
does not constitute prior registration. Selection by extreme value can bias the
estimate of run-to-run variation, and the direction and magnitude of that bias
were not quantified. With few runs per cell the precision of each deviation is low
and was not estimated.

## S6. Docking benchmark against known complexes

Before any ordering of the matrix was interpreted, the discriminative capacity of
the procedure was measured on protein–protein complexes whose bound structure is
known. Reference complexes were submitted to the protocol used for the matrix —
same sampling, same number of agents per swarm, same number of steps — and each
was run with three distinct scoring functions and five seeds per function. For
each run we recorded the root-mean-square deviation of the best sampled solution
from the deposited structure, the score assigned to that solution, its rank in
the score ordering, and the correlation between score and deviation over the
entire final population. A procedure able to recognise the correct
interface would place the lowest-deviation solution among the best-scored, and
would show a monotonic association between score and deviation. **`[REQUIRED: the`**
`sign of that association follows from the orientation of the score and must be`
`read from the score files rather than assumed. If the reported score is an`
`energy, so that lower is better, a good function gives a positive`
`score–deviation correlation; if the score is oriented so that higher is better,`
`the correlation is negative. The orientation used in these files has not been`
`verified here and no sign is asserted.]`

Six rigid-body enzyme–inhibitor complexes have complete result files: 1PPE, 2SIC,
7CEI, 1ACB, 1AVX and 1CGI. The aggregate table of this benchmark reports three of
them, with no recorded selection criterion, and case difficulty does not account
for the split since all six are rigid-body cases of the same class.

**`[REQUIRED before submission: the benchmark is not reported in this draft.`**
`Either the per-complex table for all six is produced — PDB code, scoring`
`function, seed, completion status, deviation of the best-scored solution, its`
`rank, and the score–deviation correlation over the final population — or the`
`claim is reduced to the three complexes whose results are tabulated, with the`
`omission of the other three stated. Still to be specified from the run logs: the`
`identity of the three scoring functions, whether the deviation is interface or`
`ligand RMSD and over which atom selection, the LightDock version, and how many`
`runs per complex exist. The planned design of three functions by five seeds`
`implies fifteen runs per complex, forty-five across three complexes or ninety`
`across six; which of these were executed and analysed must come from the logs`
`and is not to be inferred from the design.]`

## S7. Bilayer delimitation in membrane proteins

Composition descriptors computed over the total surface of a membrane protein
include the face directed at the bilayer, which is not available to a soluble
partner. For the five membrane-associated proteins of the panel — three delimited
as bilayer-spanning and two as monotopic, with a single plane — the region
delimited by the bilayer was taken from a membrane-orientation database and the
descriptors were recomputed over the surface lying outside that region only.

## S8. Verification of assembly expansion against the repository

The expansion step was checked against the repository rather than against earlier
output of our own. For three deposits spanning both operator conventions — 5NEU
and 3J9F, whose operators are pure rotations about the coordinate origin, and
6I2K, whose symmetry centre lies 468 Å from the origin with translations up to
935 Å — the biological assembly file served by the repository was downloaded and
measured independently.

| deposit | official assembly | our expansion |
|---|---|---|
| 5NEU | 360 chains, 103,200 Cα, 546.4–666.2 Å | 360 chains, 103,200 Cα, 546.4–666.2 Å |
| 3J9F | 780 chains, 99.5–244.7 Å | 780 chains, 99.5–244.7 Å |
| 6I2K | 660 chains, 103.4–213.5 Å | 660 chains, 103.4–213.5 Å |

Chain count, Cα count and radial extent agree to the first decimal in all three,
so the expansion step is not a source of the discrepancy reported for 5NEU in the
main article. Both figures derive from the same operator records, so the agreement
demonstrates consistent reading and not independent physical validation.

---

## S9. Exploratory docking software

The exploratory screening used LightDock/Python with active restraints
corresponding to the information-driven mode, and DFIRE for relative scoring
[S1–S4]. Supporting
analyses not used in the central inference employed IQ-TREE 2 and ModelFinder for
phylogeny and ancestral reconstruction, MAFFT v7.526 for alignment, HADDOCK3
v2026.7.0 for a reassessment, Cytoscape for networks and GROMACS 2024.2 for
topology generation and restrained minimisation [S5–S10].

The V40 record of this project documents **LightDock 3 with DFIRE** for the
historical exploratory campaign. That attribution applies to the historical layer
and does not by itself identify the engine, version or scoring function of the
**reconstructed** 169-pair matrix reported in S4, which was recomputed later on a
rebuilt preparation; the distinction is preserved deliberately.

`[REQUIRED: the docking engine version, the scoring function and the parameter
file that produced the reconstructed 169 matrix scores are not identified in the
surviving records. The engine is not present in the environment in which the remaining
analyses were re-executed. This must come from the record of the original
execution or be reported as unrecorded; it is not inferred from whatever version
would install today.]`


---

## S10. References of this supplement

S1. JIMÉNEZ-GARCÍA, B. et al. LightDock: a new multi-scale approach to protein–protein docking. *Bioinformatics* 2018;34(1):49–55. doi:10.1093/bioinformatics/btx555

S2. ROEL-TOURIS, J.; BONVIN, A. M. J. J.; JIMÉNEZ-GARCÍA, B. LightDock goes information-driven. *Bioinformatics* 2020;36(3):950–952. doi:10.1093/bioinformatics/btz642

S3. ROEL-TOURIS, J.; JIMÉNEZ-GARCÍA, B.; BONVIN, A. M. J. J. Integrative modeling of membrane-associated protein assemblies. *Nat Commun* 2020;11:6210. doi:10.1038/s41467-020-20076-5

S4. ZHOU, H.; ZHOU, Y. Distance-scaled, finite ideal-gas reference state improves structure-derived potentials of mean force for structure selection and stability prediction. *Protein Sci* 2002;11(11):2714–2726. doi:10.1110/ps.0217002

S5. MINH, B. Q. et al. IQ-TREE 2: new models and efficient methods for phylogenetic inference in the genomic era. *Mol Biol Evol* 2020;37(5):1530–1534. doi:10.1093/molbev/msaa015

S6. KALYAANAMOORTHY, S. et al. ModelFinder: fast model selection for accurate phylogenetic estimates. *Nat Methods* 2017;14(6):587–589. doi:10.1038/nmeth.4285

S7. KATOH, K.; STANDLEY, D. M. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. *Mol Biol Evol* 2013;30(4):772–780. doi:10.1093/molbev/mst010

S8. GIULINI, M. et al. HADDOCK3: a modular and versatile platform for integrative modelling of biomolecular complexes. *J Chem Inf Model* 2025;65(13):7315–7324. doi:10.1021/acs.jcim.5c00969

S9. SHANNON, P. et al. Cytoscape: a software environment for integrated models of biomolecular interaction networks. *Genome Res* 2003;13(11):2498–2504. doi:10.1101/gr.1239303

S10. ABRAHAM, M. J. et al. GROMACS: high performance molecular simulations through multi-level parallelism from laptops to supercomputers. *SoftwareX* 2015;1–2:19–25. doi:10.1016/j.softx.2015.06.001

*Records transcribed from the V40 reference list without external validation of
DOIs. The main article carries its own list of eight references.*

---

## S11. Decision pending on supplementary results and figures

This file is **Supplementary Methods** only. The historical V40 document also
carries supplementary results of the focused LightDock campaign, the HADDOCK3
reassessment, and figures S1 to S11. Those are neither reproduced nor summarized
here, so the claim in the transfer map that no valid result was eliminated holds
for the methods transferred and cannot be demonstrated for that material from this
file alone.

**`[DECISION REQUIRED: either a separate Supplementary Results and Figures document
is prepared, carrying that material with its provenance caveats, or it is removed
from the scope of this publication. The second option is legitimate for a
reduced-scope methods paper but must be a stated choice, not an omission.]`**