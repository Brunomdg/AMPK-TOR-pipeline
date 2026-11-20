# The AMPKβ carbohydrate-binding module: reference, positions and scoring rule

Source: `05_sequence_phenotypes/rbx2/rbx2_ampkb_cbm.py`. This page states what that script fixes; the script
remains the authority. Supplementary Methods subsection: sequence phenotypes of the retained subunits.

## The reference

| Item | Value | Where it is fixed |
|---|---|---|
| Reference protein | `NP_006244.2`, human PRKAB1 | `REF_ID`, line 30 |
| Reference sequence | read from `CLEAN_literature_controls.tsv`, the frozen literature controls | `reference()`, lines 40 to 46 |
| Numbering | the numbering of the rat structure 1Z0N at offset zero | the script's own note, lines 14 and 15 |
| Module envelope on the reference | residues 78 to 161 | `CBM_FROM`, `CBM_TO`, line 31; taken from the genome-wide HMMER domain table |
| Primary alignment | each copy against the reference, pairwise, MAFFT `--globalpair --maxiterate 1000 --anysymbol` | `rbx1_score.py`, line 57 |
| Declared sensitivity | one joint alignment of the reference with every copy | `main()`, line 116 |

## The two site sets

| Site set | Positions | What it is |
|---|---|---|
| `CBM_CORE5` | 5 | the binding positions; the score of record |
| `CBM_STRUCTURAL7` | 7 | the structural positions; a declared site-definition sensitivity |

**The residue positions themselves are not in this repository.** They are read at run time from
`05_ampkb_cbm/AMPKb_CBM_reference.tsv`, which carries `SITESET`, `SITE_ID`, `REF_POSITION` and `REF_AA` and is
part of the data deposit. The script refuses to run unless that file holds exactly 5 and 7 positions
(`sites()`, line 52), and refuses again unless the reference residue at every `CBM_STRUCTURAL7` position is the
`REF_AA` the file declares (the numbering gate, lines 108 to 110). The manuscript's supplementary table of
AMPKβ reference-position states is the published form of this panel.

## The scoring rule

For each retained AMPKβ copy, under each alignment (`PAIRWISE`, the primary, and `JOINT`, the sensitivity):

1. At every position of the site set, call the aligned residue. A position is **callable** when its state is one
   of the callable states of `rbx1_score.position_call`; an internal deletion is callable and scores 0, while a
   position beyond the observed sequence span is not callable and is never scored as a lost residue.
2. Count the callable positions (`N_CALLABLE`) and, among them, those identical to the reference residue
   (`N_EXACT`).
3. **The score is emitted only when every position of the set is callable**: `CBM_CORE5_SCORE_0_5` is `N_EXACT`
   when `N_CALLABLE` is 5, and `NA` otherwise; `CBM_STRUCTURAL7_SCORE_0_7` is `N_EXACT` when `N_CALLABLE` is 7,
   and `NA` otherwise. Partial information stays in `N_CALLABLE` and is not turned into a score.
4. Background: identity at every callable position of the module envelope (78 to 161) that is not one of the
   seven structural positions, on the primary alignment. `CBM_CORE5_FRACTION` is `N_EXACT / N_CALLABLE` of the
   core five on the primary alignment.
5. The pairwise and joint scores are never averaged and never chosen per protein. Their difference and whether
   they agree are both recorded.

## The species summaries

One row per species, over that species' copies: median, minimum and maximum of the pairwise core-five score,
the number of callable copies, whether the copies disagree, the median of the joint score, the median of the
structural-seven score, the median over high-provenance copies only, the median excluding model-sensitive
copies, and the median background identity. The three variants in that list are the ones the tail models swap
in as sensitivities; the AMPKγ equivalents are in `site1_score_input_variants.tsv`.
