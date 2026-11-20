# The experimental structures behind the AMPKγ AMP-contact positions

Source: `11_structure_pocket/panel_86_queries/p45b_colabfold_analysis.py`, the measurement library that both
structure scripts load after checking its checksum. This page states what that script fixes; the script remains
the authority. Supplementary Methods subsection: structure prediction and AMPKγ pocket geometry.

## The structures the repository records

| Accession | Chain | Role | Where it is fixed |
|---|---|---|---|
| `2V8Q` | E, rat γ1 | **The reference frame of the geometry measurement.** The measured site neighbourhoods are taken from it, the superposition is onto it, and the AMP ligands transplanted into each model are its own. It does not define the ten contact positions; see below | lines 3 to 10, 48, 63 to 76 |
| `4CFE` | E, human γ1 | **Structure control.** Ligands stripped, scored like a predicted model with the pLDDT column set to 100 | lines 16, 135 to 138 |

Two reference sequences are used for the numbering map, not as structures: `P54619`, canonical human γ1, which
carries the numbering everything is reported in, and `P80385`, rat γ1 (lines 48 to 51). Structures are
downloaded from RCSB as mmCIF at run time (`fetch`, lines 38 to 41).

A third accession appears elsewhere in the repository: `1Z0N` is the rat structure whose numbering the human
AMPKβ reference carries, in `05_sequence_phenotypes/rbx2/rbx2_ampkb_cbm.py`. It is a different subunit and has
no part in the AMPKγ AMP sites.

## Two different things, at two different cutoffs

The contact **positions** and the predicted-model **geometry** do not come from the same place, and the
difference matters when reading either one.

| | How it is defined | Cutoff | Where it lives |
|---|---|---|---|
| The ten site 1 contact positions | a consensus over **six crystal structures** | 4.0 Å | the frozen panel `AMPKg_contact_reference.tsv`, published as the manuscript's table of AMP-contact reference positions |
| The geometry measured in each predicted model | the γ residues around each AMP in **`2V8Q`** chain E | 5 Å | `p45b_colabfold_analysis.py`, lines 3 to 10 |

**The six crystal structures of the position consensus are not named in this repository.** `2V8Q` and `4CFE`
are the only PDB accessions in the structure code, and neither is used to define the ten positions: `2V8Q`
supplies the frame, the superposition target and the transplanted AMP ligands, and `4CFE` is a control. The
accessions of the six, and the rule that turned them into a consensus, are in the manuscript's own table of
AMP-contact reference positions and in the Supplementary Methods. They are not reconstructed here.

What the repository can confirm about that panel: each of its three site sets holds exactly ten positions, and
`rbx1_gamma.py` refuses to run otherwise; the anchor numbering is canonical PRKAG1 plus 9, as the script's own
comment states and as the published panel's two numbering columns show; and site 1 is primary, with sites 3
and 4 declared in advance as site-definition sensitivities.

## How the geometry is measured against `2V8Q`

1. The sites are the γ residues within 5 Å of each AMP in `2V8Q` chain E, numbered on canonical human γ1
   (`P54619`).
2. Site 1 is the AMP whose site overlaps the ten official site-1 contacts. That list of ten positions is read at
   run time from `sites.json` (`SITES["official_site1"]`, lines 81 and 112), which is part of the data deposit,
   not of this repository. The three site sets scored per copy (site 1, site 3, site 4) are likewise read from
   the frozen reference panel by `05_sequence_phenotypes/rbx1_chain/rbx1_gamma.py`, which refuses to run unless
   each set holds exactly ten positions.
3. Per predicted model, a site is assessable only when every one of its residues has pLDDT at or above 70;
   otherwise it is not scored.
4. Geometry: the site Cα atoms are superposed on `2V8Q` (RMSD), the AMP is transplanted with the same transform,
   and the minimum heavy-atom distance to it is measured. Cavity volume is measured by pyKVFinder in a box fixed
   in the `2V8Q` frame, the site residues' heavy-atom extent plus 3.5 Å of padding, in ligand mode with a 5 Å
   cutoff around the transplanted AMP.
