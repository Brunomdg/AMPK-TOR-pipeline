# The experimental structures behind the AMPKγ AMP-contact positions

Two separate things rest on crystal structures, at two different cutoffs, and only the second is computed by
code in this repository. Conflating them is easy and gives the wrong reading of both.

| | What it fixes | Cutoff | Structures | Computed where |
|---|---|---|---|---|
| **The position consensus** | the three ten-position site definitions scored on every copy | 4.0 Å, heavy atom | **six**, listed below | not in this repository; the result is the frozen reference panel |
| **The geometry measurement** | what is measured in each predicted model | 5.0 Å, neighbour search | `2V8Q`, with `4CFE` as control | `amp_site_measurement.py` |

Source for the first row: the Supplementary Methods, sequence phenotypes of the present subunits, and the
manuscript's table of AMP-contact reference positions. Source for the second: the script named, which remains
the authority for what it does.

## The six structures of the position consensus

All six are the mammalian γ1 subunit.

| Accession | Species | Reference |
|---|---|---|
| `4CFE` | human | Xiao et al. 2013 |
| `4CFF` | human | Xiao et al. 2013 |
| `4RER` | human | Li et al. 2015 |
| `4ZHX` | human | Langendorf et al. 2016 |
| `4EAI` | rat | Chen et al. 2012 |
| `2V8Q` | rat | Xiao et al. 2007 |

Contact positions were computed from these six at a 4.0 Å heavy-atom cutoff. The computation returns a
consensus of 32 residues over the three occupied sites, from which three sets of ten were selected, one per
site. Those 30 selections are the rows of the reference table. Site 1 is primary; sites 3 and 4 were scored
under identical rules as sensitivities declared in advance, and the three sets are reported separately and
never pooled.

**This repository does not perform that computation.** It consumes its result: `ampkg_contact_scores.py` reads the
frozen panel `AMPKg_contact_reference.tsv`, which is in the data deposit, and refuses to run unless each of
its three site sets holds exactly ten positions.

## Numbering, and how the offset was checked

Published site positions are numbered on canonical human γ1, UniProt `P54619`, 331 residues. The scoring code
indexes into RefSeq `NP_001193638.1`, which is nine residues longer, so the anchor index is the canonical
position plus 9. `ampkg_contact_scores.py` states that offset in its own comment (`ANCHOR`, line 29).

The offset was verified by recovering 29 of the 30 published site residues with it, against 1 of 30 without
it. The reference table carries that check per row, which is why exactly one of its 30 rows is marked
unverified.

## The geometry measurement, in this repository

`6_structure_prediction_and_pocket/lib/amp_site_measurement.py` is the measurement library that both
structure scripts load after checking its checksum.

1. Three AMP sites are defined on chain E of `2V8Q` by a 5.0 Å neighbour search around every atom of the bound
   AMP, at the positions occupied by `AMP1327` (CBS3), `AMP1328` (CBS1, the primary site, which carries all
   ten site 1 contacts) and `AMP1329` (CBS4). The script iterates over the AMP ligands of that chain and names
   each one's extracted ligand file by its residue number (lines 73 to 76).
2. Residues are reported on canonical human γ1 (`P54619`); `P80385`, rat γ1, is used to build the numbering
   map onto the structure (lines 48 to 51).
3. Per predicted model, a site is assessable only when every one of its residues has pLDDT at or above 70.
4. The site Cα atoms are superposed on `2V8Q`, the reference AMP is transferred with the same transform, and
   the script measures the Cα RMSD after superposition, the shortest heavy-atom distance from the model to the
   transferred ligand, the minimum and mean pLDDT over the site residues, and how many site residues mapped.
5. Cavity volume is measured by pyKVFinder in a box fixed in the `2V8Q` frame: the site residues' heavy-atom
   extent plus 3.5 Å of padding, in ligand mode with a 5 Å cutoff around the transferred AMP.
6. A ligand-stripped `4CFE` is the known-answer control for the site transfer, scored like a predicted model
   with the pLDDT column set to 100 (lines 135 to 138).

Measurements are summarised first within each query, as the median over its 25 models, then across queries. No
confidence threshold excludes a query, and a measurement that could not be made is reported as unavailable,
never as zero.

Note that `4CFE` and `2V8Q` appear in both roles: each is one of the six structures of the position consensus
and also serves the geometry measurement, `2V8Q` as the frame and `4CFE` as the control. The other four
structures of the consensus play no part in the code here.

## A third accession, for a different subunit

`1Z0N` is the rat structure whose numbering the human AMPKβ reference carries, in
`4_sequence_phenotypes/ampkb_module_score.py`. It has no part in the AMPKγ AMP sites.
