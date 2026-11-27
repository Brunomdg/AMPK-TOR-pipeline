# The TOR N-terminal HEAT structural screen: references and acceptance rule

Sources: `4_sequence_phenotypes/tor_heat_structural_screen.py` for the rule, which is the authority for what
the code does, and the Supplementary Methods, sequence phenotypes of the present subunits, for the identity
of the reference structures. Supplementary Methods subsection: sequence phenotypes of the present subunits.

## The reference structures

Each copy's N-terminal crop is aligned to these by sequence-independent structural alignment with US-align
20241201 in protein mode.

| Accession | What it is | Role |
|---|---|---|
| `5H64` | human mTOR | positive reference |
| `5FVM` | *Kluyveromyces marxianus* TOR | positive reference |
| `7PQH` | *Saccharomyces cerevisiae* TOR2 | positive reference |
| `1JDH` | β-catenin armadillo repeats | negative control on the geometric criterion |
| `1A17` | PP5 tetratricopeptide repeats | negative control |
| `2ZTA` | GCN4 coiled coil | negative control |
| `1C1G` | tropomyosin coiled coil | negative control |

The script reads them from `CONTROL_MANIFEST_EXTENDED.tsv`, which carries each reference's `CONTROL` label,
its `STRUCTURE` path and its `EXPECTED_CLASS`, and which is in the data deposit. Two self-comparisons are
skipped by name: a *Homo sapiens* query is never scored against `HUMAN_TOR_HEAT`, and a *Saccharomyces
cerevisiae* query is never scored against `SCEREVISIAE_TOR2_HEAT` (lines 53 and 54).

## The acceptance rule

Two independent gates, both fixed in code.

**Geometric gate** (line 75). Both length-normalised TM-scores at or above 0.5, and the aligned coverage at or
above 0.5 on *both* structures:

```text
min(tm) >= 0.5 and min(aligned/len_query, aligned/len_reference) >= 0.5
```

Taking the minimum of the two normalisations is what makes the gate symmetric: a short crop cannot pass by
matching a fragment of a long reference.

**Confidence gate.** At least 100 close pairs, where close is an aligned pair under 5 Å, with median pLDDT at
or above 70 over those pairs and at least 70 percent of them at or above 70. A copy counts as structurally
supported only when both gates pass (`supported_match`, line 77).

Both rules are also written verbatim into the pilot's own output,
`direct_alignment_v1/PRODUCTION_CROP_PREFLIGHT.json` (line 92), which is how they survive outside the code.

## The preflight

The screen refuses to run in production until a three-copy pilot passes: every positive reference must pass
both the geometric and the pLDDT rules, and no negative control may pass the geometric rule (lines 25, 89 and
91). The production run asserts that `PRODUCTION_CROP_PREFLIGHT.json` records `passed`, so a failed pilot
stops the screen rather than being reported alongside it.

## What this gate is not

The crop itself is chosen by transferring human HEAT alignment coordinates and capping before the observed FAT
domain. The script records that in every summary it writes: region selection only, not proof of HEAT. The
structural screen tests the crop; it does not define it. The final present or absent call is made downstream
by `4_sequence_phenotypes/tor_heat_final_calls.py`, on the accessions in
`tor_heat_evidence_accessions.tsv`.
