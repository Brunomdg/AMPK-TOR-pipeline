# The limits of this deposit

`docs/dating_gap.md` states the one step whose script is lost. This page states the rest, so that none of it
has to be discovered by running something and watching it fail.

## The model engine is not the revision that produced the published fits

`5_phylogenetic_analyses/lib/model_engine.R` is the engine whose preamble the AMPKα tail models load:
`lib/tail_model_lib.R` reads it, cuts it at the `## ---------------- Analysis 1` marker, evaluates what comes
before, and uses `fit_cont`, `row_of`, the chronogram, the seed and the three species tables it defines.

It was absent from the first deposit, which is why those two model scripts could not run at all. It is here
now, with one limit.

| | sha256 |
|---|---|
| The revision the tail models were pinned against | `526dc52ae74a87bc46a374ec55c4b83325d6ea8a9706b85fb0be49055f3abe71` |
| The revision that survives, and is deposited here | `ef3350ad8eb85b33bbead953add1f94c18f94f197aa1b8a27539b7efc384bec7` |

The two differ, and the surviving file is dated two days after every other file in the folder it came from, so
it was edited after the run. This is the same situation as the `chronos` script, with one difference: here the
file survives, and what it defines still satisfies every check the tail models made of it. The seed is
488761235, `fit_cont` is the same `phylolm` call, the tree is the 209-tip chronogram, and the preamble marker
is in place. So the models run, and they run against a later revision of the engine than the one that produced
the numbers in the paper. Nothing was reconstructed to make them run.

Its own Analysis 1 is superseded. The nine prespecified state models it fits were refitted by the audited
penalised multistart pipeline, and a guard on the first line of that section stops the file from overwriting
them. The tail models never reach it: they evaluate only the preamble.

## The integrity gates are not here

Four files were removed from the deposit, because checksum verification is not part of reproducing a result:
the SHA-256 manifest of the reference inputs, the write guard on the original cluster folders, and the two
gates that held a step until the CAMKK and LKB1 calls were declared final. Removing them also removed two
checks that compared a built output against a stored copy of it: the CAMKK and LKB1 columns of the
species-state matrix, and the same two traits in the standard coding, were each compared against a frozen
state vector. Those comparisons computed nothing. A run now stops when an input is missing, and no longer when
an input differs from a stored hash.

## What the scripts name and this deposit does not hold

Some scripts name a file that is not here. None of them loads one: the references are in comments, in headers,
or in provenance strings. They are listed so that a reader does not look for them.

| Named in | What it is |
|---|---|
| `lib/tail_model_lib.R`, `5_phylogenetic_analyses/core_complement_models.R` | `rbx2_models.R` and `rbx2c_loco.R`, the sister chains of the model engine, named in a provenance line and a comment |
| `4_sequence_phenotypes/lib/phenotype_tables_chain1.py` | `rl_r1_three_state_release.py`, cited for the rule that marks a copy model sensitive. The rule is implemented in the library itself |
| `5_phylogenetic_analyses/lib/ancestral_state_lib.R` | `rl_rm_asr.R`, the script this library's functions were copied from, and the file it was validated against |
| `6_structure_prediction_and_pocket/lib/amp_site_measurement.py` | `PLAN.md` and three `LITERATURE_BASIS` documents, the design notes behind the pocket measurement |
| `lib/targets.py` | `rl_target_registry.py`, which built the target registry the module reads |

## One external command line is missing

`external_runs/` holds one file per program and step. The three AMPKγ copies recovered on 22 September 2026
were found by profile searches with `phmmer` and CDD. The CDD command line is deposited
(`2_candidate_discovery_and_orthology/external_runs/rpsblast_cdd.sh`); the `phmmer` command line is not among
the records.

## The three recovered copies

The scripts here are the ones that made the results. `6_structure_prediction_and_pocket/measure_amp_sites.py`
covers 236 AMPKγ copies and five controls, where the Supplementary Methods reports 239 copies: the difference
is the three copies recovered after that run. The reruns that carried them into the AMPKγ-dependent results
are not included, and the data deposit holds the final tables. `docs/code_map.md`, section 3, names the three
copies.
