#!/usr/bin/env python3
# >>> deposit header
# Purpose: Join two sets of protein accessions from different tools after normalisation, and raise an error
#          unless a minimum number of known-present control keys is found; run directly, test itself.
# Methods: Candidate discovery
# Inputs:  the two key sets and the control keys passed by the calling script
# Outputs: the matched keys and the number of controls found, returned; a self-test report when run directly
# Tools:   Python 3 standard library and the project module `rl_accession`; versions: see envs/
# <<< deposit header
"""A join across tools that must show it works before its result is believed.

Why this exists. On 2 September 2026 a join between the project's protein identifiers and the raw
InterProScan output returned no hits at all for ten CDD signatures, including cd14079, which a
separate table already showed to hit 27 to 30 of the 45 AMPKα controls. The zero was reported as
absence from the data and interpreted as biology. The real cause was a fourth accession convention:
InterProScan writes "Species||ACCESSION", in addition to ':' (FastTree truncation), '@' and '_'
(IQ-TREE, OrthoFinder) and '__' (the tip-name convention). Joining on the raw field matched nothing.

The rule. A join that returns zero must demonstrate that the zero is real; it may not simply report
it. Every join carries a positive control: a key known to be present, with a minimum number that
must be found. If the control is not found, the join raises an error and prints the raw and
normalised forms from both sides, because the join is then broken and its output says nothing
about the data.

This follows the same pattern as rl_done in rl_assert_output.sh (no marker without its output) and
rl_job_state (no job state that was not read from SLURM): a check that fails, not a rule to
remember.
"""
import sys


# One definition, in rl_accession. This module used to carry a second normalise() that differed
# from rl_accession's by the removal of the '||' prefix, so the project had two functions that each
# claimed to be the single accession normalisation. The behaviour kept is this module's (the
# superset); it now lives in rl_accession and is re-exported here so existing imports keep working.
from rl_accession import normalise  # noqa: F401


class JoinFailed(RuntimeError):
    pass


def proven_join(left_keys, right_keys, control, expect_min, what="join"):
    """Intersect two key sets after normalisation, refusing a zero that has not been shown real.

    control     a key that must be found, or an iterable of such keys
    expect_min  the minimum number of control keys that must be found
    """
    L = {normalise(k): k for k in left_keys}
    R = {normalise(k): k for k in right_keys}
    hit = set(L) & set(R)

    ctrl = [control] if isinstance(control, str) else list(control)
    landed = [c for c in ctrl if normalise(c) in hit]
    if len(landed) < expect_min:
        raise JoinFailed(
            "%s: positive control failed -- %d of %d control keys landed, "
            "expected at least %d.\n"
            "  example control      : %r\n"
            "  normalised to        : %r\n"
            "  a left key raw       : %r -> %r\n"
            "  a right key raw      : %r -> %r\n"
            "The join is broken; its zero says NOTHING about the data."
            % (what, len(landed), len(ctrl), expect_min,
               ctrl[0] if ctrl else None,
               normalise(ctrl[0]) if ctrl else None,
               next(iter(left_keys), None), next(iter(L), None),
               next(iter(right_keys), None), next(iter(R), None)))
    return {L[k] for k in hit}, len(landed)


if __name__ == "__main__":
    ok = 0

    def check(desc, fn, want_raise):
        global ok
        try:
            fn()
            good = not want_raise
        except JoinFailed:
            good = want_raise
        ok += good
        print("  %-62s %s" % (desc, "PASS" if good else "FAIL"))

    print("rl_join self-test")
    ipr = ["Homo_sapiens||NP_996790.3", "Mus_musculus||NP_001013385.1"]
    proj = ["NP_996790.3", "NP_001013385.1"]

    check("the real InterProScan convention joins after normalisation",
          lambda: proven_join(ipr, proj, "NP_996790.3", 1, "ipr"), False)
    check("a broken join RAISES instead of returning zero",
          lambda: proven_join(["Species||AAA"], ["ZZZ"], "AAA", 1, "broken"), True)
    check("colon truncation (FastTree) is handled",
          lambda: proven_join(["TvY486_0200210:mRNA-p1"],
                              ["TvY486_0200210_mRNA-p1"],
                              "TvY486_0200210:mRNA-p1", 1, "colon"), False)
    m, n = proven_join(ipr, proj, ["NP_996790.3", "NP_001013385.1"], 2, "both")
    ok += (len(m) == 2 and n == 2)
    print("  %-62s %s" % ("both controls land and both keys are returned",
                          "PASS" if (len(m) == 2 and n == 2) else "FAIL"))
    print("  %d of 4 passed" % ok)
    raise SystemExit(0 if ok == 4 else 1)
