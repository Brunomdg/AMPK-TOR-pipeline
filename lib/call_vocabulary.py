#!/usr/bin/env python3
# >>> deposit header
# Purpose: Define the one set of call values shared by every line of evidence (MEMBER, NOT_MEMBER and the
#          no-opinion values), with helpers that count only real calls and reject unknown values; run
#          directly, test itself.
# Methods: Candidate discovery
# Inputs:  call values passed by the calling script
# Outputs: filtered calls and counts, returned; a self-test report when run directly
# Tools:   Python 3 standard library; versions: see envs/
# <<< deposit header
"""One set of call values for every line of evidence, so that "no opinion" cannot be read as a
vote against.

Why this exists. Two errors on 2 September 2026 had the same cause, private copies of this set of
values that had drifted apart:

  1. A display truncated CALL_kinortho to two characters. "NO_LAYER" was printed as "NO", read as
     a negative call, and an argument was built on KinOrtho disagreeing on four of five Toxoplasma
     proteins. KinOrtho runs on AMPKα only; it is NO_LAYER for all seven other targets and had
     disagreed with nothing.
  2. A test requiring every line of evidence to agree excluded NO_COVERAGE, NOT_APPLICABLE and
     NO_LAYER from its set of no-opinion values but not UNRESOLVED, so "could not decide" was
     counted as "no".

Eight scripts each carried their own copy of this set, which is how the copies drift. This module
follows accessions.py, which solved the same problem for accession separators: define the set
once, import it everywhere, and provide a check that tests the assumption instead of asserting it.

The rule enforced here: a line of evidence with no opinion never counts as a vote against.

Why assert_known exists. A value that is neither a call nor a known no-opinion value is the
dangerous case: any `x == "MEMBER"` test silently treats it as not MEMBER, so it counts as a vote
against. assert_known() stops that from happening quietly: a new method that emits a new value
fails loudly when its values are checked, instead of skewing a count months later.
"""

MEMBER = "MEMBER"
NOT_MEMBER = "NOT_MEMBER"

#: Real calls. Only these two ever count in a numerator or a denominator.
VERDICT = frozenset({MEMBER, NOT_MEMBER})

#: Every way a line of evidence says "no opinion". Never a vote against.
#: UNRESOLVED belongs here: the method ran and could not decide, which is not
#: the same as deciding no.
NO_OPINION = frozenset({
    "NO_COVERAGE",      # the method never saw this protein
    "NOT_APPLICABLE",   # the method does not apply to this target
    "NO_LAYER",         # the method was not run for this target at all
    "UNRESOLVED",       # the method ran and could not decide
    "",
})


def assert_known(values, where=""):
    """Raise an error on any value that is neither a call nor a known no-opinion value.

    An unknown value is worse than a wrong one: `x == "MEMBER"` silently counts it as a
    vote against. Fail loudly instead.
    """
    unknown = {v for v in values
               if v is not None and v not in VERDICT and v not in NO_OPINION}
    if unknown:
        raise ValueError(
            "unknown layer value(s) %s%s -- neither a verdict nor a recognised "
            "no-opinion token. Add it to call_vocabulary deliberately; do not let an "
            "equality test file it as dissent."
            % (sorted(unknown), (" in " + where) if where else ""))
    return True


def has_opinion(v):
    return v is not None and v not in NO_OPINION


def is_positive(v):
    return v == MEMBER


def is_negative(v):
    """True only for an explicit NOT_MEMBER; no opinion is not negative."""
    return v == NOT_MEMBER


def verdicts(row, cols, where=""):
    """The real calls among cols, dropping every no-opinion value."""
    vals = [row.get(c) for c in cols]
    assert_known(vals, where)
    return [v for v in vals if has_opinion(v)]


def tally(row, cols, where=""):
    """(positive, total), counting only the lines of evidence that made a call."""
    v = verdicts(row, cols, where)
    return sum(1 for x in v if x == MEMBER), len(v)


if __name__ == "__main__":
    ok = 0

    def check(desc, got, want):
        global ok
        good = got == want
        ok += good
        print("  %-62s %s" % (desc, "PASS" if good else "FAIL got=%r" % (got,)))

    print("call_vocabulary self-test")
    check("NO_LAYER has no opinion (the KinOrtho error)",
          has_opinion("NO_LAYER"), False)
    check("UNRESOLVED has no opinion (the ALL-layers error)",
          has_opinion("UNRESOLVED"), False)
    check("NO_LAYER is not negative", is_negative("NO_LAYER"), False)
    check("UNRESOLVED is not negative", is_negative("UNRESOLVED"), False)
    check("NOT_MEMBER is negative", is_negative(NOT_MEMBER), True)
    r = {"a": MEMBER, "b": "NO_LAYER", "c": "UNRESOLVED", "d": NOT_MEMBER}
    check("tally counts only real verdicts",
          tally(r, ["a", "b", "c", "d"]), (1, 2))
    try:
        assert_known(["MEMBER", "PROBABLY"])
        check("assert_known rejects an unknown token", False, True)
    except ValueError:
        check("assert_known rejects an unknown token", True, True)
    print("  %d of 7 passed" % ok)
    raise SystemExit(0 if ok == 7 else 1)
