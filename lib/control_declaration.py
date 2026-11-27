#!/usr/bin/env python3
# >>> deposit header
# Purpose: Return the path and version label of the current (highest-numbered) declared-control table,
#          `RAW_identity-controls_declaration[_vN].tsv`, so that every script reads the same version.
# Methods: several steps (shared helper)
# Inputs:  the folder that holds the control tables (when run directly: the folder in `MS210_OUT`)
# Outputs: the path and version label, returned (printed when run directly)
# Tools:   Python 3 standard library; versions: see envs/
# <<< deposit header
# Ran as: rl_declaration.py
"""Which version of the declared-control table is current?

Each round of review writes a new numbered version of the control table instead of editing the
old one, so versions v1, v2, v3 and so on accumulate. Every script that reads the table must use
the same version, and writing a version number into each script is how they drift apart: one
script could rank the evidence sources against v3 while another scored tree tips against v2, with
nothing to show it.

The rule is simply the highest version present. Callers are told which file they received, so it
can be named in the header of their own output.
"""
import os
import re

BASE = "RAW_identity-controls_declaration"


def current(outdir):
    """(path, label) of the highest-numbered control table in outdir; a file without a _vN
    suffix counts as version 1."""
    best, bestn = None, -1
    for name in os.listdir(outdir):
        if not name.startswith(BASE) or not name.endswith(".tsv"):
            continue
        m = re.match(re.escape(BASE) + r"(?:_v(\d+))?\.tsv$", name)
        if not m:
            continue
        n = int(m.group(1)) if m.group(1) else 1
        if n > bestn:
            best, bestn = name, n
    if best is None:
        raise SystemExit("no control declaration found in %s" % outdir)
    return os.path.join(outdir, best), "v%d" % bestn


if __name__ == "__main__":
    p, lab = current(os.environ["MS210_OUT"])
    print("%s\t%s" % (lab, p))
