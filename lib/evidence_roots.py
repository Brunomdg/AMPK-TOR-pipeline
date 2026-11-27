# >>> deposit header
# Purpose: Read the table that declares which folders of genome-wide HMMER results a run may read, check that
#          each folder is complete (210 domtbl and 210 tblout files), and list the folders, models and files.
# Methods: Candidate discovery
# Inputs:  `03_evidence_by_method/HMMER_EVIDENCE_ROOTS.tsv`; the model-list files it names; the `domtbl/` and
#          `tblout/` subfolders of each declared folder; the target set in TARGET_SET
# Outputs: the folders, the model-to-folder map and the file lists, returned (printed when run directly)
# Tools:   Python 3 standard library and the project module `targets`; versions: see envs/
# <<< deposit header
"""Which folders of HMMER results a run reads, declared in one table.

Why. The genome-wide search (03_evidence_by_method/hmmer_genomewide) scored every protein of the
210 proteomes against the 229-model panel. The twenty models adopted for the candidate replacement
targets were all in that panel, so their scores already exist; only the PP2C models built for
layer 4 (phosphatases) needed a new search, and re-running the whole 229-model array for one
family would take ten hours and 23 GB to reproduce files that exist. The PP2C models were
therefore searched on their own, with the same options, into a supplementary folder, and the
scripts that use the scores read a declared list of folders instead of one path.

The core eight targets read exactly the same files as before. The supplementary folder is read
only when TARGET_SET is not core8 (SET_SCOPE column), so every core8 run opens exactly the files it
opened before this mechanism existed. That follows from the construction, not from a measurement.

A model must be scored under exactly one folder (model_root stops otherwise): a model scored twice
against different panels would have two different E-values, and which one a reader saw would
depend on file order.

Each declared folder must hold 210 compressed files in both domtbl/ and tblout/; a folder that has
shrunk or is still being written stops the run instead of being read as complete.

Declaration table: 03_evidence_by_method/HMMER_EVIDENCE_ROOTS.tsv
    ROOT           folder holding domtbl/ and tblout/
    MODELS_FILE    file listing the model names scored under that folder, one per line
    SET_SCOPE      all, or a set name: which TARGET_SET values may read the folder
    PANEL_MD5      MD5 checksum of the .hmm file the folder was produced from, or - if unknown
    NOTE
"""
import csv
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import targets                                                # noqa: E402

DECL_REL = "03_evidence_by_method/HMMER_EVIDENCE_ROOTS.tsv"
EXPECT_SPECIES = 210


def declaration_path(root=None):
    return os.path.join(root or targets.project_root(), DECL_REL)


def roots(explicit=None, root=None):
    """List of dicts for the folders this TARGET_SET may read, in the row order of the
    declaration table, where the primary folder (SET_SCOPE=all) is meant to come first."""
    s = targets.set_name(explicit)
    p = declaration_path(root)
    if not os.path.exists(p):
        raise SystemExit("REFUSING: %s is missing; the evidence roots are declared, not guessed" % p)
    out = []
    with open(p, errors="replace") as fh:
        for r in csv.DictReader((l for l in fh if not l.startswith("##")), delimiter="\t"):
            r = {k: (v or "").strip() for k, v in r.items()}
            scope = r.get("SET_SCOPE") or "all"
            if scope != "all" and scope != s:
                continue
            d = r["ROOT"]
            if not os.path.isabs(d):
                d = os.path.join(root or targets.project_root(), d)
            r["ROOT"] = d
            for sub in ("domtbl", "tblout"):
                n = len(glob.glob(os.path.join(d, sub, "*.gz")))
                if n != EXPECT_SPECIES:
                    raise SystemExit("REFUSED: evidence root %s has %d %s files, expected %d. "
                                     "A root that shrank or is still being written must fail "
                                     "rather than be read as complete." % (d, n, sub, EXPECT_SPECIES))
            out.append(r)
    if not out:
        raise SystemExit("REFUSING: no evidence root declared for TARGET_SET=%s in %s" % (s, p))
    return out


def _models_of(r):
    mf = r.get("MODELS_FILE") or ""
    if not mf:
        return set()
    if not os.path.isabs(mf):
        mf = os.path.join(r["ROOT"], mf)
    if not os.path.exists(mf):
        raise SystemExit("REFUSING: MODELS_FILE %s for root %s does not exist" % (mf, r["ROOT"]))
    with open(mf, errors="replace") as fh:
        return {l.strip() for l in fh if l.strip() and not l.startswith("#")}


def model_root(explicit=None, root=None):
    """model name -> folder. Stops if a model is listed under two folders."""
    out = {}
    for r in roots(explicit, root):
        for m in _models_of(r):
            if m in out and out[m] != r["ROOT"]:
                raise SystemExit("REFUSING: model %s is scored under two roots (%s, %s); a "
                                 "reader could not tell which score it saw" % (m, out[m], r["ROOT"]))
            out[m] = r["ROOT"]
    return out


def files(kind, explicit=None, root=None):
    """(species, path) pairs over every folder, kind being domtbl or tblout, sorted within
    each folder and taken folder by folder; a species appears once per folder."""
    if kind not in ("domtbl", "tblout"):
        raise SystemExit("kind must be domtbl or tblout")
    out = []
    for r in roots(explicit, root):
        for f in sorted(glob.glob(os.path.join(r["ROOT"], kind, "*.gz"))):
            out.append((os.path.basename(f).split(".")[0], f))
    return out


def tblout_files(explicit=None, root=None):
    return [f for _, f in files("tblout", explicit, root)]


def domtbl_files(explicit=None, root=None):
    return [f for _, f in files("domtbl", explicit, root)]


if __name__ == "__main__":
    for r in roots():
        print("%s\tscope=%s\tmodels=%d" % (r["ROOT"], r.get("SET_SCOPE"), len(_models_of(r))))
