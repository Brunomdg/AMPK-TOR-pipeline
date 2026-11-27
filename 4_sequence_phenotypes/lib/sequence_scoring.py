#!/usr/bin/env python3
# >>> deposit header
# Purpose: Sequence functions shared by the scripts that measure the retained subunits: FASTA reading and writing,
#          MAFFT alignment, the state of a reference position in a copy, and the species summary over copies.
# Methods: Sequence phenotypes of the present subunits
# Inputs:  none of its own (sequences and FASTA files passed in by the calling scripts)
# Outputs: none of its own (alignment files written where the calling script asks)
# Tools:   Python 3 (standard library); MAFFT 7.526 [env record] from the orthofinder3 conda environment
# <<< deposit header
# Ran as: rbx1_score.py
## Sequence functions shared by the extraction, alignment and species steps of this chain (jobs x1, x2 and x3):
## FASTA reading and writing, the MAFFT procedure fixed before scoring, whether a position can be scored, and the
## species summary, which is computed the same way whatever the trait.
import hashlib, os, statistics, subprocess, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import phenotype_tables_chain1 as L  # noqa: E402

GAP = "-"


def read_fasta(path):
    """accession (first word of the header) -> sequence, with a trailing * removed as in the original AMPKα
    C-terminal analysis (folder `P03_actd`)."""
    out, acc, buf = {}, None, []
    with open(path, errors="replace") as fh:
        for ln in fh:
            if ln.startswith(">"):
                if acc:
                    out[acc] = "".join(buf).replace(" ", "").rstrip("*")
                acc = ln[1:].split()[0] if ln[1:].split() else ln[1:].strip()
                buf = []
            else:
                buf.append(ln.strip())
    if acc:
        out[acc] = "".join(buf).replace(" ", "").rstrip("*")
    return out


def write_fasta(path, seqs):
    with open(path, "w") as fh:
        for k, v in seqs.items():
            fh.write(">%s\n" % k)
            for i in range(0, len(v), 60):
                fh.write(v[i:i + 60] + "\n")
    return path


def seq_sha(s):
    return hashlib.sha256(s.encode()).hexdigest()


def mafft(seqs, out_path, pairwise=True):
    """The alignment procedure of the original AMPKγ pocket and AMPKα C-terminal analyses (folders `P02_gamma_pocket`
    and `P03_actd`): MAFFT --globalpair --maxiterate 1000. Returns the command line with the alignment."""
    with tempfile.NamedTemporaryFile("w", suffix=".faa", delete=False) as tf:
        tmp = tf.name
    write_fasta(tmp, seqs)
    cmd = [L.MAFFT, "--globalpair", "--maxiterate", "1000", "--anysymbol", "--quiet", tmp]
    with open(out_path, "w") as fh:
        r = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE, text=True)
    os.unlink(tmp)
    if r.returncode != 0:
        sys.exit("REFUSED: mafft failed: %s" % r.stderr[-500:])
    return " ".join(cmd[:-1]) + " <input>", read_fasta(out_path)


def ref_column_map(ref_aln):
    """reference residue index (1-based) -> alignment column index (0-based)."""
    m, k = {}, 0
    for i, c in enumerate(ref_aln):
        if c != GAP:
            k += 1
            m[k] = i
    return m


def span(qry_aln):
    """first and last alignment column carrying a residue of the query, or (None, None)."""
    idx = [i for i, c in enumerate(qry_aln) if c != GAP]
    return (idx[0], idx[-1]) if idx else (None, None)


def position_call(ref_aln, qry_aln, col):
    """Position state, exact match (1 or 0) and query residue at one alignment column.

    A position inside the aligned span of the query can be scored, whether it holds a residue or an internal gap.
    A position beyond either end of that span cannot be scored and is never counted as a lost residue.
    """
    lo, hi = span(qry_aln)
    if lo is None:
        return "NOT_CALLABLE_NO_ALIGNED_RESIDUE", 0, GAP
    if col < lo:
        return "NOT_CALLABLE_N_TERMINAL", 0, GAP
    if col > hi:
        return "NOT_CALLABLE_C_TERMINAL", 0, GAP
    q = qry_aln[col]
    if q == GAP:
        return "INTERNAL_DELETION", 0, GAP
    return ("MATCH", 1, q) if q.upper() == ref_aln[col].upper() else ("SUBSTITUTION", 0, q)


CALLABLE_STATES = ("MATCH", "SUBSTITUTION", "INTERNAL_DELETION")


def aggregate(values):
    """Species summary over the ortholog copies that have a value: median, minimum, maximum, number of copies and
    number of distinct values. No copy is chosen by its own value."""
    v = sorted(float(x) for x in values)
    if not v:
        return {"MEDIAN": None, "MIN": None, "MAX": None, "N_CALLABLE": 0, "N_UNIQUE": 0}
    return {"MEDIAN": float(statistics.median(v)), "MIN": v[0], "MAX": v[-1],
            "N_CALLABLE": len(v), "N_UNIQUE": len(set(v))}


def query_pos_at_column(qry_aln, col):
    """1-based residue index of the query at an alignment column, or None for a gap."""
    if qry_aln[col] == GAP:
        return None
    return sum(1 for c in qry_aln[:col + 1] if c != GAP)


def last_query_pos_before(qry_aln, col):
    """1-based index of the last query residue at or before a column (0 if none)."""
    return sum(1 for c in qry_aln[:col + 1] if c != GAP)
