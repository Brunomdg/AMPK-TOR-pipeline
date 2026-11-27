#!/usr/bin/env python3
## New code for the old ColabFold + Fpocket analysis, step b (analysis). Design (PLAN.md row 18; NOCODE_LITERATURE_BASIS.md,
## ColabFold + Fpocket): sites = gamma residues within 5 A of each AMP in 2V8Q (rat gamma1, chain E), numbered on canonical
## human gamma1 (P54619); site 1 = the AMP whose site overlaps the ten official site-1 contacts. Per model (25 per query):
## gate = every site residue pLDDT >= 70, else "not assessable" (Tunyasuvunakool 2021); residue geometry primary: site
## Calpha superposition on 2V8Q (RMSD), the AMP transplanted with the same transform, minimum heavy-atom distance to it and,
## for site 1, each contact's residue and nearest-atom distance; volume secondary (POCKET_VOLUME_LITERATURE_BASIS.md; Bruno:
## "pyKVFinder, fixed box"; box and restriction from IMPLEMENTATION_CHOICES_LITERATURE_BASIS.md): pyKVFinder
## (deterministic grid, default probes and step) on the model moved into the 2V8Q frame by the site superposition, with the
## B-factor column set to 0 (ColabFold writes pLDDT there). The site region is the 2V8Q AMP of that site in ligand mode,
## ligand_cutoff 5 A (the package default; the author's recommendation for comparing one binding site across proteins,
## GitHub issue #199; TRAPP's 5 A), so the volume covers that pocket only; the residue box (site residues' extent + 3.5 A)
## stays as the outer grid so every model shares one lattice (the residue box alone held half or more of a neighbouring
## nucleotide site). A site with no cavity is volume 0, counted.
## Fpocket is not used: its 4.2.3 volume is clock-seeded Monte Carlo with no seed option, and it filters on B-factors.
## Structure controls: 2V8Q and 4CFE gamma chains with ligands removed, scored like models (pLDDT column set to 100).
## Environment: p68_kvf (python 3.11, pyKVFinder 0.9.5, numpy 2, biopython, pandas).
## Usage: p45b_colabfold_analysis.py sites | analyze <task> <ntasks> | determinism | summarize
import os, sys, re, csv, json, glob, subprocess, statistics, urllib.request, shutil, math
import numpy as np
from Bio.PDB import MMCIFParser, PDBParser, PDBIO, Select, NeighborSearch
from Bio import Align
from Bio.Align import substitution_matrices
S = "/gpfs2/scratch/bdigenov/ms210_phase2/TOR-AMPK_orthologs_10target_R1"; P = S + "/20_port_68taxa_analyses"
G = P + "/02_derived/gamma"; W = P + "/02_derived/p45_colabfold"; R = P + "/03_results/p45_colabfold"; ST = P + "/02_derived/structures"
O = P + "/03_results/p45_analysis"; os.makedirs(O, exist_ok=True); os.makedirs(ST, exist_ok=True)
PAD = 3.5   ## pyKVFinder residue-box padding default (data/tests/residues-box.toml)
THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
         "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
rows = lambda p: list(csv.DictReader([l for l in open(p) if not l.startswith("#")], delimiter="\t"))
AL = Align.PairwiseAligner(); AL.mode = "global"; AL.open_gap_score = -10; AL.extend_gap_score = -0.5; AL.substitution_matrix = substitution_matrices.load("BLOSUM62")
def posmap(a, b):   ## 1-based position in a -> 1-based position in b
    aln = AL.align(a, b)[0]; m = {}
    for (s1, e1), (s2, e2) in zip(*aln.aligned):
        for k in range(e1 - s1): m[int(s1 + k + 1)] = int(s2 + k + 1)
    return m
uni = lambda acc: "".join(urllib.request.urlopen("https://rest.uniprot.org/uniprotkb/%s.fasta" % acc, timeout=60).read().decode().split("\n")[1:])
def fetch(pdbid):
    f = ST + "/%s.cif" % pdbid
    if not os.path.exists(f): urllib.request.urlretrieve("https://files.rcsb.org/download/%s.cif" % pdbid, f)
    return MMCIFParser(QUIET=True).get_structure(pdbid, f)[0]
def protein_residues(chain): return [r for r in chain if r.id[0] == " " and r.get_resname() in THREE]
def kabsch(Pm, Qm):   ## rotation/translation taking Pm onto Qm
    pc, qc = Pm.mean(0), Qm.mean(0); H = (Pm - pc).T @ (Qm - qc); U, s, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T)); D = np.diag([1, 1, d]); Rm = Vt.T @ D @ U.T
    return Rm, qc - Rm @ pc
if sys.argv[1] == "sites":
    hs1 = uni("P54619"); rat = uni("P80385"); st = fetch("2V8Q"); E = st["E"]
    res = protein_residues(E); seq = "".join(THREE[r.get_resname()] for r in res)
    e2rat = posmap(seq, rat); rat2hs = posmap(rat, hs1)
    canon = {}   ## 2V8Q chain E residue number -> canonical human position
    for i, r in enumerate(res):
        p = rat2hs.get(e2rat.get(i + 1))
        if p: canon[r.id[1]] = p
    atoms = [a for r in res for a in r]; ns = NeighborSearch(atoms)
    official = [int(r["POS_PRKAG1"]) for r in rows(G + "/contact_columns_linsi.tsv")]
    sites = {}
    for lig in [r for r in E if r.get_resname() == "AMP"]:
        near = sorted({canon[a.get_parent().id[1]] for x in lig for a in ns.search(x.coord, 5.0) if a.get_parent().id[1] in canon})
        sites["AMP%d" % lig.id[1]] = {"canonical": near, "ligand_atoms": [[float(v) for v in x.coord] for x in lig if x.element != "H"],
                                     "overlap_official_site1": len(set(near) & set(official))}
    s1 = max(sites, key=lambda k: sites[k]["overlap_official_site1"])
    ## 2V8Q site Calpha coordinates keyed by canonical position (the superposition reference)
    ref_ca = {canon[r.id[1]]: [float(v) for v in r["CA"].coord] for r in res if r.id[1] in canon and "CA" in r}
    ## fixed box per site in the 2V8Q frame: site residues' heavy-atom extent +- PAD (pyKVFinder adds its probe margin)
    inv = {v: k for k, v in canon.items()}
    for k, v in sites.items():
        xyz = np.array([a.coord for c in v["canonical"] for a in E[inv[c]] if a.element != "H"])
        lo, hi = xyz.min(0) - PAD, xyz.max(0) + PAD
        v["box"] = [[float(lo[0]), float(lo[1]), float(lo[2])], [float(hi[0]), float(lo[1]), float(lo[2])], [float(lo[0]), float(hi[1]), float(lo[2])], [float(lo[0]), float(lo[1]), float(hi[2])]]
        open(O + "/box_%s.toml" % k, "w").write("[box]\n" + "".join("p%d = [%.3f, %.3f, %.3f]\n" % (i + 1, *pt) for i, pt in enumerate(v["box"])))
    ## each site's AMP as a ligand file in the 2V8Q frame (pyKVFinder ligand mode keeps cavity points within 5 A of it)
    for lig in [r for r in fetch("2V8Q")["E"] if r.get_resname() == "AMP"]:
        class One(Select):
            def accept_residue(self, r, _id=lig.id): return r.get_parent().id == "E" and r.id == _id
        io = PDBIO(); io.set_structure(fetch("2V8Q")); io.save(O + "/ligand_AMP%d.pdb" % lig.id[1], One())
    json.dump({"site1": s1, "sites": sites, "official_site1": official, "ref_ca": ref_ca, "hs1": hs1}, open(O + "/sites.json", "w"))
    for k, v in sites.items(): print(k, "site residues (canonical):", v["canonical"], "| overlap with the 10 official contacts:", v["overlap_official_site1"])
    print("site 1 =", s1)
    sys.exit(0)
SITES = json.load(open(O + "/sites.json")); HS1 = SITES["hs1"]; REFCA = {int(k): np.array(v) for k, v in SITES["ref_ca"].items()}
def model_files():
    fs = sorted(glob.glob(R + "/*_unrelaxed_rank_*_seed_*.pdb"))
    return fs
def kvf_volume(path, Rm, t, site, tag):
    ## model moved into the 2V8Q frame (x_ref = R^T (x_model - t)), B column set to 0, pyKVFinder in the fixed site box
    import pyKVFinder
    st = PDBParser(QUIET=True).get_structure("m", path)
    for a in st.get_atoms(): a.set_coord((a.coord - t) @ Rm); a.set_bfactor(0.0)
    wd = W + "/kvf"; os.makedirs(wd, exist_ok=True); f = wd + "/%s_%s.pdb" % (tag, site)
    io = PDBIO(); io.set_structure(st); io.save(f)
    r = pyKVFinder.run_workflow(f, box=O + "/box_%s.toml" % site, ligand=O + "/ligand_%s.pdb" % site, ligand_cutoff=5.0, nthreads=1)
    if r is None or not r.volume: return 0.0, 0, 0.0
    v = list(r.volume.values()); return round(float(sum(v)), 2), len(v), round(float(max(v)), 2)
def score_structure(path, query, model_id, plddt_override=None, do_volume=True, tag=None):
    st = PDBParser(QUIET=True).get_structure("m", path)[0]; ch = next(iter(st)); res = protein_residues(ch)
    seq = "".join(THREE[r.get_resname()] for r in res); m = posmap(HS1, seq)   ## canonical -> model index (1-based)
    heavy = np.array([a.coord for r in res for a in r if a.element != "H"]); out = []
    for site, info in SITES["sites"].items():
        cps = info["canonical"]; mapped = {c: res[m[c] - 1] for c in cps if c in m}
        pl = [plddt_override if plddt_override is not None else r["CA"].get_bfactor() for r in mapped.values() if "CA" in r]
        row = dict(query=query, model=model_id, site=site, site1=int(site == SITES["site1"]), site_residues=len(cps), mapped=len(mapped),
                   min_plddt=round(min(pl), 1) if pl else "", assessable=int(len(mapped) == len(cps) and pl and min(pl) >= 70))
        common = [c for c in cps if c in mapped and c in REFCA and "CA" in mapped[c]]
        if len(common) >= 3:
            Pm = np.array([REFCA[c] for c in common]); Qm = np.array([mapped[c]["CA"].coord for c in common]); Rm, t = kabsch(Pm, Qm)
            rmsd = float(np.sqrt(((Pm @ Rm.T + t - Qm) ** 2).sum(1).mean()))
            lig = np.array(info["ligand_atoms"]) @ Rm.T + t
            dmin = float(np.sqrt(((heavy[:, None, :] - lig[None, :, :]) ** 2).sum(2)).min())
            row.update(site_ca_rmsd=round(rmsd, 3), amp_min_heavy_dist=round(dmin, 3))
            if site == SITES["site1"]:
                for c in SITES["official_site1"]:
                    if c in m:
                        rr = res[m[c] - 1]; xyz = np.array([a.coord for a in rr if a.element != "H"])
                        row["c%d_aa" % c] = THREE[rr.get_resname()]; row["c%d_dist" % c] = round(float(np.sqrt(((xyz[:, None, :] - lig[None, :, :]) ** 2).sum(2)).min()), 3)
                    else: row["c%d_aa" % c] = "-"; row["c%d_dist" % c] = ""
        if do_volume and row["assessable"] and len(common) >= 3:
            row["kvf_volume"], row["kvf_cavities"], row["kvf_largest"] = kvf_volume(path, Rm, t, site, tag or re.sub(r"[^A-Za-z0-9_]", "", query + "_" + model_id))
        out.append(row)
    return out
def strip_ligands(pdbid, chain, dest):
    st = fetch(pdbid)
    class Sel(Select):
        def accept_chain(self, c): return c.id == chain
        def accept_residue(self, r): return r.id[0] == " " and r.get_resname() in THREE
    io = PDBIO(); io.set_structure(st); io.save(dest, Sel())
    return dest
def write(rowsout, f):
    keys = []
    for r in rowsout: keys += [k for k in r if k not in keys]
    with open(f, "w") as o:
        o.write("\t".join(keys) + "\n"); [o.write("\t".join(str(r.get(k, "")) for k in keys) + "\n") for r in rowsout]
if sys.argv[1] == "analyze":
    t, nt = int(sys.argv[2]), int(sys.argv[3]); fs = model_files()[t::nt]; out = []
    if t == 0:   ## structure controls: ligand-stripped crystal gamma chains (2V8Q rat gamma1 chain E; 4CFE human gamma1 chain E)
        for pid in ("2V8Q", "4CFE"):
            st = fetch(pid); gch = max(st, key=lambda c: AL.score("".join(THREE[r.get_resname()] for r in protein_residues(c)) or "X", HS1))
            f = strip_ligands(pid, gch.id, ST + "/%s_gamma_%s_noligand.pdb" % (pid, gch.id)); out += score_structure(f, "XTAL_" + pid, "chain_" + gch.id, plddt_override=100)
    for f in fs:
        b = os.path.basename(f); q = b.split("_unrelaxed_rank_")[0]; mid = re.search(r"model_\d+_seed_\d+", b).group(0)
        out += score_structure(f, q, mid)
    write(out, O + "/part_%02d.tsv" % t); print("task", t, len(fs), "models,", len(out), "rows")
if sys.argv[1] == "determinism":   ## the volume must be identical when recomputed (pyKVFinder has no random step)
    f = ST + "/2V8Q_gamma_E_noligand.pdb"; a = score_structure(f, "XTAL_2V8Q", "run1", plddt_override=100, tag="det1"); b = score_structure(f, "XTAL_2V8Q", "run2", plddt_override=100, tag="det2")
    same = all(x.get("kvf_volume") == y.get("kvf_volume") for x, y in zip(a, b))
    open(O + "/volume_determinism.txt", "w").write("2V8Q gamma, all sites, two runs: volumes identical = %s; %s\n" % (same, [(x["site"], x.get("kvf_volume")) for x in a]))
    print(open(O + "/volume_determinism.txt").read())
if sys.argv[1] == "summarize":
    import pandas as pd
    d = pd.concat([pd.read_csv(f, sep="\t") for f in sorted(glob.glob(O + "/part_*.tsv"))], ignore_index=True)
    lab = {r["SPECIES"]: r for r in rows(G + "/contact_retained_official.tsv")}
    d["species"] = d["query"].str.replace("^G_", "", regex=True)
    d["eroded_full"] = d["species"].map(lambda s: lab[s]["ERODED_FULL"] if s in lab else "")
    num = ["min_plddt", "site_ca_rmsd", "amp_min_heavy_dist", "kvf_volume", "kvf_cavities"]
    for c in num: d[c] = pd.to_numeric(d[c], errors="coerce")
    g = d.groupby(["query", "site"]).agg(models=("model", "size"), assessable=("assessable", "sum"), eroded_full=("eroded_full", "first"),
                                          min_plddt_median=("min_plddt", "median"), rmsd_median=("site_ca_rmsd", "median"),
                                          amp_dist_median=("amp_min_heavy_dist", "median"),
                                          volume_median_assessable=("kvf_volume", "median"), no_cavity_models=("kvf_cavities", lambda x: int((x == 0).sum()))).reset_index()
    ## volume and geometry only from assessable models
    ga = d[d.assessable == 1].groupby(["query", "site"]).agg(rmsd_median_assessable=("site_ca_rmsd", "median"), amp_dist_median_assessable=("amp_min_heavy_dist", "median")).reset_index()
    g = g.merge(ga, on=["query", "site"], how="left")
    with open(O + "/p45_colabfold_summary.tsv", "w") as o:
        o.write("## p45b: per query and site; gate = all site residues pLDDT >= 70; geometry and pyKVFinder volume (fixed 2V8Q-frame box) read only from assessable models; no_cavity_models = assessable models with no cavity in the box (volume 0)\n"); g.to_csv(o, sep="\t", index=False)
    s1 = g[g.site == SITES["site1"]]
    cmp = s1.assign(group=s1["query"].map(lambda q: "control" if q.startswith(("CTRL", "XTAL")) else ("eroded" if s1.set_index("query").loc[q, "eroded_full"] == "1" else "retained")))
    tab = cmp.groupby("group").agg(queries=("query", "size"), with_any_assessable=("assessable", lambda x: int((x > 0).sum())), rmsd=("rmsd_median_assessable", "median"),
                                   amp_dist=("amp_dist_median_assessable", "median"), volume=("volume_median_assessable", "median")).reset_index()
    with open(O + "/p45_colabfold_site1_groups.tsv", "w") as o:
        o.write("## p45b: site 1 by group (medians of per-query medians over assessable models)\n"); tab.to_csv(o, sep="\t", index=False)
    print(tab.to_string(index=False)); print(s1[s1["query"].str.startswith(("CTRL", "XTAL"))].to_string(index=False))
