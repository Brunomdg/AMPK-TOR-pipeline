# >>> deposit header
# Purpose: Builds the prediction sequence of each of the 236 AMPKγ copies so that every observed CBS domain is
#          kept (the full-length protein where no CBS domain was observed, since the domain boundaries are then
#          uncertain), and decides for each copy whether completed models of the identical sequence from the first
#          γ-alone run, or the MSA of the identical query from the first all-copy run, can be reused.
# Methods: Structure prediction and the pocket geometry of AMPKγ
# Inputs:  `DOMAIN_GUARDRAIL_AUDIT.tsv` (per-copy observed CBS domains and proposed spans);
#          `ALL_AMPKG_FULL_LENGTH.faa`; `all_copy_queries.fasta` of the first all-copy run; `queries.fasta` and the
#          completion markers of the first γ-alone run (`p45_colabfold`)
# Outputs: in `domain_preserving_v2`: `all_copy_queries.fasta`, `new_msa_queries.fasta`,
#          `STRUCTURAL_QUERY_MANIFEST.tsv` (one row per copy: span, length, SHA-256, reuse decision) and
#          `STRUCTURAL_PREPARATION.json` (counts); existing FASTA outputs are never overwritten
# Tools:   Python 3 standard library; versions: see envs/
# <<< deposit header
# deposit path defaults: begin
_TA210_SCRATCH = __import__("os").environ.get("TA210_SCRATCH") or "/gpfs2/scratch/bdigenov"
# deposit path defaults: end
from pathlib import Path
import csv,json,hashlib,collections
P=Path((_TA210_SCRATCH + '/ms210_phase2/TOR-AMPK_orthologs_10target_R1/20_port_68taxa_analyses'));old=P/'02_derived/all_ampkg_expansion_20260919';D=old/'domain_preserving_v2';D.mkdir(exist_ok=True)
def fasta(p):
    d={};key=None
    for l in p.read_text().splitlines():
        if l.startswith('>'):key=l[1:].split()[0];d[key]=''
        elif key:d[key]+=l.strip()
    return d
def read(p):return list(csv.DictReader(p.open(),delimiter='\t'))
rows=read(old/'DOMAIN_GUARDRAIL_AUDIT.tsv');full=fasta(old/'ALL_AMPKG_FULL_LENGTH.faa');previous=fasta(old/'all_copy_queries.fasta');original=fasta(P/'02_derived/p45_colabfold/queries.fasta')
reusable=collections.defaultdict(list)
for q,s in original.items():
    if (P/'03_results/p45_colabfold'/f'{q}.done.txt').exists():reusable[s].append(q)
allq={};newmsa={};out=[]
for r in rows:
    seq=full[r['COPY_ID']];q=r['QUERY']
    if r['STATUS']=='NO_OBSERVED_CBS_NOT_EVIDENCE_OF_ABSENCE':a,b=1,len(seq);mode='FULL_LENGTH_UNCERTAIN_DOMAIN_BOUNDARIES'
    else:a,b=int(r['PROPOSED_START']),int(r['PROPOSED_END']);mode='MSA_SPAN_EXPANDED_TO_RETAIN_OBSERVED_CBS'
    query=seq[a-1:b];assert query;allq[q]=query
    reuse=sorted(reusable.get(query,[]));reuse=reuse[0] if reuse else ''
    msa_reuse=(query==previous[q] and not reuse)
    if not reuse and not msa_reuse:newmsa[q]=query
    out.append({**r,'PREDICTION_LENGTH':len(query),'INPUT_MODE':mode,'SEQUENCE_SHA256':hashlib.sha256(query.encode()).hexdigest(),'FINAL_START':a,'FINAL_END':b,'REUSE_QUERY_EXACT_SEQUENCE':reuse,'REUSE_MSA_FROM_V1':int(msa_reuse),'ACTION':'REUSE_VERIFIED_SEQUENCE' if reuse else 'NEW_PREDICTION'})
assert len(out)==236
for filename,seqs in (('all_copy_queries.fasta',allq),('new_msa_queries.fasta',newmsa)):
    p=D/filename;assert not p.exists();p.write_text(''.join('>'+q+'\n'+s+'\n' for q,s in seqs.items()))
with (D/'STRUCTURAL_QUERY_MANIFEST.tsv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
summary={'status':'DOMAIN_PRESERVING_INPUTS_PREPARED','copies':236,'new_predictions_required':sum(r['ACTION']=='NEW_PREDICTION' for r in out),'reuse_completed_models':sum(r['ACTION']=='REUSE_VERIFIED_SEQUENCE' for r in out),'new_msas_required':len(newmsa),'reuse_exact_v1_msa':sum(r['REUSE_MSA_FROM_V1'] for r in out),'full_length_uncertain_boundaries':sum(r['INPUT_MODE']=='FULL_LENGTH_UNCERTAIN_DOMAIN_BOUNDARIES' for r in out),'max_length':max(len(s) for s in allq.values()),'toxoplasma':[r for r in out if r['SPECIES']=='Toxoplasma_gondii'],'exclusions':0,'guardrails':'Observed sequence-verified CBS envelopes retained; no human-identity cutoff, no contact-conservation cutoff, no annotation-absence exclusion. Full-length fallback is uncertainty handling, not a positive domain annotation.'}
(D/'STRUCTURAL_PREPARATION.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
