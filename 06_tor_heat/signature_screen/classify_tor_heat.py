# >>> deposit header
# Purpose: Classifies each TOR copy by whether its HEAT signature matches lie N-terminal to the FAT domain, using the
#          copy's own FAT (PF02259) start or, without one, the median FAT start column of the full-length alignment
#          (taken from copies of at least 1,800 residues whose FAT domain starts at residue 500 or later).
# Methods: TOR N-terminal HEAT repeats
# Inputs:  `TOR_manifest.tsv`, `TOR_all_retained.faa`, `TOR_full_length_MSA.faa`, `TOR_observed_domains.tsv`
# Outputs: `TOR_HEAT_classification.tsv`; `TOR_FAT_boundary_donors.tsv` (copies that set the median boundary);
#          `TOR_domains_in_MSA.tsv`; `TOR_Nterminal_region_MSA.faa` (alignment cut at the boundary);
#          `TOR_MSA_column_occupancy.tsv`; `VALIDATION_AND_SUMMARY.json`; `SHA256SUMS.txt` (checksums of the
#          folder's sequence, table, JSON and script files)
# Tools:   Python 3 standard library only; versions: see envs/
# <<< deposit header
from pathlib import Path
import csv,json,hashlib,collections,statistics
D=Path(__file__).resolve().parent
def rd(name):return list(csv.DictReader((D/name).open(),delimiter='\t'))
def fasta(name):
    out={};key=None
    for l in (D/name).open():
        if l.startswith('>'):
            key=l[1:].split()[0];assert key not in out;out[key]=''
        elif key:out[key]+=l.strip().upper()
    return out
def table(name,rows):
    with (D/name).open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def union(intervals):
    result=[]
    for a,b in sorted(intervals):
        if result and a<=result[-1][1]+1:result[-1][1]=max(b,result[-1][1])
        else:result.append([a,b])
    return result
manifest=rd('TOR_manifest.tsv');full=fasta('TOR_all_retained.faa');aln=fasta('TOR_full_length_MSA.faa')
assert set(full)==set(aln)=={r['COPY_ID'] for r in manifest}
assert len({len(v) for v in aln.values()})==1
for r in manifest:
    k=r['COPY_ID'];assert aln[k].replace('-','')==full[k]
    assert hashlib.sha256(full[k].encode()).hexdigest()==r['SHA256']
domains=collections.defaultdict(list)
for r in rd('TOR_observed_domains.tsv'):
    r['START']=int(r['START']);r['END']=int(r['END']);domains[r['COPY_ID']].append(r)
cols={k:[i for i,c in enumerate(v) if c!='-'] for k,v in aln.items()}
fat={k:min(h['START'] for h in hs if h['FAT_ANCHOR']=='1') for k,hs in domains.items() if any(h['FAT_ANCHOR']=='1' for h in hs)}
donors=[(k,cols[k][a-1]) for k,a in fat.items() if a>=500 and len(full[k])>=1800]
assert len(donors)>=20
boundary=int(statistics.median(c for k,c in donors));q1,_,q3=statistics.quantiles([c for k,c in donors],n=4)
table('TOR_FAT_boundary_donors.tsv',[{'COPY_ID':k,'FAT_START_AA':fat[k],'FAT_START_MSA_COLUMN_1BASED':c+1} for k,c in donors])
rows=[];cropped={};domainmapped=[]
for r in manifest:
    k=r['COPY_ID'];seq=full[k];hs=domains[k]
    heat=[h for h in hs if h['HEAT_EXPLICIT']=='1'];intervals=union([(h['START'],h['END']) for h in heat])
    estimate=sum(c<boundary for c in cols[k]);low=sum(c<int(q1) for c in cols[k]);high=sum(c<int(q3) for c in cols[k])
    nend=fat[k]-1 if k in fat else estimate
    nt=union([(h['START'],min(h['END'],nend)) for h in heat if h['START']<=nend])
    status='HEAT_DETECTED_N_TERMINAL' if nt else 'HEAT_DETECTED_OUTSIDE_ESTIMATED_N_TERMINUS' if heat else 'HEAT_NOT_DETECTED_UNRESOLVED'
    cropped[k]=aln[k][:boundary]
    rows.append({'COPY_ID':k,'SPECIES':r['SPECIES'],'PROTEIN_ID':r['PROTEIN_ID'],'DISPOSITION':r['DISPOSITION'],'PROTEIN_LENGTH_AA':len(seq),'HEAT_CLASS':status,'HEAT_SIGNATURE_MATCHES':len(heat),'HEAT_INTERVALS_AA':';'.join(f'{a}..{b}' for a,b in intervals),'HEAT_UNION_LENGTH_AA':sum(b-a+1 for a,b in intervals) if heat else '', 'HEAT_ENVELOPE_START_AA':intervals[0][0] if heat else '', 'HEAT_ENVELOPE_END_AA':intervals[-1][1] if heat else '', 'HEAT_ENVELOPE_LENGTH_AA':intervals[-1][1]-intervals[0][0]+1 if heat else '', 'NTERM_HEAT_UNION_LENGTH_AA':sum(b-a+1 for a,b in nt) if nt else '', 'NTERM_REGION_LENGTH_PROXY_AA':nend,'NTERM_BOUNDARY_METHOD':'OWN_OBSERVED_PF02259_FAT_START_MINUS_1' if k in fat else 'MULTISPECIES_MSA_MEDIAN_FAT_BOUNDARY_TRANSFER','MSA_NTERM_LENGTH_PROXY_AA':estimate,'MSA_DONOR_Q1_LENGTH_PROXY_AA':low,'MSA_DONOR_Q3_LENGTH_PROXY_AA':high,'OBSERVED_FAT_START_AA':fat.get(k,''),'OBSERVED_FATC':int(any(h['ACCESSION']=='PF02260' for h in hs)),'SHORT_PROTEIN_REVIEW_LT1800AA':int(len(seq)<1800),'UNKNOWN_RESIDUES':sum(seq.count(x) for x in 'BXZJUO'),'SEQUENCE_SHA256':r['SHA256']})
    for h in hs:
        domainmapped.append({**h,'MSA_START_COLUMN_1BASED':cols[k][h['START']-1]+1,'MSA_END_COLUMN_1BASED':cols[k][h['END']-1]+1})
table('TOR_HEAT_classification.tsv',rows);table('TOR_domains_in_MSA.tsv',domainmapped)
(D/'TOR_Nterminal_region_MSA.faa').write_text(''.join('>'+k+'\n'+v+'\n' for k,v in cropped.items()))
occupancy=[]
for i in range(len(next(iter(aln.values())))):
    n=sum(v[i]!='-' for v in aln.values())
    occupancy.append({'MSA_COLUMN_1BASED':i+1,'NONGAP_SEQUENCES':n,'OCCUPANCY':round(n/len(aln),5)})
table('TOR_MSA_column_occupancy.tsv',occupancy)
summary={'status':'PASS','copies':len(rows),'species':len({r['SPECIES'] for r in rows}),'MSA_columns':len(next(iter(aln.values()))),'all_sequences_preserved_exactly':True,'HEAT_classes':dict(collections.Counter(r['HEAT_CLASS'] for r in rows)),'short_protein_review':sum(r['SHORT_PROTEIN_REVIEW_LT1800AA'] for r in rows),'FAT_boundary_donors':len(donors),'median_FAT_column_1based':boundary+1,'donor_FAT_column_Q1_Q3_1based':[q1+1,q3+1],'Toxoplasma':[r for r in rows if r['SPECIES']=='Toxoplasma_gondii'],'no_HEAT_detection_is_not_absence':True}
(D/'VALIDATION_AND_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
(D/'SHA256SUMS.txt').write_text(''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.name+'\n' for p in sorted(D.iterdir()) if p.is_file() and p.suffix in ('.faa','.tsv','.json','.py','.sbatch')))
print(json.dumps(summary,indent=2))
