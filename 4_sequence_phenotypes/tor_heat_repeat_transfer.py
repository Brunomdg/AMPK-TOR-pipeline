# >>> deposit header
# Purpose: Transfers the 32 N-terminal HEAT repeats annotated in UniProt for human mTOR (P42345) to every TOR copy
#          through the full-length alignment, and reclassifies each copy by the HEAT signature matches that fall
#          inside the transferred region, with its coverage of and identity to human residues 16 to 1345.
# Methods: TOR N-terminal HEAT repeats
# Inputs:  `human_mTOR_P42345.txt` (UniProt entry, https://rest.uniprot.org/uniprotkb/P42345.txt);
#          `TOR_full_length_MSA.faa`; `TOR_all_retained.faa`; `TOR_HEAT_classification.tsv`;
#          `TOR_observed_domains.tsv`; `VALIDATION_AND_SUMMARY.json`
# Outputs: `TOR_HEAT_classification.tsv` (N-terminal columns added); `TOR_Nterminal_repeat_mapping.tsv` (one row
#          per copy and reference repeat); `TOR_Nterminal_HEAT_MSA.faa`; updated `VALIDATION_AND_SUMMARY.json`
#          and `SHA256SUMS.txt`
# Tools:   Python 3 standard library only; versions: see envs/
# <<< deposit header
# Ran as: refine_nterminal_heat.py
from pathlib import Path
import csv,json,re,hashlib,collections
D=Path(__file__).resolve().parent
def read(name):return list(csv.DictReader((D/name).open(),delimiter='\t'))
def fasta(name):
    out={}
    for l in (D/name).open():
        if l.startswith('>'):key=l[1:].split()[0];out[key]=''
        else:out[key]+=l.strip().upper()
    return out
def table(name,rows):
    with (D/name).open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
text=(D/'human_mTOR_P42345.txt').read_text()
human=re.sub('[^A-Z]','',text.split('SQ   SEQUENCE')[1].split('\n',1)[1])
repeats=[(int(a),int(b),n) for a,b,n in re.findall(r'FT   REPEAT\s+(\d+)\.\.(\d+)\nFT\s+/note="HEAT (\d+)"',text) if int(a)<1382]
assert len(repeats)==32
aln=fasta('TOR_full_length_MSA.faa');full=fasta('TOR_all_retained.faa');rows=read('TOR_HEAT_classification.tsv')
ref=[r['COPY_ID'] for r in rows if r['SPECIES']=='Homo_sapiens' and full[r['COPY_ID']]==human]
assert len(ref)==1;ref=ref[0]
hcols=[i for i,c in enumerate(aln[ref]) if c!='-']
start,end=hcols[min(a for a,b,n in repeats)-1],hcols[max(b for a,b,n in repeats)-1]
native=collections.defaultdict(list)
for r in read('TOR_observed_domains.tsv'):
    if r['HEAT_EXPLICIT']=='1':native[r['COPY_ID']].append((int(r['START']),int(r['END']),r['ACCESSION']))
mapped=[];ntaln={}
for r in rows:
    key=r['COPY_ID'];s=aln[key];cols=[i for i,c in enumerate(s) if c!='-'];aa={c:i+1 for i,c in enumerate(cols)}
    region=[c for c in cols if start<=c<=end];coverage=sum(s[c]!='-' for c in hcols[15:1345])/1330
    jointly=[c for c in hcols[15:1345] if s[c]!='-'];identity=sum(s[c]==aln[ref][c] for c in jointly)/len(jointly) if jointly else 0
    hits=[(a,b,acc) for a,b,acc in native[key] if any(start<=cols[p-1]<=end for p in range(a,b+1))]
    status='CURATED_NTERMINAL_HEAT_REFERENCE' if key==ref else 'NTERMINAL_HEAT_SIGNATURE_DETECTED' if hits else 'NTERMINAL_HEAT_NOT_DETECTED_UNRESOLVED'
    r.update({'NTERMINAL_HEAT_CLASS':status,'NTERMINAL_HEAT_NATIVE_SIGNATURES':len(hits),'NTERMINAL_HEAT_NATIVE_ACCESSIONS':';'.join(sorted({x[2] for x in hits})),'NTERMINAL_HEAT_MAPPED_START_AA':aa[region[0]] if region else '', 'NTERMINAL_HEAT_MAPPED_END_AA':aa[region[-1]] if region else '', 'NTERMINAL_HEAT_MAPPED_LENGTH_AA':len(region),'NTERMINAL_HEAT_REFERENCE_COVERAGE':round(coverage,4),'NTERMINAL_HEAT_REFERENCE_IDENTITY':round(identity,4),'NTERMINAL_HEAT_LENGTH_EVIDENCE':'CURATED_UNIPROT_P42345_16..1345' if key==ref else 'MSA_TRANSFER_FROM_P42345_16..1345_NOT_STRUCTURAL_VALIDATION'})
    ntaln[key]=s[start:end+1]
    for a,b,n in repeats:
        cs=[c for c in cols if hcols[a-1]<=c<=hcols[b-1]]
        rc=hcols[a-1:b]
        mapped.append({'COPY_ID':key,'SPECIES':r['SPECIES'],'PROTEIN_ID':r['PROTEIN_ID'],'REFERENCE_REPEAT':n,'REFERENCE_START':a,'REFERENCE_END':b,'MAPPED_START_AA':aa[cs[0]] if cs else '', 'MAPPED_END_AA':aa[cs[-1]] if cs else '', 'MAPPED_LENGTH_AA':len(cs),'REFERENCE_RESIDUES_COVERED':sum(s[c]!='-' for c in rc),'REFERENCE_REPEAT_LENGTH':b-a+1,'EVIDENCE':'CURATED_REFERENCE' if key==ref else 'ALIGNMENT_TRANSFER_NOT_DETECTED_REPEAT'})
table('TOR_HEAT_classification.tsv',rows);table('TOR_Nterminal_repeat_mapping.tsv',mapped)
(D/'TOR_Nterminal_HEAT_MSA.faa').write_text(''.join('>'+k+'\n'+s+'\n' for k,s in ntaln.items()))
summary=json.loads((D/'VALIDATION_AND_SUMMARY.json').read_text())
summary.update({'Nterminal_HEAT_classes':dict(collections.Counter(r['NTERMINAL_HEAT_CLASS'] for r in rows)),'Nterminal_reference':{'accession':'P42345','copy_id':ref,'identical_sequence_verified':True,'length':len(human),'repeat_count':len(repeats),'envelope_aa':[16,1345],'source_url':'https://rest.uniprot.org/uniprotkb/P42345.txt','sha256':hashlib.sha256((D/'human_mTOR_P42345.txt').read_bytes()).hexdigest()},'Toxoplasma':[r for r in rows if r['SPECIES']=='Toxoplasma_gondii']})
(D/'VALIDATION_AND_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
(D/'SHA256SUMS.txt').write_text(''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.name+'\n' for p in sorted(D.iterdir()) if p.is_file() and p.suffix in ('.faa','.tsv','.json','.py','.sbatch','.txt') and p.name!='SHA256SUMS.txt'))
print(json.dumps(summary,indent=2))
