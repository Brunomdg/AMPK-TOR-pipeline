# >>> deposit header
# Purpose: Crops each TOR AlphaFold model to its N-terminal HEAT region (transferred from human mTOR through the
#          alignment, ended before the observed FAT domain), aligns the crop with US-align to the control structures
#          other than the same-species TOR structure, and calls the region structure-supported when a TOR HEAT
#          reference passes the geometric and pLDDT rules and no non-HEAT control passes the geometric rule.
# Methods: TOR N-terminal HEAT repeats
# Inputs:  arguments: "pilot" (the three predicted controls TOR0099, TOR0182, TOR0183) or a batch word with task index
#          and number of tasks; `TOR_HEAT_classification.tsv`; `CONTROL_MANIFEST_EXTENDED.tsv`; `receipts/`;
#          the two control test summaries, which must show all pairs correct; `PRODUCTION_CROP_PREFLIGHT.json`
#          (batch runs need the pilot to have passed)
# Outputs: `direct_alignment_v1/<copy>/` (`query.pdb`, US-align reports, residue mappings, `SUMMARY.json`);
#          `direct_alignment_v1/PRODUCTION_CROP_PREFLIGHT.json` (pilot)
# Tools:   US-align (-mol prot); Biopython; versions: see envs/
# <<< deposit header
from pathlib import Path
import csv,json,subprocess,sys,hashlib,re,statistics
from Bio.PDB import PDBParser,PDBIO,Structure,Model,Chain
from Bio.SeqUtils import seq1
D=Path(__file__).resolve().parent;O=D/'direct_alignment_v1';O.mkdir(exist_ok=True)
pilot=sys.argv[1]=='pilot';idx=int(sys.argv[2]) if not pilot else 0;chunks=int(sys.argv[3]) if not pilot else 1
prior={r['COPY_ID']:r for r in csv.DictReader((D.parent/'tor_heat_20260919/TOR_HEAT_classification.tsv').open(),delimiter='\t')}
refs=list(csv.DictReader((D/'CONTROL_MANIFEST_EXTENDED.tsv').open(),delimiter='\t'))
assert json.loads((D/'control_validation_v2/VALIDATION_SUMMARY.json').read_text())['all_correct']
assert json.loads((D/'predicted_control_validation/SUMMARY.json').read_text())['all_correct']
if not pilot:assert json.loads((O/'PRODUCTION_CROP_PREFLIGHT.json').read_text())['passed']
inputs=[json.loads(p.read_text()) for p in sorted((D/'receipts').glob('TOR*.json'))]
inputs=[r for r in inputs if r['MODEL_FILES']]
if pilot:inputs=[r for r in inputs if r['COPY_ID'] in ('TOR0099','TOR0182','TOR0183')]
inputs=[r for i,r in enumerate(inputs) if i%chunks==idx]
def residues(path):
    obj=PDBParser(QUIET=True).get_structure('x',path);ch=max(obj[0],key=lambda c:sum('CA' in r for r in c))
    return [r for r in ch if r.id[0]==' ' and 'CA' in r]
def write_tsv(path,rows,fields=None):
    with path.open('w') as f:
        w=csv.DictWriter(f,fieldnames=fields or list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
preflight=[]
for record in inputs:
    key=record['COPY_ID'];r=prior[key];out=O/key;out.mkdir(exist_ok=True)
    if (out/'SUMMARY.json').exists() and not pilot:continue
    m=record['MODEL_FILES'][0];path=Path(m['path']);rr=residues(path);seq=''.join(seq1(x.resname) for x in rr)
    assert seq==m['metadata']['sequence'] and hashlib.sha256(seq.encode()).hexdigest()==r['SEQUENCE_SHA256'],key
    a=int(r['NTERMINAL_HEAT_MAPPED_START_AA'] or 0);b=int(r['NTERMINAL_HEAT_MAPPED_END_AA'] or 0)
    fat=int(r['OBSERVED_FAT_START_AA'] or 0)
    if fat:b=min(b,fat-1)
    if not (1<=a<=b<=len(rr)):
        (out/'SUMMARY.json').write_text(json.dumps({'COPY_ID':key,'status':'UNRESOLVED_NO_USABLE_REGION','pairs':[]},indent=2));continue
    cropped=rr[a-1:b];query_seq=seq[a-1:b]
    obj=Structure.Structure(key);model=Model.Model(0);obj.add(model);ch=Chain.Chain('A');model.add(ch)
    for x in cropped:ch.add(x.copy())
    qp=out/'query.pdb';io=PDBIO();io.set_structure(obj);io.save(str(qp))
    results=[]
    for ref in refs:
        if r['SPECIES']=='Homo_sapiens' and ref['CONTROL']=='HUMAN_TOR_HEAT':continue
        if r['SPECIES']=='Saccharomyces_cerevisiae' and ref['CONTROL']=='SCEREVISIAE_TOR2_HEAT':continue
        k=ref['CONTROL'];report=out/(k+'.usalign.txt');cmd=['USalign',str(qp),ref['STRUCTURE'],'-mol','prot']
        if not report.exists():
            tmp=report.with_suffix('.partial')
            with tmp.open('w') as f:run=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,timeout=900)
            assert run.returncode==0,(key,k);tmp.replace(report)
        txt=report.read_text();tm=[float(v) for v in re.findall(r'TM-score=\s*([0-9.]+) \(normalized by length of Structure_[12]',txt)]
        assert len(tm)==2,(key,k,'scores')
        lengths=[int(v) for v in re.findall(r'Length of Structure_[12]:\s*(\d+) residues',txt)];assert len(lengths)==2
        ali=int(re.search(r'Aligned length=\s*(\d+)',txt).group(1));lines=txt.splitlines();start=next(i for i,l in enumerate(lines) if 'denotes residue pairs' in l)
        qa,labels,ta=lines[start+1:start+4];assert len(qa)==len(labels)==len(ta),(key,k,'alignment lengths')
        tr=residues(Path(ref['STRUCTURE']));assert qa.replace('-','')==query_seq and ta.replace('-','')==''.join(seq1(x.resname) for x in tr),(key,k,'alignment sequence')
        qi=ti=0;paired=[]
        for x,label,y in zip(qa,labels,ta):
            if x!='-':qi+=1
            if y!='-':ti+=1
            if x!='-' and y!='-' and label in ':.':
                paired.append({'QUERY_RESIDUE':a+qi-1,'REFERENCE_RESIDUE':tr[ti-1].id[1],'GEOMETRY':'LT5A' if label==':' else 'OTHER_ALIGNED','PLDDT':cropped[qi-1]['CA'].bfactor})
        assert len(paired)==ali,(key,k,len(paired),ali)
        write_tsv(out/(k+'_RESIDUE_MAPPING.tsv'),paired,['QUERY_RESIDUE','REFERENCE_RESIDUE','GEOMETRY','PLDDT'])
        close=[v for v in paired if v['GEOMETRY']=='LT5A'];cf=[v['PLDDT'] for v in close]
        geometry=min(tm)>=0.5 and min(ali/lengths[0],ali/lengths[1])>=0.5
        confidence=len(cf)>=100 and statistics.median(cf)>=70 and sum(v>=70 for v in cf)/len(cf)>=0.7
        results.append({'reference':k,'expected_family':ref['EXPECTED_CLASS'],'tm_query':tm[0],'tm_reference':tm[1],'aligned_residues':ali,'query_coverage':ali/lengths[0],'reference_coverage':ali/lengths[1],'geometric_match':geometry,'close_residues':len(cf),'close_median_plddt':statistics.median(cf) if cf else None,'close_fraction_plddt_ge70':sum(v>=70 for v in cf)/len(cf) if cf else None,'confidence_pass':confidence,'supported_match':geometry and confidence})
    positives=[v for v in results if v['expected_family']=='HEAT' and v['supported_match']]
    negative_matches=[v for v in results if v['expected_family']!='HEAT' and v['geometric_match']]
    status='STRUCTURE_SUPPORTED_NTERMINAL_TOR_HEAT_ARCHITECTURE' if positives and not negative_matches else ('UNRESOLVED_CONFLICTING_ARCHITECTURE' if negative_matches else 'UNRESOLVED_ALIGNMENT_OR_CONFIDENCE')
    positions=set();close_positions=set()
    for v in positives:
        for p in csv.DictReader((out/(v['reference']+'_RESIDUE_MAPPING.tsv')).open(),delimiter='\t'):
            positions.add(int(p['QUERY_RESIDUE']))
            if p['GEOMETRY']=='LT5A' and float(p['PLDDT'])>=70:close_positions.add(int(p['QUERY_RESIDUE']))
    summary={'COPY_ID':key,'status':status,'sequence_verified':True,'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'crop_start':a,'crop_end':b,'crop_length':b-a+1,'crop_method':'Existing human-HEAT MSA coordinate transfer, capped before observed FAT; region selection only, not proof of HEAT','supported_references':[v['reference'] for v in positives],'aligned_union_aa':len(positions),'aligned_envelope_start':min(positions) if positions else None,'aligned_envelope_end':max(positions) if positions else None,'aligned_envelope_length':max(positions)-min(positions)+1 if positions else None,'close_confident_union_aa':len(close_positions),'complete_heat_domain_length_aa':None,'pairs':results}
    (out/'SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(key,status,summary['aligned_union_aa'],flush=True)
    if pilot:preflight.append({'key':key,'passed':bool(positives) and not negative_matches,'positive_matches':len(positives),'negative_matches':len(negative_matches)})
if pilot:
    assert len(preflight)==3
    (O/'PRODUCTION_CROP_PREFLIGHT.json').write_text(json.dumps({'passed':all(r['passed'] for r in preflight),'controls':preflight,'confidence_rule':'At least 100 close (<5A) pairs, median pLDDT>=70 and >=70% of close pairs pLDDT>=70','tm_rule':'Both normalized TM>=0.5, both coverage>=0.5'},indent=2)+'\n')
    assert all(r['passed'] for r in preflight),preflight
