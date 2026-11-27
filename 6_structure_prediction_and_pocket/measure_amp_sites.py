# >>> deposit header
# Purpose: Measures the three reference AMP sites in all 25 models of each γ-alone prediction of the 236 AMPKγ
#          copies and the five controls (241 queries), one sixteenth of the queries per array task or two queries
#          in pilot mode, with no confidence filter: site pLDDT, superposition RMSD, distance to the transferred
#          AMP, site-1 contact distances and pyKVFinder pocket volume.
# Methods: Structure prediction and the pocket geometry of AMPKγ
# Inputs:  `STRUCTURAL_QUERY_MANIFEST.tsv` (the table of the 236 copy queries); the ColabFold models in
#          `all_ampkg_domain_preserving_20260919` and, for reused sequences and the controls, in `p45_colabfold`;
#          `queries.fasta` of the first γ-alone run; `sites.json`; the functions of `p45b_colabfold_analysis.py`
# Outputs: `part_NN.tsv` (or `pilot.tsv`), one row per model and site; `part_NN_PROVENANCE.json` (checksums of
#          the model files, the function library and `sites.json`, with row and error counts)
# Tools:   Python 3 with numpy; Biopython and pyKVFinder through `p45b_colabfold_analysis.py` (pyKVFinder 0.9.5
#          is held only in the `p68_kvf` environment [env record]); versions: see envs/
# <<< deposit header
# Ran as: analyze_all_ampkg_structures.py
# deposit path defaults: begin
_TA210_SCRATCH = __import__("os").environ.get("TA210_SCRATCH") or "/gpfs2/scratch/bdigenov"
# deposit path defaults: end
from pathlib import Path
import csv, json, os, sys, re, hashlib, importlib.util
import numpy as np
S=Path((_TA210_SCRATCH + '/ms210_phase2/TOR-AMPK_orthologs_10target_R1'));P=S/'20_port_68taxa_analyses'
D=P/'02_derived/all_ampkg_expansion_20260919/domain_preserving_v2'
O=P/'03_results/all_ampkg_structural_evidence_20260919';O.mkdir(exist_ok=True)
mode=sys.argv[1] if len(sys.argv)>1 else 'array'
sys.argv=[sys.argv[0],'library_only']
source=Path(os.path.dirname(os.path.abspath(__file__)))/'lib'/'amp_site_measurement.py'
spec=importlib.util.spec_from_file_location('frozen_pocket_analysis',source);lib=importlib.util.module_from_spec(spec);spec.loader.exec_module(lib)
lib.W=str(O/'volume_work')
def rd(path):return list(csv.DictReader((l for l in path.open() if not l.startswith('#')),delimiter='\t'))
manifest=rd(D/'STRUCTURAL_QUERY_MANIFEST.tsv');assert len(manifest)==236
queries=[]
for row in manifest:
    reuse=row['ACTION']=='REUSE_VERIFIED_SEQUENCE'
    queries.append({'query':row['QUERY'],'source_query':row['REUSE_QUERY_EXACT_SEQUENCE'] if reuse else row['QUERY'],
                    'directory':P/'03_results'/('p45_colabfold' if reuse else 'all_ampkg_domain_preserving_20260919'),
                    'species':row['SPECIES'],'protein':row['PROTEIN_ID'],'expected_sha':row['SEQUENCE_SHA256'],
                    'expected_length':int(row['PREDICTION_LENGTH']),'kind':'RETAINED_COPY','reuse':reuse})
old_queries={};k=None
for line in (P/'02_derived/p45_colabfold/queries.fasta').read_text().splitlines():
    if line.startswith('>'):k=line[1:].split()[0];old_queries[k]=''
    elif k:old_queries[k]+=line.strip()
for q,seq in old_queries.items():
    if q.startswith('CTRL'):
        queries.append({'query':q,'source_query':q,'directory':P/'03_results/p45_colabfold','species':'CONTROL','protein':'',
                        'expected_sha':hashlib.sha256(seq.encode()).hexdigest(),'expected_length':len(seq),'kind':'PREDICTED_CONTROL','reuse':True})
assert len(queries)==241
task=int(os.environ.get('SLURM_ARRAY_TASK_ID','0'))
selected=([next(q for q in queries if q['kind']=='RETAINED_COPY' and q['reuse']),next(q for q in queries if q['kind']=='PREDICTED_CONTROL')]
          if mode=='pilot' else queries[task::16])
records=[];inputs=[]
for query in selected:
    files=sorted(query['directory'].glob(query['source_query']+'_unrelaxed_rank_*_seed_*.pdb'))
    assert len(files)==25,(query['query'],len(files))
    expected={f'model_{m}_seed_{seed:03d}' for m in range(1,6) for seed in range(5)}
    observed={re.search(r'model_\d+_seed_\d+',p.name).group() for p in files};assert observed==expected,(query['query'],observed)
    for path in files:
        model=re.search(r'model_\d+_seed_\d+',path.name).group()
        structure=lib.PDBParser(QUIET=True).get_structure('model',str(path))[0]
        residues=lib.protein_residues(next(iter(structure)))
        sequence=''.join(lib.THREE[r.get_resname()] for r in residues)
        assert len(sequence)==query['expected_length'],(query['query'],'length')
        assert hashlib.sha256(sequence.encode()).hexdigest()==query['expected_sha'],(query['query'],'sequence mismatch')
        mapping=lib.posmap(lib.HS1,sequence)
        rows=lib.score_structure(str(path),query['query'],model,do_volume=False)
        for row in rows:
            site=row['site'];info=lib.SITES['sites'][site]
            mapped={c:residues[mapping[c]-1] for c in info['canonical'] if c in mapping}
            values=[r['CA'].get_bfactor() for r in mapped.values() if 'CA' in r]
            common=[c for c in mapped if c in lib.REFCA and 'CA' in mapped[c]]
            row.update(SPECIES=query['species'],PROTEIN_ID=query['protein'],INPUT_KIND=query['kind'],
                       SOURCE_QUERY=query['source_query'],CONFIDENCE_FILTER_APPLIED=0,
                       mean_site_plddt=float(np.mean(values)) if values else '',
                       median_site_plddt=float(np.median(values)) if values else '',
                       n_site_plddt_below70=sum(v<70 for v in values),
                       geometry_status='AVAILABLE' if len(common)>=3 else 'INSUFFICIENT_MAPPED_REFERENCE_ATOMS')
            if len(common)>=3:
                R,t=lib.kabsch(np.array([lib.REFCA[c] for c in common]),np.array([mapped[c]['CA'].coord for c in common]))
                try:
                    row['kvf_volume'],row['kvf_cavities'],row['kvf_largest']=lib.kvf_volume(str(path),R,t,site,query['query']+'_'+model)
                    row['volume_status']='CALCULATED_WITHOUT_CONFIDENCE_FILTER'
                except Exception as error:
                    row['volume_status']='CALCULATION_ERROR';row['volume_error']=repr(error)
            else:row['volume_status']='UNAVAILABLE_MAPPING'
            records.append(row)
        inputs.append({'QUERY':query['query'],'MODEL':model,'PDB':str(path),'SHA256':hashlib.sha256(path.read_bytes()).hexdigest()})
    print(query['query'],'25 models retained',flush=True)
name='pilot' if mode=='pilot' else f'part_{task:02d}'
lib.write(records,str(O/(name+'.tsv')))
receipt={'status':'CALCULATED_REQUIRES_VALIDATION','queries':len(selected),'models':len(inputs),'site_rows':len(records),
         'confidence_filter_applied':False,'volume_errors':sum(r.get('volume_status')=='CALCULATION_ERROR' for r in records),
         'library_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
         'sites_sha256':hashlib.sha256((P/'03_results/p45_analysis/sites.json').read_bytes()).hexdigest(),'inputs':inputs}
(O/(name+'_PROVENANCE.json')).write_text(json.dumps(receipt,indent=2)+'\n')
print(name,len(records),'site rows; volume errors',receipt['volume_errors'])
