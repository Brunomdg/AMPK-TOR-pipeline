# >>> deposit header
# Purpose: Collects per-copy AMPKγ measurements (CBS-domain identity away from the AMP contacts, site-1 contact score,
#          agreement of pairwise and joint alignments at site 1, CBS domain evidence, predicted-structure confidence
#          and RMSD, length) and compares the coccidian copies with the others; descriptive statistics only.
# Methods: Sequence phenotypes of the retained subunits
# Inputs:  a JSON configuration file (first argument) naming the output folder, the number of copies and the input
#          tables: per-copy AMPKγ table, per-position residues, InterPro domain coordinates, per-site structure
#          summaries, the check of the residue correspondence, taxonomy, species states, agreement with the earlier
#          residue calls, and (after the rescue) the correspondence check for the two rescued Eimeria copies
# Outputs: `METRICS_PER_COPY.tsv`, `LOSS_PATTERN_APICOMPLEXA.tsv`, `CONCORDANCE_PER_COPY.tsv`, `GROUP_COMPARISONS.tsv`,
#          `ANALYSIS_SUMMARY.json` (with the SHA-256 of every input); measurements summarised in Supplementary
#          Table S28
# Tools:   Python 3 with pandas, NumPy and SciPy (scipy.stats.mannwhitneyu); versions: see envs/
# <<< deposit header
"""Descriptive measurements of the coccidian AMPKγ copies: extraction from existing tables and summary statistics.

No new alignment, refit or model. All inputs are existing tables, and the SHA-256 checksum of each is written to
ANALYSIS_SUMMARY.json. The script was run once on the published inputs and once on the inputs after the rescue
(see coccidian.sbatch). Outputs: METRICS_PER_COPY.tsv, LOSS_PATTERN_APICOMPLEXA.tsv, CONCORDANCE_PER_COPY.tsv,
GROUP_COMPARISONS.tsv (used in the report) and ANALYSIS_SUMMARY.json.
"""
import json, hashlib
from pathlib import Path
import pandas as pd
import numpy as np
from scipy import stats

import sys
CFG = json.loads(Path(sys.argv[1]).read_text())
OUT = Path(CFG['out']); OUT.mkdir(parents=True, exist_ok=True)
TS = CFG['ts']
N_COPIES = int(CFG['n_copies']); N_POS = N_COPIES * 60
COHORT_EXTRA = dict(CFG.get('cohort_extra', {}))
SRC = {k: Path(v) for k, v in CFG['inputs'].items()}
assert set(SRC) - {'register_extension'} == {'copy_analysis', 'position_source', 'domain_coords', 'site_summaries', 'phase1_validation', 'taxonomy', 'presence_long', 'legacy_reproduction'}

def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

INPUT_HASHES = {k: sha256(v) for k, v in SRC.items()}

# ---- Read the inputs ----
copies = pd.read_csv(SRC['copy_analysis'], sep='\t')
pos = pd.read_csv(SRC['position_source'], sep='\t')
dom = pd.read_csv(SRC['domain_coords'], sep='\t')
sites = pd.read_csv(SRC['site_summaries'], sep='\t')
val = pd.read_csv(SRC['phase1_validation'], sep='\t')
tax = pd.read_csv(SRC['taxonomy'], sep='\t', dtype=str)
t2 = pd.read_csv(SRC['presence_long'], sep='\t', comment='#')
leg = pd.read_csv(SRC['legacy_reproduction'], sep='\t', comment='#')

assert len(copies) == N_COPIES, len(copies)
assert len(pos) == N_POS, len(pos)
assert len(tax) == 210

# ---- Mark the coccidian and the apicomplexan copies ----
COCCIDIAN = {  # the coccidian set, chosen in advance: 6 eucoccidia (Sarcocystidae) and 3 Cryptosporidium
    'Toxoplasma_gondii': 'AG0212', 'Neospora_caninum': 'AG0130', 'Hammondia_hammondi': 'AG0084',
    'Sarcocystis_neurona': None, 'Besnoitia_besnoiti': None, 'Cystoisospora_suis': None,
    'Cryptosporidium_hominis': None, 'Cryptosporidium_muris': None, 'Cryptosporidium_parvum': None,
}
COCCIDIAN.update(COHORT_EXTRA)  # from 23 September 2026 the set also holds the two rescued Eimeria copies (configuration file)
api_species = set(tax.loc[tax['OFFICIAL_LINEAGE_WITH_RANKS'].str.contains('Apicomplexa', na=False), 'SPECIES_KEY'])
assert len(api_species) == 41
copies['APICOMPLEXAN'] = copies['SPECIES'].isin(api_species)
copies['COCCIDIAN'] = copies['SPECIES'].isin(COCCIDIAN)
assert copies['COCCIDIAN'].sum() == len(COCCIDIAN), copies['COCCIDIAN'].sum()
# Each coccidian species has exactly one copy
assert copies.loc[copies['COCCIDIAN'], 'SPECIES'].value_counts().max() == 1

# ---- Measurement 2: per-position agreement of the pairwise and joint alignments at site 1 (per-position table) ----
s1 = pos[pos['SITESET'] == 'site1']
piv = s1.pivot_table(index=['COPY_ID', 'CANONICAL_POS_PRKAG1'], columns='ALIGNMENT_METHOD',
                     values=['OBSERVED_AA', 'EXACT_MATCH'], aggfunc='first')
assert set(piv.columns.get_level_values(1).unique()) >= {'PAIRWISE', 'JOINT'}
agree = piv.groupby('COPY_ID').apply(
    lambda d: pd.Series({
        'SITE1_N_POSITIONS': len(d),
        'SITE1_PW_JOINT_SAME_OBSERVED_AA': int((d[('OBSERVED_AA', 'PAIRWISE')] == d[('OBSERVED_AA', 'JOINT')]).sum()),
        'SITE1_EXACT_PAIRWISE': int(pd.to_numeric(d[('EXACT_MATCH', 'PAIRWISE')], errors='coerce').fillna(0).sum()),
        'SITE1_EXACT_JOINT': int(pd.to_numeric(d[('EXACT_MATCH', 'JOINT')], errors='coerce').fillna(0).sum()),
    }), include_groups=False)
copies['COPY_ID'] = copies['SPECIES'] + '|' + copies['PROTEIN_ID']
copies = copies.merge(agree, on='COPY_ID', how='left')

# ---- Measurement 3: number of InterPro hits that are evidence of a CBS domain (covered length and fraction are
# already in the copy table) ----
cbs_hits = dom[dom['CBS_EVIDENCE'] == 1].groupby('QUERY').size().rename('N_CBS_EVIDENCE_HITS')
copies = copies.merge(cbs_hits, on='QUERY', how='left')
copies['N_CBS_EVIDENCE_HITS'] = copies['N_CBS_EVIDENCE_HITS'].fillna(0).astype(int)

# ---- Measurement 4: predicted structure per copy (median over its three sites) ----
stru = sites.groupby('QUERY').agg(
    STRUCT_N_SITES=('SITE', 'size'),
    MIN_PLDDT_MEDIAN_PERCOPY=('MIN_PLDDT_MEDIAN', 'median'),
    RMSD_MEDIAN_PERCOPY=('RMSD_MEDIAN', 'median'),
    ORIGINAL_GATE_PASS_MIN=('ORIGINAL_GATE_PASS', 'min'),
    ORIGINAL_GATE_PASS_SITES=('ORIGINAL_GATE_PASS', lambda s: int((s > 0).sum())),
)
copies = copies.merge(stru, on='QUERY', how='left')

# ---- Agreement of each copy with the earlier residue calls (input `legacy_reproduction`), in four categories ----
def leg_cat(r):
    if r['PAIRWISE_MATCHES_LEGACY'] == 'yes' and r['JOINT_MATCHES_LEGACY'] == 'yes':
        return 'AGREE_BOTH'
    if r['PAIRWISE_MATCHES_LEGACY'] == 'yes':
        return 'AGREE_PAIRWISE_ONLY'
    if r['JOINT_MATCHES_LEGACY'] == 'yes':
        return 'AGREE_JOINT_ONLY'
    return 'DISAGREE_BOTH'
leg['LEGACY_AGREEMENT_CATEGORY'] = leg.apply(leg_cat, axis=1)
legacy_counts = leg['LEGACY_AGREEMENT_CATEGORY'].value_counts().to_dict()
copies = copies.merge(leg[['SPECIES', 'PROTEIN_ID', 'LEGACY_AGREEMENT_CATEGORY']],
                      on=['SPECIES', 'PROTEIN_ID'], how='left')

# ---- Percentile ranks among all copies: the percentage of the other copies with a strictly lower value ----
for col, out in [('CBS_NONCONTACT_IDENTITY', f'IDENTITY_PERCENTILE_IN_{N_COPIES}'),
                 ('MIN_PLDDT_MEDIAN_PERCOPY', f'PLDDT_PERCENTILE_IN_{N_COPIES}')]:
    v = pd.to_numeric(copies[col], errors='coerce')
    copies[out] = (v.rank(method='min') - 1) / (v.notna().sum() - 1) * 100

METRIC_COLS = ['QUERY', 'SPECIES', 'PROTEIN_ID', 'COCCIDIAN', 'APICOMPLEXAN', 'FULL_LENGTH',
               'CBS_NONCONTACT_N_CALLABLE', 'CBS_NONCONTACT_N_IDENTICAL', 'CBS_NONCONTACT_IDENTITY',
               f'IDENTITY_PERCENTILE_IN_{N_COPIES}',
               'AMP_CONTACT_SCORE_SITE1_PAIRWISE', 'N_CALLABLE_SITE1_PAIRWISE',
               'SITE1_N_POSITIONS', 'SITE1_PW_JOINT_SAME_OBSERVED_AA', 'SITE1_EXACT_PAIRWISE', 'SITE1_EXACT_JOINT',
               'PAIRWISE_JOINT_IDENTICAL_SITE1',
               'N_CBS_EVIDENCE_HITS', 'OBSERVED_CBS_UNION_LENGTH', 'OBSERVED_CBS_UNION_FRACTION',
               'STRUCT_N_SITES', 'MIN_PLDDT_MEDIAN_PERCOPY', f'PLDDT_PERCENTILE_IN_{N_COPIES}',
               'RMSD_MEDIAN_PERCOPY', 'ORIGINAL_GATE_PASS_MIN', 'ORIGINAL_GATE_PASS_SITES',
               'LEGACY_AGREEMENT_CATEGORY']
met = copies[METRIC_COLS].sort_values(['COCCIDIAN', 'SPECIES'], ascending=[False, True])
met.to_csv(OUT / 'METRICS_PER_COPY.tsv', sep='\t', index=False)
assert len(met) == N_COPIES

# ---- LOSS_PATTERN_APICOMPLEXA.tsv: the AMPKγ state and copy number of each of the 41 apicomplexan species ----
g = t2[t2['TARGET'] == 'AMPKg'].copy()
assert len(g) == 210
copy_counts = copies.groupby('SPECIES').size()
loss_rows = []
for sp in sorted(api_species):
    row = g[g['SPECIES'] == sp]
    assert len(row) == 1, sp
    r = row.iloc[0]
    lin = tax.loc[tax['SPECIES_KEY'] == sp, 'OFFICIAL_LINEAGE_WITH_RANKS'].iloc[0]
    genus = sp.split('_')[0]
    note = ''
    if genus == 'Eimeria' and r['STATE'] == 'ABSENT':
        n_eimeria_absent = int((g.loc[g['SPECIES'].str.startswith('Eimeria_'), 'STATE'] == 'ABSENT').sum())
        note = 'Eimeria absence block (8 species)' if n_eimeria_absent == 8 else 'Eimeria absence block (%d of 8 species)' % n_eimeria_absent
    elif genus == 'Eimeria':
        note = 'Eimeria genus, AMPKg retained after the owner-approved EG02 rescue; coccidian cohort' if sp in COCCIDIAN else 'Eimeria genus, AMPKg retained'
    elif sp == 'Cyclospora_cayetanensis':
        note = 'Cyclospora absence (sister to Eimeria block)'
    elif sp in COCCIDIAN:
        note = 'coccidian cohort, single retained copy'
    loss_rows.append({
        'SPECIES': sp, 'GENUS': genus,
        'AMPKG_STATE': r['STATE'], 'ABSENCE_CALLABLE': r['ABSENCE_CALLABLE'],
        'DATA_GRADE': r['DATA_GRADE'],
        'N_RETAINED_GAMMA_COPIES': int(copy_counts.get(sp, 0)),
        'IN_COCCIDIAN_COHORT': sp in COCCIDIAN,
        'LINEAGE_HAS_COCCIDIA': ('Coccidia' in lin),
        'NOTE': note,
    })
loss = pd.DataFrame(loss_rows)
loss.to_csv(OUT / 'LOSS_PATTERN_APICOMPLEXA.tsv', sep='\t', index=False)
assert len(loss) == 41

# ---- CONCORDANCE_PER_COPY.tsv: the coccidian copies (9, or 11 with the two rescued Eimeria copies) ----
VAL_NAME = {'Sarco': 'Sarcocystis_neurona', 'Besno': 'Besnoitia_besnoiti', 'Cysto': 'Cystoisospora_suis',
            'C_mur': 'Cryptosporidium_muris', 'C_hom': 'Cryptosporidium_hominis', 'C_par': 'Cryptosporidium_parvum'}
val_rate = val.assign(sp=val['copy'].map(VAL_NAME)).groupby('sp').agg(
    PHASE1_N=('verdict', 'size'),
    PHASE1_REPRODUCED=('verdict', lambda s: int((s == 'MATCH').sum())),
)
val_rate['PHASE1_REPRODUCTION_RATE'] = val_rate['PHASE1_REPRODUCED'] / val_rate['PHASE1_N']
conc = copies[copies['COCCIDIAN']].copy()
conc = conc.merge(val_rate, left_on='SPECIES', right_index=True, how='left')
CONC_COLS = ['QUERY', 'SPECIES', 'PROTEIN_ID',
             'SITE1_PW_JOINT_SAME_OBSERVED_AA', 'PAIRWISE_JOINT_IDENTICAL_SITE1',
             'PHASE1_REPRODUCED', 'PHASE1_N', 'PHASE1_REPRODUCTION_RATE',
             'LEGACY_AGREEMENT_CATEGORY']
# From 23 September 2026: whether the residue register reproduces for the two rescued Eimeria copies, measured in
# the 20-sequence extension of the fixed 18-sequence alignment, with a table layout declared in advance.
if 'register_extension' in SRC:
    EXT_NAME = dict(VAL_NAME)
    EXT_NAME.update({'AG0237': 'Eimeria_necatrix', 'AG0238': 'Eimeria_tenella'})
    ext = pd.read_csv(SRC['register_extension'], sep='\t')
    ext_rate = ext.assign(sp=ext['copy'].map(EXT_NAME)).groupby('sp').agg(
        EXT20_N=('verdict', 'size'),
        EXT20_REPRODUCED=('verdict', lambda s: int((s == 'MATCH').sum())),
    )
    assert ext['copy'].map(EXT_NAME).notna().all(), 'an extension row does not map to a cohort species'
    ext_rate['EXT20_REPRODUCTION_RATE'] = ext_rate['EXT20_REPRODUCED'] / ext_rate['EXT20_N']
    conc['PHASE1_SOURCE_ALIGNMENT'] = np.where(conc['PHASE1_N'].notna(), '18seq_frozen', '')
    conc = conc.merge(ext_rate, left_on='SPECIES', right_index=True, how='left')
    for sp in ('Eimeria_necatrix', 'Eimeria_tenella'):
        m = (conc['SPECIES'] == sp) & conc['PHASE1_N'].isna() & conc['EXT20_N'].notna()
        conc.loc[m, 'PHASE1_N'] = conc.loc[m, 'EXT20_N']
        conc.loc[m, 'PHASE1_REPRODUCED'] = conc.loc[m, 'EXT20_REPRODUCED']
        conc.loc[m, 'PHASE1_REPRODUCTION_RATE'] = conc.loc[m, 'EXT20_REPRODUCTION_RATE']
        conc.loc[m, 'PHASE1_SOURCE_ALIGNMENT'] = '20seq_extension'
    CONC_COLS = CONC_COLS[:-1] + ['PHASE1_SOURCE_ALIGNMENT', 'EXT20_REPRODUCED', 'EXT20_N',
                                  'EXT20_REPRODUCTION_RATE'] + CONC_COLS[-1:]
conc_out = conc[CONC_COLS].sort_values('SPECIES')
conc_out.to_csv(OUT / 'CONCORDANCE_PER_COPY.tsv', sep='\t', index=False)

# ---- Group comparisons: median and quartiles of each measurement ----
coc = copies[copies['COCCIDIAN']]
other_api = copies[copies['APICOMPLEXAN'] & ~copies['COCCIDIAN']]
noncoc = copies[~copies['COCCIDIAN']]
nonapi = copies[~copies['APICOMPLEXAN']]
groups = {f'coccidian_n{len(coc)}': coc, 'other_apicomplexan': other_api, f'non_coccidian_n{len(noncoc)}': noncoc,
          'non_apicomplexan': nonapi, f'all_{len(copies)}': copies}
metrics = {'CBS_NONCONTACT_IDENTITY': None, 'AMP_CONTACT_SCORE_SITE1_PAIRWISE': None,
           'SITE1_PW_JOINT_SAME_OBSERVED_AA': None, 'N_CBS_EVIDENCE_HITS': None,
           'OBSERVED_CBS_UNION_FRACTION': None, 'MIN_PLDDT_MEDIAN_PERCOPY': None,
           'RMSD_MEDIAN_PERCOPY': None, 'FULL_LENGTH': None}
rows = []
for m in metrics:
    r = {'METRIC': m}
    for gname, gdf in groups.items():
        v = pd.to_numeric(gdf[m], errors='coerce').dropna()
        r[f'{gname}_n'] = len(v)
        r[f'{gname}_median'] = round(v.median(), 4) if len(v) else None
        r[f'{gname}_q1'] = round(v.quantile(0.25), 4) if len(v) else None
        r[f'{gname}_q3'] = round(v.quantile(0.75), 4) if len(v) else None
    rows.append(r)
comp = pd.DataFrame(rows)
comp.to_csv(OUT / 'GROUP_COMPARISONS.tsv', sep='\t', index=False)

# Exploratory Mann-Whitney U tests, coccidian against non-coccidian copies, two-sided, not adjusted for multiple tests
mw = {}
for m in ['CBS_NONCONTACT_IDENTITY', 'OBSERVED_CBS_UNION_FRACTION', 'MIN_PLDDT_MEDIAN_PERCOPY']:
    x = pd.to_numeric(coc[m], errors='coerce').dropna()
    y = pd.to_numeric(noncoc[m], errors='coerce').dropna()
    u = stats.mannwhitneyu(x, y, alternative='two-sided')
    mw[m] = {'n_coc': len(x), 'n_noncoc': len(y), 'U': float(u.statistic), 'p_nominal': float(u.pvalue)}

# Counts of AMPKγ states
nonapi_g = g[~g['SPECIES'].isin(api_species)]
summary = {
    'TS': TS,
    'legacy_agreement_counts': legacy_counts,
    'legacy_rows': len(leg),
    't2_gamma_states_all210': g['STATE'].value_counts().to_dict(),
    't2_gamma_states_apicomplexa41': g[g['SPECIES'].isin(api_species)]['STATE'].value_counts().to_dict(),
    't2_gamma_states_nonapicomplexa169': nonapi_g['STATE'].value_counts().to_dict(),
    'apicomplexan_species_with_gamma': int((loss['N_RETAINED_GAMMA_COPIES'] > 0).sum()),
    'other_apicomplexan_copy_rows': int(len(other_api)),
    'other_apicomplexan_species': sorted(other_api['SPECIES'].unique().tolist()),
    'eimeria_absent': loss[(loss['GENUS'] == 'Eimeria')]['AMPKG_STATE'].value_counts().to_dict(),
    'cyclospora_state': loss.loc[loss['SPECIES'] == 'Cyclospora_cayetanensis', 'AMPKG_STATE'].iloc[0],
    'coccidian_single_copy_all': bool((loss.loc[loss['IN_COCCIDIAN_COHORT'], 'N_RETAINED_GAMMA_COPIES'] == 1).all()),
    'mannwhitney_exploratory': mw,
    'coccidian_metrics': coc[['QUERY', 'SPECIES', 'CBS_NONCONTACT_IDENTITY', f'IDENTITY_PERCENTILE_IN_{N_COPIES}',
                              'OBSERVED_CBS_UNION_FRACTION', 'MIN_PLDDT_MEDIAN_PERCOPY', f'PLDDT_PERCENTILE_IN_{N_COPIES}',
                              'RMSD_MEDIAN_PERCOPY', 'AMP_CONTACT_SCORE_SITE1_PAIRWISE',
                              'SITE1_PW_JOINT_SAME_OBSERVED_AA', 'N_CBS_EVIDENCE_HITS', 'FULL_LENGTH',
                              'ORIGINAL_GATE_PASS_MIN', 'LEGACY_AGREEMENT_CATEGORY']].to_dict('records'),
    'group_comparisons': rows,
    'input_hashes': INPUT_HASHES,
}
(OUT / 'ANALYSIS_SUMMARY.json').write_text(json.dumps(summary, indent=2, default=str))
print('legacy agreement:', legacy_counts)
print('other apicomplexan copies:', len(other_api), sorted(other_api['SPECIES'].unique()))
print('t2 apicomplexa:', summary['t2_gamma_states_apicomplexa41'])
print('eimeria:', summary['eimeria_absent'], 'cyclospora:', summary['cyclospora_state'])
print('MW:', json.dumps(mw, indent=1))
print(comp[['METRIC'] + [c for c in comp.columns if c.endswith('_median')]].to_string(index=False))
