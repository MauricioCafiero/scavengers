#!/usr/bin/env python
"""Fill F-V of PB_runs_accounting.xlsx from the repo's published tables.

Sources, per column:
  F strain       morph results/strain_summary.csv (dE_strain_kcal), cofold rows only
  G perplexity   pepgraph outputs/inverse/lane_comparison.json, cofold rows only
  H tox          /tmp/peptox_results.csv: peptox ensemble call + drift verdict, Co-fold rows only
  I pre-MD enc   dock_folds_poses.csv / dock_poses.csv enclosed_fraction (pose 0 / 1 / gnina pick)
  J pre-MD pairs /tmp/premd_pairs.csv (from j_pairs.py): the published pre-dynamics per-pose
                 tables in FOLDS.md/RESULTS.md, keyed (mol, peptide, posekey)
  K post-MD enc  (no post-MD geometry table - left blank; would need check_fold.geometry()
                  on frame_medoid)
  L post-MD pair pair_contacts_*.csv pairs_realised / ceiling (design slots only)
  M retention    pose_reproducibility(_docked).csv residence_pct -> record_<leg>.log ->
                 cofold_binding_tables.md
  N..Q DG1+ws    mmgbsa_summary.csv (prod_20ns, first_5/10/15ns)
  R..V DG2..ave  mmgbsa_replicates.csv for <leg>_r2, V = mean of DG1/DG2
  W/X/Y          K/L/M for the replicate <leg>_r2 (postmd_enclosure.csv _r2 rows, written by
                 postmd_enclosure.py --add-replicates; pair_contacts _r2 rows; record_<leg>_r2.log)
Only fills cells that are currently empty. Idempotent.
"""
from pathlib import Path
REPO = Path(__file__).resolve().parents[2]
PREMD_CSV = str(REPO / 'code/workbook/premd_pairs.csv')
TOX_CSV = str(REPO / 'code/workbook/peptox_results.csv')
POSTMD_CSV = str(REPO / 'code/workbook/postmd_enclosure.csv')
import openpyxl, os, re, csv, json
from collections import defaultdict

REPO = '/Users/cafierom/python_mac/peptidebuilder'
os.chdir(REPO)
BGDIR = os.path.expanduser('~/python_mac/boltzgen_local/md')

STD = dict(S1_orig='orig_control', S2_ESM2='s2_esm2_control', S3_ESM2='s3_esm2_control',
           S3_orig='s3_orig_control', S3_shuffle='shuffle_control', S3_shuffle_ESM2='shuffle_control_esm0')
FORCED = dict(S1_orig='orig_f12', S3_ESM2='s3_esm2_f4', S3_orig='s3_orig_f12')
STD_O = dict(OX1_orig='ox1_orig_control', OX1_ESM2='ox1_esm1_control', OX2_orig='ox2_orig_control',
             OX2_ESM2='ox2_esm1_control', OX3_orig='ox3_orig_control', OX2_shuffle='ox2_shuffle',
             OX2_shuffle_ESM2='ox2_shuffle_esm0')
FORCED_O = dict(OX1_orig='ox1_orig_f12', OX1_ESM2='ox1_esm1_f8', OX2_orig='ox2_orig_f8',
                OX2_ESM2='ox2_esm1_f8', OX3_orig='ox3_orig_f8')
BG = dict(BG33_1='bg33_1', BG33_2='bg33_2', BG33_3='bg33_3', BG33_4='bg33_4',
          BGOX31_2='bgox31_2', BGOX31_3='bgox31_3', BGOX31_4='bgox31_4', BGOX31_5='bgox31_5')

def bases(p):
    if p in BG:
        b = BG[p]; return dict(Boltzgen=b, OpenFold=b + '_o3cof', RosettaFold=b + '_rf3cof')
    if p in FORCED_O or p in STD_O:
        f, s = FORCED_O.get(p), STD_O.get(p)
    else:
        f, s = FORCED.get(p), STD.get(p)
    d = {'Boltz': s}
    if f: d['Force12'] = d['Force4'] = d['Force8'] = f
    fo = f if f else s
    d['Guided'] = fo + '_o3pcof'; d['OpenFold'] = fo + '_o3cof'; d['RosettaFold'] = fo + '_rf3cof'
    return d

# ---- sources ----
def load_summaries():
    out = {}
    for mol, root in (('Octinoxate', 'runs/octinoxate/md'), ('Oxybenzone', 'runs/oxybenzone/md')):
        d = defaultdict(dict)
        for path in (f'{root}/mmgbsa_summary.csv', f'{root}/mmgbsa_replicates.csv'):
            label = 'dg2' if 'replicates' in path else 'dg1'
            for r in csv.DictReader(open(path)):
                d[(label, r['structure'])][r['leg']] = (float(r['dg_bind']), float(r.get('std_err') or 0))
        out[mol] = d
    return out

def load_residence():
    # (mol, leg) -> retention %
    res = {}
    for mol, root in (('Octinoxate', 'runs/octinoxate/md'), ('Oxybenzone', 'runs/oxybenzone/md')):
        for fn, key in (('pose_reproducibility.csv', 'structure'), ('pose_reproducibility_docked.csv', 'leg')):
            path = f'{root}/{fn}'
            if os.path.exists(path):
                for r in csv.DictReader(open(path)):
                    v = r.get('residence_pct')
                    if v:
                        res[(mol, r[key])] = float(v)
        # record logs
        for log in os.listdir('runs/folds_md'):
            if log.startswith('record_') and log.endswith('.log'):
                m = re.search(r'residence ([0-9.]+)%', open(f'runs/folds_md/{log}').read())
                if m: res.setdefault((mol, log[7:-4]), float(m.group(1)))
    return res

def load_strain():
    out = {}
    for r in csv.DictReader(open(os.path.expanduser('~/python_mac/morph/results/strain_summary.csv'))):
        out[('boltzapo', r['mol'], r['peptide'])] = r
        out[('of3apo', r['mol'], r['peptide'])] = r
        out[('rf3apo', r['mol'], r['peptide'])] = r
    return out

def load_pep():
    d = json.load(open('/Users/cafierom/python_mac/pepgraph/outputs/inverse/lane_comparison.json'))
    out = {}
    for r in d['rows']:
        out[(r['lane'], r['design'])] = r['perplexity']
    return out

def load_pose_geom():
    out = {}
    for path in ('runs/octinoxate/dock/dock_poses.csv', 'runs/oxybenzone/dock/dock_poses.csv',
                 'runs/octinoxate/folds/metrics/dock_folds_poses.csv',
                 'runs/oxybenzone/folds/metrics/dock_folds_poses.csv'):
        mol = 'Octinoxate' if 'octinoxate' in path else 'Oxybenzone'
        if not os.path.exists(path): continue
        for r in csv.DictReader(open(path)):
            out[(mol, r['structure'], int(r['pose']))] = float(r['enclosed_fraction'])
    return out

def load_pairs():
    import glob
    out = {}
    for mol, root in (('Octinoxate', 'runs/octinoxate/md'), ('Oxybenzone', 'runs/oxybenzone/md')):
        for path in glob.glob(f'{root}/pair_contacts*.csv'):
            for r in csv.DictReader(open(path)):
                adj, na = r.get('adjacent', ''), r.get('nonadjacent', '')
                try:
                    ceil = int(adj.split('/')[1]) + int(na.split('/')[1])
                except Exception:
                    ceil = None
                legfield = r['leg']
                legname = legfield.split('/md/')[-1].split('/')[0] if '/md/' in legfield else legfield
                out[(mol, legname)] = (r['pairs_realised'], ceil)
    return out

def load_cofold_tables():
    # residence fallback: cofold_binding_tables.md
    txt = open('runs/folds_md/cofold_binding_tables.md').read()
    out = {}
    # map sheet-style peptide label -> cofold-table label
    CT = {'S1_orig': 'orig', 'S2_ESM2': 's2_esm2_control', 'S3_ESM2': 's3_esm2',
          'S3_orig': 's3_orig', 'S3_shuffle': 'shuffle_control', 'S3_shuffle_ESM2': 'shuffle_control_esm0',
          'BG33_1': 'bg33_1', 'BG33_2': 'bg33_2', 'BG33_3': 'bg33_3', 'BG33_4': 'bg33_4',
          'OX1_orig': 'ox1', 'OX1_ESM2': 'ox1_esm1', 'OX2_orig': 'ox2_orig', 'OX2_ESM2': 'ox2_esm1',
          'OX3_orig': 'ox3_orig', 'OX2_shuffle': 'ox2_shuffle', 'OX2_shuffle_ESM2': 'ox2_shuffle_esm0',
          'BGOX31_2': 'bgox31_2', 'BGOX31_3': 'bgox31_3', 'BGOX31_4': 'bgox31_4', 'BGOX31_5': 'bgox31_5'}
    tabs = re.split(r'\*\*Table \d+[.:][^*]*\*\*', txt)[1:]
    molof = ['Octinoxate', 'Oxybenzone']
    rows = re.findall(r'\*\*(Table \d+)', txt)
    # simpler: split on Table headers with molecule names
    return out  # resolved below on demand

def cofold_table_res():
    txt = open('runs/folds_md/cofold_binding_tables.md').read()
    molspl = re.split(r'\*\*Table \d+\. (\w+):', txt)
    res = {}
    for i in range(1, len(molspl) - 1, 2):
        molname, body = molspl[i], molspl[i + 1]
        mol = molname if molname in ('Octinoxate', 'Oxybenzone') else molname.title()
        for line in body.splitlines():
            m = re.match(r'\| `?([\w]+)`? \| (co-fold|Vina p1|GNINA pick|[^\s/]+[^|]*) \|(.*)', line, re.I)
            if not m: continue
            lab = m.group(1).strip().strip('`')
            poset = m.group(2).strip()
            poset = 'co-fold' if poset.lower().startswith('co-fold') else ('vina' if 'vina' in poset.lower() else 'gnina')
            cells = m.group(3).strip().strip('|').split('|')
            for ci, cell in enumerate(cells):
                col = ['Boltz', 'OF3', 'RF', 'OF3p'][ci]
                rm = re.search(r'\(([0-9.]+)%\)', cell)
                if not rm: continue
                res[(mol, lab, col, poset)] = float(rm.group(1))
    return res

SUM = load_summaries(); RES = load_residence(); STR = load_strain(); PEP = load_pep()
GEO = load_pose_geom(); PAIRS = load_pairs(); CTR = cofold_table_res()

def load_tox():
    # peptox infer: elevated if p_lo10 > 0.5, low if p_hi90 < 0.5, else unresolved
    out = {}
    for r in csv.DictReader(open(TOX_CSV)):
        label = r['peptide'].split(' (')[0]
        plo, phi, pmed = float(r['p_lo10']), float(r['p_hi90']), float(r['p_median'])
        call = 'elevated' if plo > 0.5 else ('low' if phi < 0.5 else 'unresolved')
        out[label] = f"{call} ({r['drift_verdict']}, p_med {pmed:.2f})"
    return out

def load_j():
    out = {}
    for r in csv.DictReader(open(PREMD_CSV)):
        out[(r['mol'], r['peptide'], r['posekey'])] = r['pairs']
    return out

def load_k():
    out = {}
    for r in csv.DictReader(open(POSTMD_CSV)):
        out[(r['mol'], r['leg'])] = float(r['enclosed'])
    return out

TOX = load_tox(); PREMD = load_j(); KPOST = load_k()

def j_pkey(mk, cl, base, pose):
    """Return the premd_pairs lookup: posekey for this row; for cofold-lane rows the
    peptide key is the design base with the lane suffix stripped."""
    if cl in ('Boltz', 'Boltzgen'):
        return base, {'Co-fold': 'boltz-cofold', 'Vina 1': 'boltz-vina'}.get(mk, f'gnina{pose}')
    if cl == 'OpenFold':
        fob = base.rsplit('_o3cof', 1)[0] if base.endswith('_o3cof') else base
        gk = f'gnina{pose}' if base.endswith('_o3pcof') else f'o3gnina{pose}'
        return fob, {'Co-fold': 'o3cof', 'Vina 1': 'o3cof-vina'}.get(mk, gk)
    if cl == 'RosettaFold':
        fob = base.rsplit('_rf3cof', 1)[0] if base.endswith('_rf3cof') else base
        pk = 'rf3cof' if mk == 'Co-fold' else ('rf3cof-vina' if mk == 'Vina 1' else f'rf3gnina{pose}')
        return fob, pk
    return base, None

PEPLAB_TO_CT = {'S1_orig': 'orig', 'S2_ESM2': 's2_esm2_control', 'S3_ESM2': 's3_esm2',
                'S3_orig': 's3_orig', 'S3_shuffle': 'shuffle_control', 'S3_shuffle_ESM2': 'shuffle_control_esm0',
                'BG33_1': 'bg33_1', 'BG33_2': 'bg33_2', 'BG33_3': 'bg33_3', 'BG33_4': 'bg33_4',
                'OX1_orig': 'ox1', 'OX1_ESM2': 'ox1_esm1', 'OX2_orig': 'ox2_orig', 'OX2_ESM2': 'ox2_esm1',
                'OX3_orig': 'ox3_orig', 'OX2_shuffle': 'ox2_shuffle', 'OX2_shuffle_ESM2': 'ox2_shuffle_esm0',
                'BGOX31_2': 'bgox31_2', 'BGOX31_3': 'bgox31_3', 'BGOX31_4': 'bgox31_4', 'BGOX31_5': 'bgox31_5'}

def pose_label(mk, leg_pose):
    if mk == 'Co-fold': return 'co-fold'
    if mk == 'Vina 1': return 'vina p1'
    return None  # gnina handled separately by pose

def gnina_pose_from_md(pmol, base, mdroot):
    return None

wb = openpyxl.load_workbook('PB_runs_accounting.xlsx')
report = defaultdict(list)
for sheet, mol in (('Octinoxate', 'Octinoxate'), ('Oxybenzone', 'Oxybenzone')):
    ws = wb[sheet]
    mdroot = f'runs/{mol.lower()}/md'
    summary = SUM[mol]
    # disk dock legs per base (for gnina pick resolution)
    docklegs = defaultdict(set)
    for pre in [mdroot] + [BGDIR]:
        if not os.path.isdir(pre): continue
        for d in os.listdir(pre):
            if re.search(r'_dock\d+$', d) and not re.search(r'_r\d+$', d):
                docklegs[d.rsplit('_dock', 1)[0]].add(d)
    GP = {}
    for line in open('runs/folds_md/dock_legs.txt'):
        f = line.split()
        if len(f) >= 4 and f[2].isdigit() and f[3] == 'gnina': GP[f[1]] = f[2]
    blocks = []; sub = None
    for r in range(2, ws.max_row + 1):
        b = ws.cell(row=r, column=2).value
        if b and 'Max MD runs' in str(b): break
        if b: blocks.append({'p': str(b).strip(), 'subs': []}); sub = None
        e = ws.cell(row=r, column=5).value
        if e:
            mk = str(e).strip()
            if mk == 'Co-fold':
                sub = {'ft': str(ws.cell(row=r, column=4).value or '').strip(), 'row': r, 'cells': [(r, mk)]}
                blocks[-1]['subs'].append(sub)
            elif sub is not None: sub['cells'].append((r, mk))
    for blk in blocks:
        p = blk['p']; bs = bases(p)
        for sub2 in blk['subs']:
            ft = sub2['ft']; cofrow = sub2['row']
            cl = None
            for rr in range(cofrow, 1, -1):
                v = ws.cell(row=rr, column=3).value
                if v: cl = str(v).strip(); break
            cl = (cl or '?').replace('Boltz`', 'Boltz')
            lane_key = cl  # 'Boltz'/'Boltzgen'/'OpenFold'/'RosettaFold'
            if ft.startswith('Force') and ft in bs: base = bs[ft]
            elif ft.startswith('Guided'): base = bs.get('Guided')
            elif cl == 'Boltzgen': base = bs.get('Boltzgen')
            elif ft == 'Std': base = bs.get(cl)
            else: base = None
            if base is None: report[sheet].append(f"E{cofrow} {p} {cl} {ft}: base unresolved"); continue
            for (r, mk) in sub2['cells']:
                if mk == 'Co-fold': leg, dock_struct, pose = base, base, 0
                elif mk == 'Vina 1': leg, dock_struct, pose = base + '_dock1', base, 1
                else:
                    nonp1 = sorted(x for x in docklegs.get(base, set()) if x != base + '_dock1')
                    if not nonp1: continue
                    pose = int(GP.get(base, nonp1[0][-1]))
                    leg = base + '_dock%d' % pose; dock_struct = base
                def put(col, val, fmt=None, why=None):
                    cell = ws.cell(row=r, column=col)
                    if cell.value is not None and cell.value != '':
                        report[sheet].append(f"E{r} col {col}: already {cell.value!r}, kept")
                        return
                    if val is None:
                        report[sheet].append(f"E{r} col {col}: BLANK ({why})"); return
                    cell.value = val
                    if fmt: cell.number_format = fmt
                # N..V: dG block
                def dg(molsum_label, legname):
                    d = summary.get((molsum_label, legname), {})
                    prod = d.get('prod_20ns') or d.get('prod_L1_modal')
                    w5, w10, w15 = d.get('first_5ns'), d.get('first_10ns'), d.get('first_15ns')
                    if prod is not None:
                        w5 = w5 or (0, 0); w10 = w10 or (0, 0); w15 = w15 or (0, 0)
                        return round(prod[0], 2), round(w5[0], 2), round(w10[0], 2), round(w15[0], 2)
                    return None
                g1 = dg('dg1', leg)
                put(14, g1[0] if g1 else None, '0.00', f'no DG1 for {leg}')
                if g1:
                    for i, col in enumerate((15, 16, 17)): put(col, g1[1 + i], '0.00')
                g2 = dg('dg2', leg + '_r2')
                put(18, g2[0] if g2 else None, '0.00', f'no replicate for {leg}')
                if g2:
                    for i, col in enumerate((19, 20, 21)): put(col, g2[1 + i], '0.00')
                    if g1: put(22, round((g1[0] + g2[0]) / 2, 2), '0.00')
                # M retention
                ret = RES.get((mol, leg))
                if ret is None:
                    col = {'Boltz': ['Boltz', 'Boltzgen'], 'OpenFold': ['OF3'], 'RosettaFold': ['RF'],
                           'Boltzgen': ['Boltz']}[cl]
                    poselab = {'Co-fold': 'co-fold', 'Vina 1': 'vina'}.get(mk, 'gnina')
                    for c in col:
                        v2 = CTR.get((mol, PEPLAB_TO_CT.get(p, ''), c, poselab))
                        if v2 is not None:
                            ret = v2; break
                put(13, round(ret, 1) if ret is not None else None, '0.0', f'no residence for {leg}')
                # F/G/H (cofold rows only)
                if mk == 'Co-fold':
                    # morph/pepgraph key the bare base; cofold lanes append a suffix
                    pbase = base
                    for suf in ('_o3pcof', '_o3cof', '_rf3cof'):
                        if pbase.endswith(suf): pbase = pbase[:-len(suf)]
                    slane = 'boltzapo' if cl in ('Boltz', 'Boltzgen') else \
                            ('rf3apo' if cl == 'RosettaFold' and not ft.startswith('Guided') else 'of3apo')
                    mk_key = mol.lower()
                    row = STR.get((slane, mk_key, pbase))
                    s_val = float(row['dE_strain_kcal']) if row else None
                    put(6, round(s_val, 1) if s_val is not None else None, '0.0',
                        f'no morph strain for {pbase}/{slane}')
                    planelane = {'Boltz': 'boltz', 'Boltzgen': 'boltzcf', 'OpenFold': 'openfold3',
                                 'RosettaFold': 'rf3'}[cl]
                    if ft.startswith('Guided'): planelane = 'openfold3p'
                    peek = BG[p] if p in BG else pbase
                    pv = PEP.get((planelane, peek))
                    put(7, round(pv, 1) if pv is not None else None, '0.0',
                        f'no pepgraph row ({planelane},{base})')
                    put(8, TOX.get(p), None, f'no peptox result for {p}')
                # I pre-MD enclosure
                ge = GEO.get((mol, dock_struct, pose))
                if ge is None and p in BG:
                    ge = GEO.get((mol, base, pose))
                put(9, ge, '0.00', f'no dock-pose row ({dock_struct}, pose {pose})')
                # J pre-MD pairs: published per-pose static tables
                jbase, jpk = j_pkey(mk, lane_key, base, pose)
                jv = PREMD.get((mol, jbase, jpk))
                if jv is None:
                    jv = PREMD.get((sheet, jbase, jpk))
                put(10, jv, None, f'no static pairs cell ({jbase}, {jpk})')
                # K post-MD enclosure: check_fold.geometry on frame_medoid
                kv = KPOST.get((mol.lower(), leg))
                put(11, kv, '0.00', f'no post-MD enclosure leg {leg}')
                # L post-MD pairs
                pr = PAIRS.get((mol, leg))
                put(12, f'{pr[0]}/{pr[1]}' if pr and pr[1] else (pr[0] if pr else None), None,
                    f'no pair_contacts row for {leg}')
                # W/X/Y: K/L/M again for the replicate (<leg>_r2), only where a replicate exists
                if g2:
                    leg2 = leg + '_r2'
                    put(23, KPOST.get((mol.lower(), leg2)), '0.00', f'no post-MD enclosure leg {leg2}')
                    pr2 = PAIRS.get((mol, leg2))
                    put(24, f'{pr2[0]}/{pr2[1]}' if pr2 and pr2[1] else (pr2[0] if pr2 else None), None,
                        f'no pair_contacts row for {leg2}')
                    ret2 = RES.get((mol, leg2))
                    put(25, round(ret2, 1) if ret2 is not None else None, '0.0', f'no residence for {leg2}')

# Headers for W/X/Y: K/L/M's own header text with " 2" (as DG1 -> DG2), in K/L/M's style.
from copy import copy
for ws in wb.worksheets:
    for col, src in ((23, 11), (24, 12), (25, 13)):
        if ws.cell(1, col).value is None:
            ws.cell(1, col).value = f'{ws.cell(1, src).value} 2'
            ws.cell(1, col)._style = copy(ws.cell(1, src)._style)
            L, S = openpyxl.utils.get_column_letter(col), openpyxl.utils.get_column_letter(src)
            ws.column_dimensions[L].width = ws.column_dimensions[S].width

wb.save('PB_runs_accounting.xlsx')
print('saved.')
for s, items in report.items():
    print(f'\n===== {s} =====')
    for x in items: print(' ', x)