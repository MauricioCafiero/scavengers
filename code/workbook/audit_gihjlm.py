#!/usr/bin/env python
"""Independent audit of workbook columns G,H,I,J,L,M vs raw sources."""
from pathlib import Path
REPO = Path(__file__).resolve().parents[2]
PREMD_CSV = str(REPO / 'code/workbook/premd_pairs.csv')
TOX_CSV = str(REPO / 'code/workbook/peptox_results.csv')
POSTMD_CSV = str(REPO / 'code/workbook/postmd_enclosure.csv')
import openpyxl, csv, json, os, re
from collections import defaultdict

os.chdir('/Users/cafierom/python_mac/peptidebuilder')
wb = openpyxl.load_workbook('PB_runs_accounting.xlsx')

STD = dict(S1_orig='orig_control', S2_ESM2='s2_esm2_control', S3_ESM2='s3_esm2_control',
           S3_orig='s3_orig_control', S3_shuffle='shuffle_control', S3_shuffle_ESM2='shuffle_control_esm0')
FR  = dict(S1_orig='orig_f12', S3_ESM2='s3_esm2_f4', S3_orig='s3_orig_f12')
STDO = dict(OX1_orig='ox1_orig_control', OX1_ESM2='ox1_esm1_control', OX2_orig='ox2_orig_control',
            OX2_ESM2='ox2_esm1_control', OX3_orig='ox3_orig_control', OX2_shuffle='ox2_shuffle',
            OX2_shuffle_ESM2='ox2_shuffle_esm0')
FRO = dict(OX1_orig='ox1_orig_f12', OX1_ESM2='ox1_esm1_f8', OX2_orig='ox2_orig_f8',
           OX2_ESM2='ox2_esm1_f8', OX3_orig='ox3_orig_f8')
BGM = dict(BG33_1='bg33_1', BG33_2='bg33_2', BG33_3='bg33_3', BG33_4='bg33_4',
           BGOX31_2='bgox31_2', BGOX31_3='bgox31_3', BGOX31_4='bgox31_4', BGOX31_5='bgox31_5')
BGSET = set(BGM.values())

def bases(p):
    if p in BGM:
        b = BGM[p]; return dict(Boltzgen=b, OpenFold=b+'_o3cof', RosettaFold=b+'_rf3cof', Guided=None)
    f, s = (FRO.get(p), STDO.get(p)) if (p in FRO or p in STDO) else (FR.get(p), STD.get(p))
    d = {'Boltz': s}
    if f: d.update({'Force12': f, 'Force4': f, 'Force8': f})
    fo = f if f else s
    d['Guided'] = fo + '_o3pcof'; d['OpenFold'] = fo + '_o3cof'; d['RosettaFold'] = fo + '_rf3cof'
    return d

# pepgraph seq-label alias for the octinoxate shell-3 control (dir named by sequence)
SEQNAME = {}
for r in csv.DictReader(open('runs/octinoxate/folds/manifest.csv')):
    SEQNAME[r['peptide']] = r['sequence']
SEQNAME = {STD[p]: SEQNAME[p] + '_control' for p in FR if p in SEQNAME}   # orig_control, s3_*
SEQNAME.update({STD.get(p, ''): '' for p in ()})

pepgraph = json.load(open('/Users/cafierom/python_mac/pepgraph/outputs/inverse/lane_comparison.json'))
PEP = {}
for r in pepgraph['rows']:
    PEP[(r['lane'], r['design'])] = round(r['perplexity'], 1)

# tox
TOX = {}
for r in csv.DictReader(open(TOX_CSV)):
    label = r['peptide'].split(' (')[0]
    plo, phi, pmed = float(r['p_lo10']), float(r['p_hi90']), float(r['p_median'])
    call = 'elevated' if plo > 0.5 else ('low' if phi < 0.5 else 'unresolved')
    TOX[label] = f"{call} ({r['drift_verdict']}, p_med {pmed:.2f})"

# pre-MD enclosure static tables
GEO = {}
for path, mol in (('runs/octinoxate/dock/dock_poses.csv', 'Octinoxate'),
                  ('runs/oxybenzone/dock/dock_poses.csv', 'Oxybenzone'),
                  ('runs/octinoxate/folds/metrics/dock_folds_poses.csv', 'Octinoxate'),
                  ('runs/oxybenzone/folds/metrics/dock_folds_poses.csv', 'Oxybenzone')):
    if not os.path.exists(path): continue
    for r in csv.DictReader(open(path)):
        GEO[(mol, r['structure'], int(r['pose']))] = float(r['enclosed_fraction'])

# static pairs
PREMD = {}
for r in csv.DictReader(open(PREMD_CSV)):
    PREMD[(r['mol'], r['peptide'], r['posekey'])] = r['pairs']

# post-MD pairs
PAIRS = {}
for mol, root in (('Octinoxate', 'runs/octinoxate/md'), ('Oxybenzone', 'runs/oxybenzone/md')):
    for path in os.listdir(root):
        if not path.startswith('pair_contacts'): continue
        for r in csv.DictReader(open(os.path.join(root, path))):
            legfield = r['leg']
            legname = legfield.split('/md/')[-1].split('/')[0] if '/md/' in legfield else legfield
            adj, na = r.get('adjacent', ''), r.get('nonadjacent', '')
            try: ceil = int(adj.split('/')[1]) + int(na.split('/')[1])
            except Exception: ceil = None
            PAIRS[(mol, legname)] = f"{r['pairs_realised']}/{ceil}" if ceil else r['pairs_realised']

# residence
RES = {}
for mol, root in (('Octinoxate', 'runs/octinoxate/md'), ('Oxybenzone', 'runs/oxybenzone/md')):
    for fn, key in (('pose_reproducibility.csv', 'structure'), ('pose_reproducibility_docked.csv', 'leg')):
        p = os.path.join(root, fn)
        if not os.path.exists(p): continue
        for r in csv.DictReader(open(p)):
            v = r.get('residence_pct')
            if v: RES[(mol, r[key])] = float(v)
for log in os.listdir('runs/folds_md'):
    if log.startswith('record_') and log.endswith('.log'):
        molguess = None
        m = re.search(r'residence ([0-9.]+)%', open(os.path.join('runs/folds_md', log)).read())
        leg = log[7:-4]
        if not m: continue
        molguess = 'Oxybenzone' if log.startswith(('record_ox', 'record_bgox')) else 'Octinoxate'
        RES.setdefault((molguess, leg), float(m.group(1)))

GP = {}
for line in open('runs/folds_md/dock_legs.txt'):
    f = line.split()
    if len(f) >= 4 and f[2].isdigit() and f[3] == 'gnina':
        GP.setdefault(f[1], int(f[2]))
DL = defaultdict(set)
for mol, root in (('Octinoxate', 'runs/octinoxate/md'), ('Oxybenzone', 'runs/oxybenzone/md')):
    for d in os.listdir(root):
        m = re.match(r'(.+)_dock(\d+)$', d)
        if m and int(m.group(2)) > 1: DL[(mol, m.group(1))].add(int(m.group(2)))

def fnum(x):
    try: return float(x)
    except Exception: return None

PEP_LB = {'S1_orig': 'orig', 'S2_ESM2': 's2_esm2_control', 'S3_ESM2': 's3_esm2',
          'S3_orig': 's3_orig', 'S3_shuffle': 'shuffle_control', 'S3_shuffle_ESM2': 'shuffle_control_esm0',
          'OX1_orig': 'ox1', 'OX1_ESM2': 'ox1_esm1', 'OX2_orig': 'ox2_orig', 'OX2_ESM2': 'ox2_esm1',
          'OX3_orig': 'ox3_orig', 'OX2_shuffle': 'ox2_shuffle', 'OX2_shuffle_ESM2': 'ox2_shuffle_esm0',
          'BGOX31_2': 'bgox31_2', 'BGOX31_3': 'bgox31_3', 'BGOX31_4': 'bgox31_4', 'BGOX31_5': 'bgox31_5'}

prob = []
for sheet, mol in (('Octinoxate', 'Octinoxate'), ('Oxybenzone', 'Oxybenzone')):
    ws = wb[sheet]
    molkey = mol.lower()
    blocks, sub = [], None
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
        for s2 in blk['subs']:
            ft = s2['ft']; cofrow = s2['row']
            cl = None
            for rr in range(cofrow, 1, -1):
                v = ws.cell(row=rr, column=3).value
                if v: cl = str(v).strip().replace('Boltz`', 'Boltz').strip(); break
            if ft.startswith('Force') and ft in bs: base = bs[ft]
            elif ft.startswith('Guided'): base = bs.get('Guided')
            elif cl == 'Boltzgen': base = bs.get('Boltzgen')
            elif ft == 'Std': base = bs.get(cl)
            else: base = None
            if base is None: prob.append((sheet, cofrow, 'base unresolved')); continue
            for (r, mk) in s2['cells']:
                leg = {'Co-fold': base, 'Vina 1': base + '_dock1'}.get(mk)
                if mk not in ('Co-fold', 'Vina 1'):
                    poseN = GP.get(base) or (min(DL[(mol, base)]) if DL.get((mol, base)) else None)
                    leg = base + ('_dock%d' % poseN) if poseN else None
                # ---- G ----
                if mk == 'Co-fold':
                    lane = {'Boltz': 'boltz', 'Boltzgen': 'boltzcf', 'OpenFold': 'openfold3',
                            'RosettaFold': 'rf3'}[cl]
                    if ft.startswith('Guided'): lane = 'openfold3p'
                    peek = BGM.get(p) or (base.rsplit('_o3pcof', 1)[0] if base.endswith('_o3pcof')
                        else base.rsplit('_o3cof', 1)[0] if base.endswith('_o3cof')
                        else base.rsplit('_rf3cof', 1)[0] if base.endswith('_rf3cof') else base)
                    gv = PEP.get((lane, peek))
                    if gv is None and cl in ('Boltz', 'Boltzgen'):
                        gv = PEP.get((lane, SEQNAME.get(peek, '\0')))
                    c = ws.cell(row=r, column=7).value
                    if gv is None:
                        if c is not None: prob.append((sheet, r, 'G', f'filled {c} but no PEP({lane},{peek})'))
                    elif fnum(c) is None or abs(fnum(c) - gv) > 0.05:
                        prob.append((sheet, r, 'G', f'expected {gv} ({lane},{peek})', str(c)))
                if mk == 'Co-fold':
                    t = TOX.get(p)
                    c = ws.cell(row=r, column=8).value
                    if t and c is not None and not str(c).startswith(t.split(' (')[0]):
                        prob.append((sheet, r, 'H', f'expected {t}', str(c)))
                # ---- I ----
                poseI = 0 if mk == 'Co-fold' else (1 if mk == 'Vina 1' else (GP.get(base) or (min(DL[(mol, base)]) if DL.get((mol, base)) else None)))
                giv = GEO.get((mol, base, poseI)) if poseI is not None else None
                c = ws.cell(row=r, column=9).value
                if giv is None:
                    if c is not None: prob.append((sheet, r, 'I', f'filled {c} but no GEO row'))
                elif fnum(c) is None or abs(fnum(c) - giv) > 0.005:
                    prob.append((sheet, r, 'I', f'expected {giv}', str(c)))
                # ---- J ----
                jb, jpk = None, None
                if cl in ('Boltz', 'Boltzgen'):
                    jb = base
                    jpk = {'Co-fold': 'boltz-cofold', 'Vina 1': 'boltz-vina'}.get(mk, f'gnina{poseI}')
                elif cl == 'OpenFold':
                    jb = base.rsplit('_o3cof', 1)[0]
                    jpk = {'Co-fold': 'o3cof', 'Vina 1': 'o3cof-vina'}.get(
                        mk, f"gnina{poseI}" if jb.endswith('_o3pcof') else f"o3gnina{poseI}")
                elif cl == 'RosettaFold':
                    jb = base.rsplit('_rf3cof', 1)[0]
                    jpk = 'rf3cof' if mk == 'Co-fold' else ('rf3cof-vina' if mk == 'Vina 1' else f'rf3gnina{poseI}')
                jv = PREMD.get((mol, jb, jpk)) or PREMD.get((sheet, jb, jpk))
                c = ws.cell(row=r, column=10).value
                if jv is None:
                    if c is not None: prob.append((sheet, r, 'J', f'filled {c} but no PREMD({jb},{jpk})'))
                elif str(c) != str(jv):
                    prob.append((sheet, r, 'J', f'expected {jv} ({jb},{jpk})', str(c)))
                # ---- K spot-existence: recomputed geom values live in /tmp/postmd_enclosure.csv ----
                kv = None
                if leg:
                    for r2 in csv.DictReader(open(POSTMD_CSV)):
                        if r2['mol'] == molkey and r2['leg'] == leg: kv = float(r2['enclosed']); break
                c = ws.cell(row=r, column=11).value
                if kv is None:
                    if c is not None: prob.append((sheet, r, 'K', f'filled {c} but no enclosure row'))
                elif fnum(c) is None or abs(fnum(c) - kv) > 0.005:
                    prob.append((sheet, r, 'K', f'expected {kv}', str(c)))
                # ---- L ----
                lv = PAIRS.get((mol, leg)) if leg else None
                c = ws.cell(row=r, column=12).value
                if lv is None:
                    if c is not None: prob.append((sheet, r, 'L', f'filled {c} but no pair row ({leg})'))
                elif str(c) != str(lv):
                    prob.append((sheet, r, 'L', f'expected {lv}', str(c)))
                # ---- M ----
                mv = RES.get((mol, leg)) if leg else None
                c = ws.cell(row=r, column=13).value
                if mv is None:
                    if c is not None: prob.append((sheet, r, 'M', f'filled {c} but no residence ({leg})'))
                elif abs(fnum(c) - mv) > 0.05:
                    prob.append((sheet, r, 'M', f'expected {mv}', str(c)))

from collections import Counter
print(Counter(p[2] for p in prob))
print(len(prob), 'problems')
for x in prob: print(x)