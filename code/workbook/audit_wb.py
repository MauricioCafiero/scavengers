#!/usr/bin/env python
"""Independent audit of PB_runs_accounting.xlsx columns F..V against raw sources.

Re-derives every cofold-row cell from the primary tables (NOT reusing the fill
scripts' intermediates /tmp/passcsv except where that file IS the source-of-truth
parse) and diffs against the workbook. Column index map: workbook col c -> dump row[c+1].
"""
import openpyxl, csv, json, os, re
from collections import defaultdict

os.chdir('/Users/cafierom/python_mac/peptidebuilder')
wb = openpyxl.load_workbook('PB_runs_accounting.xlsx')
fnum = lambda x: (float(x) if str(x).replace('.','',1).replace('-','',1).isdigit() else None)

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

def bases(p):
    if p in BGM:
        b = BGM[p]; return dict(Boltzgen=b, OpenFold=b+'_o3cof', RosettaFold=b+'_rf3cof', Guided=None)
    f, s = (FRO.get(p), STDO.get(p)) if (p in FRO or p in STDO) else (FR.get(p), STD.get(p))
    d = {'Boltz': s}
    if f: d.update({'Force12': f, 'Force4': f, 'Force8': f})
    fo = f if f else s
    d['Guided'] = fo + '_o3pcof'; d['OpenFold'] = fo + '_o3cof'; d['RosettaFold'] = fo + '_rf3cof'
    return d

# ---- raw sources ----
def load_dg(path):
    d = defaultdict(dict)
    for r in csv.DictReader(open(path)):
        d[r['structure']][r['leg']] = float(r['dg_bind'])
    return d
SUM = {'Octinoxate': load_dg('runs/octinoxate/md/mmgbsa_summary.csv'),
       'Oxybenzone': load_dg('runs/oxybenzone/md/mmgbsa_summary.csv')}
REP = {m: load_dg(f'runs/{m.lower()}/md/mmgbsa_replicates.csv') for m in SUM}
ALL = defaultdict(dict)
for m in SUM:
    for s, ld in SUM[m].items(): ALL[m][s] = dict(ld)
    for s, ld in REP[m].items(): ALL[m].setdefault(s, {}).update(ld)

GP = {}   # base -> gnina pose from dock_legs.txt
for line in open('runs/folds_md/dock_legs.txt'):
    f = line.split()
    if len(f) >= 4 and f[2].isdigit() and f[3] == 'gnina':
        GP.setdefault(f[1], set()).add(int(f[2]))

def prod(tbl):
    return tbl.get('prod_20ns') or tbl.get('prod_L1_modal')

prob, checked = [], 0
for sheet, molkey in (('Octinoxate', 'octinoxate'), ('Oxybenzone', 'oxybenzone')):
    ws = wb[sheet]
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
        bs = bases(blk['p'])
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
            if base is None:
                prob.append((sheet, cofrow, 'base unresolved')); continue
            for (r, mk) in s2['cells']:
                leg = {'Co-fold': base, 'Vina 1': base + '_dock1'}.get(mk)
                if mk not in ('Co-fold', 'Vina 1'):
                    poses = GP.get(base, set())
                    nonp1 = sorted(int(d.rsplit('_dock', 1)[1]) for d in os.listdir(f"runs/{molkey}/md")
                                   if d.startswith(base + '_dock') and not d.endswith('_r2')
                                   and int(d.rsplit('_dock', 1)[1]) > 1) if False else None
                    if poses: leg = base + '_dock%d' % (max(poses))
                    else:
                        dl = [d for d in os.listdir(f"runs/{molkey}/md")
                              if d.startswith(base + '_dock') and d[len(base):].count('_') in (0, 1)
                              and not d.endswith('_r2') and int(d.rsplit('_dock', 1)[1]) > 1]
                        leg = (base + '_dock%d' % min(int(d.rsplit('_dock', 1)[1]) for d in dl)) if dl else None
                if leg is None:
                    # gnina row with no dock leg: every DG cell should be empty
                    for c in range(14, 23):
                        if fnum(ws.cell(row=r, column=c).value):
                            prob.append((sheet, r, col, 'DG filled but leg unresolved'))
                    continue
                tbl = ALL[sheet].get(leg) or {}
                p1 = prod(tbl)
                if p1 is None:
                    for c in range(14, 23):
                        if fnum(ws.cell(row=r, column=c).value):
                            prob.append((sheet, r, c, f'DG filled but {leg} has no prod row'))
                    continue
                checked += 1
                exp = {14: round(p1, 2)}
                for wk, c in ('first_5ns', 15), ('first_10ns', 16), ('first_15ns', 17):
                    w = tbl.get(wk); exp[c] = round(w, 2) if w is not None else 0.0 if tbl else None
                rp = prod(REP[sheet].get(leg + '_r2') or {})
                if rp is not None:
                    exp[18] = round(rp, 2)
                    rt = REP[sheet].get(leg + '_r2')
                    for wk, c in ('first_5ns', 19), ('first_10ns', 20), ('first_15ns', 21):
                        w = rt.get(wk); exp[c] = round(w, 2) if w is not None else 0.0
                    exp[22] = round((round(p1, 2) + round(rp, 2)) / 2, 2)
                for c in range(14, 23):
                    g = ws.cell(row=r, column=c).value
                    ev = exp.get(c)
                    if ev is None:
                        if fnum(g) is not None:
                            prob.append((sheet, r, c, 'cell filled, expected blank', g, leg))
                    elif abs((fnum(g) if fnum(g) is not None else 1e9) - ev) > 0.005:
                        prob.append((sheet, r, c, f'expected {ev} ({leg})', str(g)))

print(f'{checked} DG blocks checked, {len(prob)} problems')
for p in prob: print(p)