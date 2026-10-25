#!/usr/bin/env python
"""J Pre-MD pairs source: parse the published pre-dynamics per-pose tables
(FOLDS.md 992/1211/1731/2008, RESULTS.md 285/684) into (mol, peptide, posekey) -> 'n/ceil'.
Molecule inferred from the peptide label (ox*/bgox31* -> Oxybenzone). Writes premd_pairs.csv beside itself (repo: code/workbook/).
"""
import re, csv, os
REPO = '/Users/cafierom/python_mac/peptidebuilder'
os.chdir(REPO)

def posekey(p):
    p = p.strip()
    pu = p.upper()
    if 'GNINA' in pu:
        m = re.search(r'[pP]\s*(\d+)', pu)
        return 'GNINA' + (m.group(1) if m else '?')   # lane attributed by context in the table loop
    if 'CNNSCORE' in pu:
        return 'cnnscore'
    if 'Boltz co-fold' in p or p == 'co-folded':
        return 'boltz-cofold'
    if 'Vina p1' in p:
        if 'RF3' in p: return 'rf3cof-vina'
        if 'OF3' in p: return 'o3cof-vina'
        return 'boltz-vina'
    if 'RF3 re-fold' in p: return 'rf3cof'
    if 'OF3 re-fold' in p: return 'o3cof'
    return None

rows = {}
conflicts = []
for doc in ('FOLDS.md', 'RESULTS.md'):
    lines = open(doc).read().splitlines()
    for i, l in enumerate(lines):
        if not re.match(r'^\|\s*peptide\s*\|\s*pose\s*\|', l):
            continue
        hdr = [c.strip() for c in l.strip('|').split('|')]
        pj = next(k for k, h in enumerate(hdr) if h.lower() == 'pose')
        pi_hit = [k for k, h in enumerate(hdr) if 'pairs' in h.lower()]
        if not pi_hit:
            continue  # window/residence/Rg tables carry no pairs cell
        pi = pi_hit[0]
        j = i + 2
        cur = None
        seen = []
        while j < len(lines) and lines[j].startswith('|'):
            cells = [c.strip() for c in lines[j].strip().strip('|').split('|')]
            if len(cells) >= 4:
                lab = cells[0].strip().strip('`')
                if lab:
                    cur = re.sub(r'\s*\(.*$', '', lab).strip().strip('`').strip()
                    seen = []
                pk = posekey(cells[pj])
                if pk and pk.startswith('GNINA'):
                    # a GNINA row belongs to the docked lane the table block was discussing:
                    # the re-fold tables (FOLDS.md) group OF3/RF3 rows under one peptide, the
                    # RESULTS.md per-pose table lists the Boltz lane's docks. Key them apart
                    # so the OF/RF picks never collide with the Boltz pick's pose number.
                    ctx = 'rf3' if any('RF3' in s for s in seen) else \
                          'o3' if any('OF3' in s for s in seen) else ''
                    pk = ctx + pk.lower()
                seen.append(cells[pj])
                if cur and cur.startswith(('ox', 'bgox31')):
                    mol = 'Oxybenzone'
                elif cur and cur.startswith(('bg33',)):
                    mol = 'Octinoxate'
                else:
                    mol = 'Oxybenzone' if ('ox' in (cur or '')) else 'Octinoxate'
                key = (mol, cur, pk)
                v = re.sub(r'\*', '', cells[pi]).strip()
                if key[2] and cur and re.fullmatch(r'\d+/\d+', v):
                    if key in rows and rows[key] != v:
                        conflicts.append((key, rows[key], v))
                    rows.setdefault(key, v)
            j += 1

# The committed pair_contacts tables are the canonical per-pose pair source -- the
# markdown tables in FOLDS.md/RESULTS.md were parsed from these, and several poses
# (guided cofolds, re-fold Vina-p1, re-fold GNINA, rf3 refolds) never got a markdown
# row at all. Parse every lane CSV, skip replicate (_r2) legs, and map each leg to the
# key the fill looks up: base has _o3cof/_rf3cof stripped (guided _o3pcof kept), posekey
# derived from the leg suffix. Generator arms (bg33_*/bgox31_*) stay unkeyed -- no
# slots, 2026-10-22 ruling -- so they produce no rows and their cells stay blank.
import os
for mol in ('octinoxate', 'oxybenzone'):
    root = f'runs/{mol}/md'
    for fn in sorted(os.listdir(root)):
        if not fn.startswith('pair_contacts'):
            continue
        for r in csv.DictReader(open(os.path.join(root, fn))):
            leg = r['leg'].split('/md/')[-1].split('/')[0]
            if leg.endswith('_r2') or re.match(r'bg(ox)?\d+_', leg):
                continue  # replicates are not run-1 pose rows; generator arms have no slots
            if leg.endswith('_o3pcof'):
                base, pk = leg, 'o3cof'
            elif leg.endswith('_o3cof'):
                base, pk = leg[:-len('_o3cof')], 'o3cof'
            elif leg.endswith('_rf3cof'):
                base, pk = leg[:-len('_rf3cof')], 'rf3cof'
            elif '_dock' not in leg:
                base, pk = leg, 'boltz-cofold'
            else:
                base = leg.rsplit('_dock', 1)[0]
                pose = leg.rsplit('_dock', 1)[1]
                dock1 = pose == '1'
                if base.endswith('_o3pcof'):
                    pk = 'o3cof-vina' if dock1 else 'gnina' + pose
                elif base.endswith('_o3cof'):
                    base = base[:-len('_o3cof')]
                    pk = 'o3cof-vina' if dock1 else 'o3gnina' + pose
                elif base.endswith('_rf3cof'):
                    base = base[:-len('_rf3cof')]
                    pk = 'rf3cof-vina' if dock1 else 'rf3gnina' + pose
                else:
                    pk = 'boltz-vina' if dock1 else 'gnina' + pose
            val = f"{int(r['adjacent'].split('/')[0]) + int(r['nonadjacent'].split('/')[0])}/{r['ceiling']}"
            key = (mol.capitalize(), base, pk)
            if key in rows and rows[key] != val:
                conflicts.append(('csv', key, rows[key], val))
            rows.setdefault(key, val)
        if key in rows and rows[key] != val:
            conflicts.append((key, rows[key], val))
        rows.setdefault(key, val)

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'premd_pairs.csv'), 'w', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['mol', 'peptide', 'posekey', 'pairs'])
    for k in sorted(rows, key=lambda t: (t[0], t[1], t[2])):
        w.writerow([k[0], k[1], k[2], rows[k]])
print('keys:', len(rows), '| conflicts:', conflicts if conflicts else 'none')