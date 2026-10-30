#!/usr/bin/env python
"""Post-MD enclosure per leg: check_fold.geometry() on the leg's prod_20ns/frame_medoid.pdb.

Writes the rows of code/workbook/postmd_enclosure.csv that fill_accounting.py reads for workbook
column K (run 1) and the run-2 column (legs named <leg>_r2). The run-1 rows were first produced
ad hoc on 2026-10-25 and this script was written afterwards to extend the table to replicates;
`--check` recomputes every existing row and reports any that do not reproduce.

Peptide = every non-LIG heavy atom of the medoid frame, ligand = residue LIG's heavy atoms, in the
frame's own coordinates (both come from the same wrapped solute, so no re-imaging is needed).

    .venv/bin/python code/workbook/postmd_enclosure.py --check            # verify existing rows
    .venv/bin/python code/workbook/postmd_enclosure.py --add-replicates   # append missing _r2 rows

Appends only; never rewrites an existing row.
"""
import csv, os, sys
import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, 'code'))
from check_fold import geometry

CSV = os.path.join(REPO, 'code/workbook/postmd_enclosure.csv')
FIELDS = ['mol', 'leg', 'row', 'enclosed', 'centroid_sep', 'engaged', 'lig_heavy', 'wrapped', 'pdb']


def medoid_path(mol, leg):
    for root in (f'{REPO}/runs/{mol}/md', os.path.expanduser('~/python_mac/boltzgen_local/md')):
        for prod in ('prod_20ns', 'prod_L1_modal'):   # Modal legs keep theirs under prod_L1_modal
            p = f'{root}/{leg}/{prod}/frame_medoid.pdb'
            if os.path.exists(p):
                return p
    return None


def measure(pdb):
    pep, lig = [], []
    for line in open(pdb):
        if not line.startswith(('ATOM', 'HETATM')):
            continue
        el = line[76:78].strip() or line[12:16].strip()[0]
        if el == 'H':
            continue
        xyz = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        (lig if line[17:20] == 'LIG' else pep).append(xyz)
    return geometry(np.array(pep), np.array(lig))


def row_for(mol, leg, wbrow=''):
    p = medoid_path(mol, leg)
    if p is None:
        return None
    g = measure(p)
    return dict(mol=mol, leg=leg, row=wbrow, enclosed=g['enclosed_fraction'],
                centroid_sep=g['centroid_separation'], engaged=g['engaged'],
                lig_heavy=g['ligand_heavy_atoms'], wrapped=g['wrapped_fraction'],
                pdb=p.split('/md/', 1)[1])


def main():
    rows = list(csv.DictReader(open(CSV)))
    if '--check' in sys.argv:
        bad = 0
        for r in rows:
            n = row_for(r['mol'], r['leg'])
            if n is None or abs(float(r['enclosed']) - n['enclosed']) > 0.0005 \
                    or int(r['engaged']) != n['engaged'] or int(r['lig_heavy']) != n['lig_heavy']:
                bad += 1
                print('MISMATCH', r['mol'], r['leg'], r['enclosed'], n and n['enclosed'])
        print(f'{len(rows)} rows checked, {bad} do not reproduce')
    if '--add-replicates' in sys.argv:
        have = {(r['mol'], r['leg']) for r in rows}
        new, missing = [], []
        for r in rows:
            key = (r['mol'], r['leg'] + '_r2')
            if r['leg'].endswith('_r2') or key in have:
                continue
            n = row_for(*key, wbrow=r['row'])
            (new if n else missing).append(n or key)
        with open(CSV, 'a', newline='') as f:
            csv.DictWriter(f, fieldnames=FIELDS).writerows(new)
        print(f'appended {len(new)} _r2 rows; no medoid for {len(missing)}: {missing}')


if __name__ == '__main__':
    main()
