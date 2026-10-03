"""Score nesso against MM/GBSA dG on the six structures that have dynamics.

The README ranks every metric this way (static sum -0.80, UMA -0.80, enclosed
+0.40, VDWAALS +0.80, mean engagement +0.83), so nesso goes in the same column.
"""
import json, os, csv

def spearman(x, y):
    def rank(v):
        s = sorted(range(len(v)), key=lambda i: v[i]); r = [0.0]*len(v); i = 0
        while i < len(s):
            j = i
            while j+1 < len(s) and v[s[j+1]] == v[s[i]]: j += 1
            for k in range(i, j+1): r[s[k]] = (i+j)/2 + 1
            i = j+1
        return r
    rx, ry = rank(x), rank(y); n = len(x)
    mx, my = sum(rx)/n, sum(ry)/n
    num = sum((a-mx)*(b-my) for a, b in zip(rx, ry))
    den = (sum((a-mx)**2 for a in rx) * sum((b-my)**2 for b in ry))**0.5
    return num/den if den else float('nan')

# structure -> (nesso input name, dG bind, kind).  orig_f12 and s1_design are the
# same sequence: shell 1's flagship. s3_orig_f12/s3_orig_control likewise collapse,
# since they differ only in forced contacts, which nesso never sees.
SET = [
    ("s3_orig_f12",          "s3_design",      -24.33, "design"),
    ("s3_esm2_f4",           "s3_esm2",        -21.08, "design"),
    ("s2_esm2_control",      "s2_esm2_control",-16.25, "design"),
    ("shuffle_control",      "s3_shuffle",     -15.13, "null"),
    ("shuffle_control_esm0", "s3_shuffle_esm", -14.32, "null"),
    ("orig_f12",             "s1_design",      -13.71, "design"),
]
base = os.path.dirname(os.path.abspath(__file__))
rows = []
for struct, inp, dg, kind in SET:
    a = json.load(open(f"{base}/predictions/{inp}/affinity.json"))
    rows.append(dict(structure=struct, kind=kind, dg=dg,
                     nesso=a["affinity_pred_value"],
                     ens=abs(a["affinity_pred_value1"]-a["affinity_pred_value2"]),
                     pbind=a["affinity_probability_binary"],
                     ent=a["entropy_crop_pl"]))

hdr = f"{'structure':22}{'kind':8}{'dG':>8}{'nesso':>8}{'ens':>6}{'P(bind)':>9}{'ent_pl':>8}"
print(hdr); print("-"*len(hdr))
for r in sorted(rows, key=lambda r: r["dg"]):
    print(f"{r['structure']:22}{r['kind']:8}{r['dg']:>8.2f}{r['nesso']:>8.3f}"
          f"{r['ens']:>6.2f}{r['pbind']:>9.3f}{r['ent']:>8.2f}")

dg = [r["dg"] for r in rows]; nv = [r["nesso"] for r in rows]
# dG: more negative = stronger.  nesso: lower = stronger (assumed convention).
# Agreement therefore shows as POSITIVE rho.
rho = spearman(dg, nv)
print(f"\nSpearman(nesso, dG) = {rho:+.2f}   n=6   [+ve = agrees, -ve = ranks backwards]")
print("README comparators vs dG: static sum -0.80, UMA -0.80, enclosed +0.40,")
print("                          MM/GBSA VDWAALS +0.80, mean engagement +0.83 (n=6)")
print("\nby dG rank:   ", [r["structure"] for r in sorted(rows, key=lambda r: r["dg"])])
print("by nesso rank:", [r["structure"] for r in sorted(rows, key=lambda r: r["nesso"])])

des = [r for r in rows if r["kind"] == "design"]; nul = [r for r in rows if r["kind"] == "null"]
print(f"\ndesigns  nesso {min(r['nesso'] for r in des):.2f}..{max(r['nesso'] for r in des):.2f}"
      f"   nulls {min(r['nesso'] for r in nul):.2f}..{max(r['nesso'] for r in nul):.2f}"
      "   -> separable?" ,
      max(r['nesso'] for r in nul) < min(r['nesso'] for r in des)
      or min(r['nesso'] for r in nul) > max(r['nesso'] for r in des))
print(f"mean ensemble spread {sum(r['ens'] for r in rows)/len(rows):.2f} "
      f"vs total nesso range {max(nv)-min(nv):.2f}")
with open(f"{base}/comparison.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print("\nwrote comparison.csv")
