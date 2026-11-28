#!/usr/bin/env python3
# DCF L4065-L4081 grader (separate). Amendment a (parser only): float tolerance at the inclusive +-0.02 boundary (0.73 - 0.75 = 0.0200000000000000018 in floating point). argv: L4065, L4069, L4073, L4081 outputs. Gary NOMEREMORTAL Mason. MIT.
import re, sys
o = [open(p).read() for p in sys.argv[1:5]]; res = {}
k = {n: (float(t), float(c)) for n, t, c in re.findall(r"(K\d) tone 1530 rf_true ([\d.]+) rf_calibrated ([\d.]+)", o[0])}
k3 = float(re.search(r"K3 .* peak_hz ([\d.]+)", o[0]).group(1))
w = re.search(r"W median_ref masked ([\d.]+) unmasked ([\d.]+) \| single_ref masked ([\d.]+) unmasked ([\d.]+)", o[0])
res['K1'] = round(abs(k['K1'][1] - 0.9), 9) <= 0.02; res['K2'] = round(abs(k['K2'][1] - 0.75), 9) <= 0.02; res['K3'] = abs(k3 - 1130) <= 3
res['W1'] = float(w.group(1)) >= 0.8; res['W2'] = float(w.group(3)) < 0.8
d = re.search(r"false_pos_upper (\d+) false_pos_lower (\d+) mean_zlength_gap ([-\d.]+)", o[1])
res['D1'] = int(d.group(1)) <= 10; res['D2'] = int(d.group(2)) >= 16; res['D3'] = float(d.group(3)) >= 0.05
F = re.findall(r"FS e (\S+) seed \d+ ran (\d) reported (\d+) inside_cds (\d+) hits150 (\d+) reported_true150 (\d+) est (\S+) lo (\S+) hi (\S+)", o[2])
res['F1'] = all(r[1] == '1' for r in F)
E = [r for r in F if float(r[0]) > 0]; Z = [r for r in F if float(r[0]) == 0]
ins = sum(int(r[3]) for r in E); hit = sum(int(r[4]) for r in E); rep = sum(int(r[2]) for r in E); tr = sum(int(r[5]) for r in E)
try:
    lo = sum(float(r[7]) for r in E) / len(E); hi = sum(float(r[8]) for r in E) / len(E); rec = 100 * hit / max(ins, 1)
    res['F2'] = res['F1'] and lo <= rec <= hi; f2s = f"recall {rec:.1f}% ({hit}/{ins}) vs interval {lo:.0f}-{hi:.0f}%"
except ValueError:
    res['F2'] = False; f2s = "no estimate"
res['F3'] = res['F1'] and rep > 0 and tr / rep >= 0.5; res['F4'] = res['F1'] and all(int(r[2]) <= 5 for r in Z)
c = {int(a): float(b) for a, b in re.findall(r"COND occluder (\d) cond ([\d.e+]+)", o[3])}
er = {float(n): (float(a), float(b), float(r)) for n, a, b, r in re.findall(r"ERR noise ([\d.]+) no_occluder ([\d.]+) occluder ([\d.]+) ratio ([\d.]+)", o[3])}
res['P1'] = c[0] / c[1] >= 100; res['P2'] = er[0.01][2] >= 2; res['P3'] = er[0.01][1] < 0.5; res['P4'] = er[0.1][2] >= 1.5
print(f"L4065: K1 {k['K1'][1]} K2 {k['K2'][1]} K3 {k3} Hz; W median+mask {w.group(1)}, single+mask {w.group(3)}")
print(f"L4069: upper {d.group(1)}/200, lower {d.group(2)}/200, gap {d.group(3)}")
print(f"L4073: {f2s}; precision {tr}/{rep}; error-free reports {[int(r[2]) for r in Z]}")
print(f"L4081: cond ratio {c[0]/c[1]:.3g}; 1% err {er[0.01][0]:.3f} -> {er[0.01][1]:.3f}; 10% ratio {er[0.1][2]:.2f}")
for key in res: print(f"{key} -> {'PASS' if res[key] else 'FAIL'}")
for rd, keys in (("L4065", "K1 K2 K3 W1 W2"), ("L4069", "D1 D2 D3"), ("L4073", "F1 F2 F3 F4"), ("L4081", "P1 P2 P3 P4")):
    ks = keys.split(); print(f"RESULT {rd}: {sum(res[x] for x in ks)} of {len(ks)} pass")
