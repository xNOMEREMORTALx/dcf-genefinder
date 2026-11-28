#!/usr/bin/env python3
# DCF L4029 grader (separate from the tool). Gary NOMEREMORTAL Mason. MIT.
import re, sys
txt = open(sys.argv[1]).read()
g = re.search(r"G1 gain ([\d.]+) resid_px [\d.]+ ratio_to_sigma ([\d.]+)", txt)
G1 = 0.9 <= float(g.group(1)) <= 1.1 and abs(float(g.group(2)) - 1) <= 0.2
neu = {m[0]: (int(m[1]), float(m[2])) for m in re.findall(r"NEURON (\S+) rest \S+ peak \S+ fired (\d) halfwidth_ms ([\d.]+)", txt)}
G2 = neu["HH"][0] == 1 and neu["fastK"][0] == 1 and neu["fastK"][1] < neu["HH"][1]
print(f"G1 pipeline gain {g.group(1)}, noise ratio {g.group(2)} -> {'PASS' if G1 else 'FAIL'}")
print(f"G2 neurons fire; halfwidth HH {neu['HH'][1]} ms, fastK {neu['fastK'][1]} ms -> {'PASS' if G2 else 'FAIL'}")
thr = {}
for m in re.findall(r"THR (\S+) (\S+) (\S+) single_nm ([\d.]+) avg1000_nm ([\d.]+) avg10000_nm ([\d.]+) pdet_1nm ([\d.]+)", txt):
    thr[(m[0], m[1], m[2])] = tuple(float(x) for x in m[3:])
names = ("HH", "fastK")
N1 = all(thr[(n, "full", "global")][3] <= 0.05 for n in names)
N2 = all(20 <= thr[(n, "full", "global")][0] <= 1000 for n in names)
N3 = all(thr[(n, "full", "global")][2] > 3 for n in names)
sta = {(a, b): float(c) for a, b, c in re.findall(r"STA (\S+) (\S+) global amp_x_thr 10 halfwidth_ms ([\d.]+)", txt)}
N4 = all(sta[(n, "full")] >= 4.0 for n in names)
true_ratio = neu["HH"][1] / neu["fastK"][1]
rec_ratio = sta[("HH", "short")] / sta[("fastK", "short")]
N5 = abs(rec_ratio / true_ratio - 1) <= 0.25
N6 = all(0.5 <= thr[(n, "short", "rolling")][0] / thr[(n, "short", "global")][0] <= 2 for n in names)
for n in names:
    print(f"  {n}: single full-global {thr[(n,'full','global')][0]:.1f} nm; 10k-avg {thr[(n,'full','global')][2]:.3f} nm; "
          f"pdet(1nm) {thr[(n,'full','global')][3]:.4f}; rolling/global short {thr[(n,'short','rolling')][0]/thr[(n,'short','global')][0]:.2f}")
print(f"N1 single 1 nm detection <= 5% -> {'PASS' if N1 else 'FAIL'}")
print(f"N2 single-spike threshold 20-1000 nm -> {'PASS' if N2 else 'FAIL'}")
print(f"N3 10k-spike threshold > 3 nm -> {'PASS' if N3 else 'FAIL'}")
print(f"N4 full-exposure STA halfwidth >= 4 ms both: HH {sta[('HH','full')]} fastK {sta[('fastK','full')]} -> {'PASS' if N4 else 'FAIL'}")
print(f"N5 short-exposure ratio recovered {rec_ratio:.3f} vs true {true_ratio:.3f} (within 25%) -> {'PASS' if N5 else 'FAIL'}")
print(f"N6 rolling vs global within 2x -> {'PASS' if N6 else 'FAIL'}")
k = sum([N1, N2, N3, N4, N5, N6])
print(f"RESULT L4029: {k} of 6 pass" + ("" if G1 and G2 else " (VOID: gate failed)"))
