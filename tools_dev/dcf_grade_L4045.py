#!/usr/bin/env python3
# DCF L4045 grader (separate). Gary NOMEREMORTAL Mason. MIT.
import re, sys
t = open(sys.argv[1]).read()
gain = float(re.search(r"EDGE sigma_global_px \S+ sigma_row_px \S+ gain ([\d.]+)", t).group(1)); gate = 0.8 <= gain <= 1.2
E = {(a, b): (float(c), float(d)) for a, b, c, d in re.findall(r"EDGE_THR HH (\S+) (\S+) single_nm ([\d.]+) avg10000_nm ([\d.]+)", t)}
neu = {a: (int(b), float(c)) for a, b, c in re.findall(r"NEURON (\S+) peak \S+ fired (\d) fwhm_ms \S+ eqw_ms ([\d.]+)", t)}
S = {(a, b, c): float(d) for a, b, c, d in re.findall(r"SHAPE (\S+) (\S+) (\S+) eqw_ms ([\d.]+)", t)}
V1 = gate and E[("full", "global")][0] > 155.8 and E[("full", "global")][1] > 3
V2 = gate and E[("short", "rolling")][0] < 0.5 * E[("short", "global")][0]
true = neu["slowK"][1] / neu["HH"][1]; rec = S[("slowK", "short", "rolling")] / S[("HH", "short", "rolling")]
V3 = neu["slowK"][0] == 1 and abs(rec / true - 1) <= 0.25
V4 = all(S[(n, "full", "global")] >= 3.0 for n in ("HH", "slowK"))
V5 = "TREMOR_SELFTEST PASS" in t
print(f"gate cell-edge gain {gain:.3f} -> {'PASS' if gate else 'FAIL (V1, V2 void)'}")
print(f"V1 edge single {E[('full','global')][0]:.1f} nm, 10k {E[('full','global')][1]:.3f} nm -> {'PASS' if V1 else 'FAIL'}")
print(f"V2 edge rolling/global short {E[('short','rolling')][0]/E[('short','global')][0]:.3f} -> {'PASS' if V2 else 'FAIL'}")
print(f"V3 width ratio recovered {rec:.3f} vs true {true:.3f} -> {'PASS' if V3 else 'FAIL'}")
print(f"V4 full-exposure widths HH {S[('HH','full','global')]:.3f} ms, slowK {S[('slowK','full','global')]:.3f} ms -> {'PASS' if V4 else 'FAIL'}")
print(f"V5 tremor self-test -> {'PASS' if V5 else 'FAIL'}")
print(f"RESULT L4045: {sum([V1,V2,V3,V4,V5])} of 5 pass")
