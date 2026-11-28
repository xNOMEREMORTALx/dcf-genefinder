#!/usr/bin/env python3
# DCF L4037 + L4049 grader (separate). Gary NOMEREMORTAL Mason. MIT.
import re, sys
a = re.search(r"SELFTEST planted_detected (\d+)/20 matched_false_pos (\d+)/100 naive_false_pos (\d+)/100", open(sys.argv[1]).read())
pw, fp, nv = map(int, a.groups())
r = [pw >= 18, fp <= 10, nv >= 30]
print(f"L4037 A1 power {pw}/20 (>=18) {'PASS' if r[0] else 'FAIL'}; A2 false pos {fp}/100 (<=10) {'PASS' if r[1] else 'FAIL'}; A3 naive {nv}/100 (>=30) {'PASS' if r[2] else 'FAIL'} -> {sum(r)} of 3")
if len(sys.argv) > 2:
    B = {int(k): (int(p), int(f), float(g)) for k, p, f, g in re.findall(r"SELFTEST2 k (\d+) planted_detected (\d+)/20 matched_false_pos (\d+)/100 mean_zlength_gap ([-\d.]+)", open(sys.argv[2]).read())}
    b = [B[50][1] >= 8, B[10][1] <= 8, B[10][0] >= 16, B[50][2] >= 0.05 and B[10][2] <= 0.5 * B[50][2]]
    print(f"L4049 k50: power {B[50][0]}/20, fp {B[50][1]}/100, gap {B[50][2]:.4f}; k10: power {B[10][0]}/20, fp {B[10][1]}/100, gap {B[10][2]:.4f}")
    for n, v in zip(("B1", "B2", "B3", "B4"), b): print(f"{n} -> {'PASS' if v else 'FAIL'}")
    print(f"RESULT L4049: {sum(b)} of 4 pass")
