#!/usr/bin/env python3
# DCF L4085-L4093 grader (separate). Gary NOMEREMORTAL Mason. MIT.
import re, sys, math
a, b, c = (open(p).read() for p in sys.argv[1:4]); R = {}
cal = [tuple(map(float, x)) for x in re.findall(r"CAL rf_true ([\d.]+) coarse ([\d.]+) refined ([\d.]+) err_coarse ([-+\d.]+) err_refined ([-+\d.]+)", a)]
R['B1'] = all(round(abs(x[4]), 9) <= 0.005 for x in cal); R['B2'] = round(abs(sum(x[4] for x in cal) / len(cal)), 9) <= 0.003
R['B3'] = sum(round(abs(x[4]), 9) <= round(abs(x[3]), 9) for x in cal) >= 3
S = {(float(n), int(o)): int(s) for n, o, s in re.findall(r"SMIN noise ([\d.]+) occluder (\d) s_min (\d+)", b)}
R['R1'] = S[(0.01, 1)] <= 8; R['R2'] = S[(0.01, 0)] > 20; R['R3'] = S[(0.1, 1)] <= 16
def pc(k, l):
    s, t = 0.0, math.exp(-l)
    for i in range(k + 1):
        if i: t *= l / i
        s += t
    return s
def interval(l):
    lo = 0
    while pc(lo, l) < 0.025: lo += 1
    hi = lo
    while pc(hi, l) < 0.975: hi += 1
    return lo, hi
gate = all(x == '1' for x in re.findall(r"GATE \S+ membership_exact (\d)", c))
rows = {}
for line in c.splitlines():
    m = re.match(r"N (\d+) size", line)
    if m and 25 <= int(m.group(1)) <= 40:
        for k, o, e in re.findall(r"(tri|star|hex) O (\d+) E ([\d.]+)", line):
            rows.setdefault(k, [0, 0.0]); rows[k][0] += int(o); rows[k][1] += float(e)
for k, q in (('star', 'Q1'), ('hex', 'Q2')):
    lo, hi = interval(rows[k][1]); R[q] = gate and lo <= rows[k][0] <= hi
    print(f"{q} {k}: O {rows[k][0]} E {rows[k][1]:.3f} interval [{lo},{hi}]")
R['Q3'] = gate and 0.5 <= rows['tri'][1] / 5.542 <= 2; print(f"Q3 tri E {rows['tri'][1]:.3f} vs L4001 5.542 (O {rows['tri'][0]})")
print("calibration:", cal); print("s_min:", S); print("gate", gate)
for k in R: print(f"{k} -> {'PASS' if R[k] else 'FAIL'}")
for rd, ks in (("L4085", "B1 B2 B3"), ("L4089", "R1 R2 R3"), ("L4093", "Q1 Q2 Q3")):
    print(f"RESULT {rd}: {sum(R[x] for x in ks.split())} of 3 pass")
