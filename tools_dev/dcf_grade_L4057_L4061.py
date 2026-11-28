#!/usr/bin/env python3
# DCF L4057 + L4061 grader (separate). Gary NOMEREMORTAL Mason. MIT.
import re, sys
a = open(sys.argv[1]).read(); b = open(sys.argv[2]).read()
r = re.search(r"A1 r_raw ([\d.]+) r_wiener ([\d.]+)", a); A1 = float(r.group(2)) - float(r.group(1)) >= 0.20
m = re.search(r"aware_minus_dry_dB ([-\d.]+) wiener_minus_aware_dB ([-\d.]+)", a)
G = float(m.group(2)) <= 0.01; A2 = float(m.group(1)) >= 3
p = re.search(r"on_vs_off_dB ([-\d.]+) on_vs_free_dB ([-\d.]+)", a); A3 = float(p.group(1)) >= 10 and float(p.group(2)) >= 10
print(f"Gate no-free-lunch: Wiener-then-dry minus aware {m.group(2)} dB -> {'PASS' if G else 'FAIL (void)'}")
print(f"A1 r raw {r.group(1)} -> Wiener {r.group(2)} -> {'PASS' if A1 else 'FAIL'}")
print(f"A2 aware minus dry {m.group(1)} dB -> {'PASS' if A2 else 'FAIL'}")
print(f"A3 pipe on vs off {p.group(1)} dB, on vs free {p.group(2)} dB -> {'PASS' if A3 else 'FAIL'}")
R = {n: (float(pk), c) for n, pk, c in re.findall(r"(R\d) f \d+ rolling \d rf_assumed \S+ rate \d+ peak_hz ([\d.]+)(?: corr_with_truth ([\d.]+))?", b)}
R1 = abs(R['R1'][0] - 530) <= 3 and float(R['R1'][1]) >= 0.8
R2 = abs(R['R2'][0] - 50) <= 3; R3 = abs(R['R3'][0] - 1530) <= 3; R4 = abs(R['R4'][0] - 1530) > 3
for n, v, s in (("R1", R1, f"peak {R['R1'][0]} Hz, corr {R['R1'][1]}"), ("R2", R2, f"peak {R['R2'][0]} Hz"), ("R3", R3, f"peak {R['R3'][0]} Hz"), ("R4", R4, f"peak {R['R4'][0]} Hz")):
    print(f"{n} {s} -> {'PASS' if v else 'FAIL'}")
print(f"RESULT L4057: {sum([A1, A2, A3])} of 3 pass{'' if G else ' (VOID)'}; RESULT L4061: {sum([R1, R2, R3, R4])} of 4 pass")
