#!/usr/bin/env python3
# DCF L4101-L4109 grader (separate). Gary NOMEREMORTAL Mason. MIT.
import re, sys
a, b, c = (open(p).read() for p in sys.argv[1:4]); R = {}
cal = [tuple(map(float, x)) for x in re.findall(r"CAL rf_true ([\d.]+) two_tone_2s \S+ one_tone_2s \S+ one_tone_10s \S+ err_two ([-+\d.]+) err_one ([-+\d.]+) err_long ([-+\d.]+)", a)]
R['T1'] = all(round(abs(x[1]), 9) <= 0.005 for x in cal); R['T2'] = all(round(abs(x[3]), 9) <= 0.005 for x in cal)
R['T3'] = sum(abs(x[1]) for x in cal) < sum(abs(x[2]) for x in cal)
gate = all(v == '1' for v in re.findall(r"membership_ok (\d)", b))
H = re.findall(r"HIT m (\d+) k (\d+) d (\d+) family (\S+) value \S+(?: \S+)? big (\d)", b)
fam3 = {int(k) for m, k, d, f, big in H if m == '3' and d == '3' and f == 'centred_hexagonal'}
R['C1'] = gate and fam3 == set(range(1, 21))
others = [(m, k, d, f) for m, k, d, f, big in H if big == '1' and not (m == '3' and d == '3' and f == 'centred_hexagonal')]
R['C2'] = gate and not others
s = re.search(r"lamps_found_exact_geometry (\d+)/10 lamps_found_occluder_off_2cm (\d+)/10", c)
R['K1'] = int(s.group(1)) >= 8; R['K2'] = int(s.group(2)) >= 6
print("calibration errors (two, one, 10 s):", [(x[0], x[1], x[2], x[3]) for x in cal])
print("gate", gate, "| m3/d3 hex k:", sorted(fam3), "| other big hits:", others)
print("corner:", s.group(1), s.group(2))
for k in R: print(f"{k} -> {'PASS' if R[k] else 'FAIL'}")
for rd, ks in (("L4101", "T1 T2 T3"), ("L4105", "C1 C2"), ("L4109", "K1 K2")):
    print(f"RESULT {rd}: {sum(R[x] for x in ks.split())} of {len(ks.split())} pass")
