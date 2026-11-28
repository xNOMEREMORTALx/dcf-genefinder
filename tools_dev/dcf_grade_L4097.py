#!/usr/bin/env python3
# DCF L4097 grader (separate). Gary NOMEREMORTAL Mason. MIT.
import re, sys
c = [tuple(map(float, x)) for x in re.findall(r"CAL rf_true ([\d.]+) energy ([\d.]+) sideband ([\d.]+) err_energy ([-+\d.]+) err_sideband ([-+\d.]+)", open(sys.argv[1]).read())]
S1 = all(round(abs(x[4]), 9) <= 0.005 for x in c); S2 = round(abs(sum(x[4] for x in c) / len(c)), 9) <= 0.003
S3 = sum(round(abs(x[4]), 9) <= round(abs(x[3]), 9) for x in c) >= 3
for x in c: print(f"true {x[0]}: energy err {x[3]:+.4f}, sideband err {x[4]:+.4f}")
for n, v in (("S1", S1), ("S2", S2), ("S3", S3)): print(f"{n} -> {'PASS' if v else 'FAIL'}")
print(f"RESULT L4097: {sum([S1, S2, S3])} of 3 pass")
