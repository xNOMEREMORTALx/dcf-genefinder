#!/usr/bin/env python3
# DCF L4041 grader (separate). Gary NOMEREMORTAL Mason. MIT.
import re, sys
from collections import defaultdict
rows = re.findall(r"RUN e (\S+) seed (\d+) table (\S+) npred (\d+) found ([\d.]+) correct ([\d.]+) exact ([\d.]+) errors \d+ "
                  r"indels (\d+) indels_in_cds (\d+) fs_reported (\S+) fs_hits (\d+) fs_est (\S+) fs_lo (\S+) fs_hi (\S+)", open(sys.argv[1]).read())
R = defaultdict(list)
for r in rows: R[float(r[0])].append(r)
mean = lambda e, i: sum(float(r[i]) for r in R[e]) / len(R[e])
b = R[0.0][0]
print("rate | runs | found | correct | exact start | predicted genes | tables")
for e in sorted(R):
    print(f"{e:g} | {len(R[e])} | {mean(e,4):.2f} | {mean(e,5):.2f} | {mean(e,6):.2f} | {mean(e,3):.1f} | {','.join(r[2] for r in R[e])}")
f0 = float(b[4])
LR0 = f0 >= 75 and b[2] == '11'
LR1 = f0 - mean(1e-4, 4) <= 3
LR2 = f0 - mean(1e-3, 4) <= 15
LR3 = f0 - mean(5e-3, 4) >= 20
LR4 = mean(5e-3, 3) >= 1.3 * float(b[3])
LR5 = all(r[2] == '11' for e in R for r in R[e])
rs = R[1e-3]
inside = sum(int(r[8]) for r in rs); hits = sum(int(r[10]) for r in rs)
try:
    lo = sum(float(r[12]) for r in rs) / len(rs); hi = sum(float(r[13]) for r in rs) / len(rs); est = sum(float(r[11]) for r in rs) / len(rs)
    rec = 100 * hits / max(inside, 1); LR6 = lo <= rec <= hi
    s6 = f"recall {rec:.1f}% ({hits}/{inside}) vs self-estimate {est:.0f}% [{lo:.0f}-{hi:.0f}]"
except ValueError:
    LR6 = False; s6 = "detector printed no estimate"
for k, v, s in (("LR0", LR0, f"baseline found {f0:.2f}%, table {b[2]}"), ("LR1", LR1, f"drop {f0-mean(1e-4,4):.2f}"),
                ("LR2", LR2, f"drop {f0-mean(1e-3,4):.2f}"), ("LR3", LR3, f"drop {f0-mean(5e-3,4):.2f}"),
                ("LR4", LR4, f"{mean(5e-3,3):.1f} vs {b[3]}"), ("LR5", LR5, "code call"), ("LR6", LR6, s6)):
    print(f"{k} {s} -> {'PASS' if v else 'FAIL'}")
print(f"RESULT L4041: {sum([LR0,LR1,LR2,LR3,LR4,LR5,LR6])} of 7 pass")
