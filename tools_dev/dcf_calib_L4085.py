#!/usr/bin/env python3
# DCF L4085 - row-timing calibration: coarse (0.005) vs refined (0.001) scan on four simulated phones. Gary NOMEREMORTAL Mason. MIT.
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import dcf_rolling_L4065 as T, dcf_visualmic as vm
print("DCF L4085 tool output")
for rf, seed in ((0.6, 4085), (0.75, 4086), (0.8, 4087), (0.9, 4088)):
    R = T.rows(1530, rf, seed, 31)
    c = vm.calibrate_row_fraction(R, T.F, 1530); f = vm.calibrate_row_fraction(R, T.F, 1530, refine=True)
    print(f"CAL rf_true {rf} coarse {c:.4f} refined {f:.4f} err_coarse {c - rf:+.4f} err_refined {f - rf:+.4f}")
print("END")
