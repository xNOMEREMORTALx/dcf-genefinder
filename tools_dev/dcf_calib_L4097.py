#!/usr/bin/env python3
# DCF L4097 - calibration criterion: tone energy (L4065) vs sideband leakage, four NEW simulated phones. Gary NOMEREMORTAL Mason. MIT.
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import dcf_rolling_L4065 as T, dcf_visualmic as vm
print("DCF L4097 tool output")
for rf, seed in ((0.65, 4097), (0.7, 4098), (0.85, 4099), (0.95, 4100)):
    R = T.rows(1530, rf, seed, 31)
    e = vm.calibrate_row_fraction(R, T.F, 1530, refine=True, criterion="energy")
    s = vm.calibrate_row_fraction(R, T.F, 1530, refine=True, criterion="sideband")
    print(f"CAL rf_true {rf} energy {e:.4f} sideband {s:.4f} err_energy {e - rf:+.4f} err_sideband {s - rf:+.4f}")
print("END")
