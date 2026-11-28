#!/usr/bin/env python3
# DCF L4101 - two-tone vs one-tone calibration, and a longer clip. Gary NOMEREMORTAL Mason. MIT.
import os, sys, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import dcf_visualmic as vm
F, ROWS, NOISE = 240, 120, 2.0
def frames(freqs, amps, rf, secs, seed):
    rng = np.random.default_rng(seed); base = rng.random((ROWS, 320)).astype(np.float32)
    if vm.cv2 is not None: base = vm.cv2.GaussianBlur(base, (0, 0), 2.0)
    xs = np.arange(320, dtype=np.float32)
    for k in range(int(secs * F)):
        img = np.empty((ROWS, 160), np.float32)
        for r in range(ROWS):
            t = k / F + rf * r / (F * ROWS)
            img[r] = np.interp(np.arange(160) + 40 - sum(a * math.sin(2 * math.pi * f * t) for f, a in zip(freqs, amps)), xs, base[r])
        yield np.clip(np.round(img * 255 + rng.normal(0, NOISE, img.shape)), 0, 255).astype(np.uint8)
def rows(freqs, amps, rf, secs, seed):
    vm.REF_FRAMES = 31; vm.MAX_WIDTH = 320; vm.ROI = None
    return vm.rolling_rows(frames(freqs, amps, rf, secs, seed))
print("DCF L4101 tool output")
for rf, seed in ((0.62, 4101), (0.78, 4102), (0.88, 4103), (0.93, 4104)):
    R2 = rows([530, 1530], [0.2, 0.2], rf, 2.0, seed)
    two = vm.calibrate_row_fraction(R2, F, [530, 1530], refine=True)
    one = vm.calibrate_row_fraction(R2, F, 1530, refine=True)
    R10 = rows([1530], [0.3], rf, 10.0, seed + 100)
    lng = vm.calibrate_row_fraction(R10, F, 1530, refine=True)
    print(f"CAL rf_true {rf} two_tone_2s {two:.4f} one_tone_2s {one:.4f} one_tone_10s {lng:.4f} "
          f"err_two {two - rf:+.4f} err_one {one - rf:+.4f} err_long {lng - rf:+.4f}")
print("END")
