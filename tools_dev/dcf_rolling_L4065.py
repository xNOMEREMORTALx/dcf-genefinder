#!/usr/bin/env python3
# DCF L4065 - rolling-shutter auto-calibration and waveform fix. Gary NOMEREMORTAL Mason. MIT.
import os, sys, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import dcf_visualmic as vm
F, ROWS, AMP, NOISE, SECS = 240, 120, 0.3, 2.0, 2.0

def frames(freq, rf, seed):
    rng = np.random.default_rng(seed)
    base = rng.random((ROWS, 320)).astype(np.float32)
    if vm.cv2 is not None: base = vm.cv2.GaussianBlur(base, (0, 0), 2.0)
    xs = np.arange(320, dtype=np.float32)
    for k in range(int(SECS * F)):
        img = np.empty((ROWS, 160), np.float32)
        for r in range(ROWS):
            img[r] = np.interp(np.arange(160) + 40 - AMP * math.sin(2 * math.pi * freq * (k / F + rf * r / (F * ROWS))), xs, base[r])
        yield np.clip(np.round(img * 255 + rng.normal(0, NOISE, img.shape)), 0, 255).astype(np.uint8)

def rows(freq, rf, seed, ref_frames):
    vm.REF_FRAMES = ref_frames; vm.MAX_WIDTH = 320; vm.ROI = None
    return vm.rolling_rows(frames(freq, rf, seed))

def corr(R, freq, rf):
    x, rate, grid, mask = vm.resample_rows(R, F, rf); x = vm.clean(x, rate)
    truth = np.sin(2 * np.pi * freq * grid)
    return abs(np.corrcoef(x[mask], truth[mask])[0, 1]), abs(np.corrcoef(x, truth)[0, 1])

def main():
    print("DCF L4065 tool output")
    for name, rf_true, seed in (("K1", 0.9, 4065), ("K2", 0.75, 4066)):
        R = rows(1530, rf_true, seed, 1)
        print(f"{name} tone 1530 rf_true {rf_true} rf_calibrated {vm.calibrate_row_fraction(R, F, 1530):.3f}")
    R2 = rows(1530, 0.75, 4066, 1); rf_cal = vm.calibrate_row_fraction(R2, F, 1530)
    R3 = rows(1130, 0.75, 4067, 1)
    x, rate, _, _ = vm.resample_rows(R3, F, rf_cal); x = vm.clean(x, rate); fr, X = vm.spectrum(x, rate); sel = fr >= 20
    print(f"K3 second video 1130 Hz analysed with rf {rf_cal:.3f} peak_hz {fr[sel][np.argmax(X[sel])]:.2f}")
    Rm = rows(530, 0.9, 4061, 31); R1 = rows(530, 0.9, 4061, 1)
    cm, cma = corr(Rm, 530, 0.9); c1, c1a = corr(R1, 530, 0.9)
    print(f"W median_ref masked {cm:.4f} unmasked {cma:.4f} | single_ref masked {c1:.4f} unmasked {c1a:.4f}")
    print("END")

if __name__ == '__main__':
    main()
