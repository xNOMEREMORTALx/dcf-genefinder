#!/usr/bin/env python3
# DCF L4061 - recovering what happens BETWEEN frames from rolling-shutter row timing (dcf_visualmic ROLLING mode, first test).
# Gary NOMEREMORTAL Mason. MIT.
import os, sys, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import dcf_visualmic as vm
F, ROWS, RF, AMP, NOISE, SECS = 240, 120, 0.9, 0.3, 2.0, 2.0

def frames(freq, rng):
    base = rng.random((ROWS, 160 * 2)).astype(np.float32)
    if vm.cv2 is not None: base = vm.cv2.GaussianBlur(base, (0, 0), 2.0)
    xs = np.arange(320, dtype=np.float32)
    for k in range(int(SECS * F)):
        img = np.empty((ROWS, 160), np.float32)
        for r in range(ROWS):
            tr = k / F + RF * r / (F * ROWS)
            img[r] = np.interp(np.arange(160) + 40 - AMP * math.sin(2 * math.pi * freq * tr), xs, base[r])
        yield np.clip(np.round(img * 255 + rng.normal(0, NOISE, img.shape)), 0, 255).astype(np.uint8)

def analyse(freq, rolling, rf_assumed, seed):
    vm.ROLLING = rolling; vm.CAPTURE_FPS = F; vm.ROW_FRACTION = rf_assumed; vm.MAX_WIDTH = 320; vm.ROI = None
    x, rate = vm.measure(frames(freq, np.random.default_rng(seed))); x = vm.clean(x, rate)
    fr, X = vm.spectrum(x, rate); sel = fr >= 20; pk = fr[sel][np.argmax(X[sel])]
    return x, rate, float(pk)

def main():
    print("DCF L4061 tool output")
    for name, freq, roll, rf in (("R1", 530, True, 0.9), ("R2", 530, False, 0.9), ("R3", 1530, True, 0.9), ("R4", 1530, True, 0.8)):
        x, rate, pk = analyse(freq, roll, rf, 4061)
        line = f"{name} f {freq} rolling {int(roll)} rf_assumed {rf} rate {rate:.0f} peak_hz {pk:.2f}"
        if name == "R1":
            g = np.arange(len(x)) / rate; truth = np.sin(2 * np.pi * freq * g)
            line += f" corr_with_truth {abs(np.corrcoef(x, truth)[0, 1]):.4f}"
        print(line)
    print("END")

if __name__ == '__main__':
    main()
