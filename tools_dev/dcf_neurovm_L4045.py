#!/usr/bin/env python3
# DCF L4045 - neuron visual microphone, round 2: realistic cell-edge scene, K+ channel speed from spike shape, tremor mode.
# Reuses the frozen L4029 module unchanged. Gary NOMEREMORTAL Mason. MIT. Not for medical use.
import os, sys, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import importlib.util
spec = importlib.util.spec_from_file_location('L4029', os.path.join(HERE, 'dcf_neurovm_L4029.py'))
M = importlib.util.module_from_spec(spec); spec.loader.exec_module(M)
vm = M.vm
SEED = 4045

def edge_base(rng):      # one smooth vertical membrane edge, 10 DN contrast over ~3 px, flat elsewhere (values in 0-1 units)
    x = np.arange(M.W * 2, dtype=np.float32)
    prof = 0.5 + (10 / 255) / (1 + np.exp(-(x - (M.W + 0)) / 1.0))
    return np.tile(prof, (M.H, 1)).astype(np.float32)

def calibrate(base, rng):
    ref = vm.to_gray_small(M.frame(base, 0.0, rng)); meter = vm.MotionMeter(ref)
    gx, rows = [], []
    for _ in range(300):
        g = vm.to_gray_small(M.frame(base, 0.0, rng)); gx.append(meter.shift(g)[0]); rows.append(meter.row_shift(g))
    k = np.arange(400); truth = 0.3 * np.sin(2 * np.pi * 7 * k / M.F)
    est = np.array([meter.shift(vm.to_gray_small(M.frame(base, d, rng)))[0] for d in truth]); est -= est.mean()
    gain = float(np.dot(est, truth) / np.dot(truth, truth))
    return float(np.std(gx)), float(np.median(np.std(np.array(rows), axis=0))), gain

def eq_width(t, y):      # equivalent rectangular width of the main lobe: area of y>0 within +-2 ms of the peak / peak height
    p = np.argmax(y); w = (t >= t[p] - 2e-3) & (t <= t[p] + 2e-3); yy = np.clip(y[w], 0, None)
    return float(np.trapezoid(yy, t[w]) / y[p]) if hasattr(np, 'trapezoid') else float(np.trapz(yy, t[w]) / y[p])

def sta_profile(tpl, mode, expo, sigma, nrows, nspk, A, rng, binw=0.05e-3, smooth_bins=5):
    edges = np.arange(tpl.t[0], tpl.t[-1] + binw, binw); acc = np.zeros(len(edges) - 1); cnt = np.zeros(len(edges) - 1)
    for _ in range(nspk):
        ph = rng.random() / M.F; st = M.sample_times(mode, nrows, ph, tpl.t[0], tpl.t[-1])
        y = A * tpl.boxavg(st, st + expo) + rng.normal(0, sigma, st.shape); tc = (st + expo / 2).ravel()
        idx = np.digitize(tc, edges) - 1; ok = (idx >= 0) & (idx < len(acc))
        np.add.at(acc, idx[ok], y.ravel()[ok]); np.add.at(cnt, idx[ok], 1)
    prof = np.convolve(np.where(cnt > 0, acc / np.maximum(cnt, 1), 0.0), np.ones(smooth_bins) / smooth_bins, 'same')
    return (edges[:-1] + edges[1:]) / 2, prof

def main():
    rng = np.random.default_rng(SEED)
    print("DCF L4045 tool output")
    tpls = {}
    for name, kf in (("HH", 1.0), ("slowK", 0.25)):
        t, v, vr = M.hh_spike(kf); y = (v - vr) / (v.max() - vr)
        tpls[name] = M.Template(t, y)
        print(f"NEURON {name} peak {v.max():.2f} fired {int(v.max() > 0)} fwhm_ms {M.fwhm(t, v - vr) * 1e3:.4f} eqw_ms {eq_width(t, y) * 1e3:.4f}")
    # V1, V2: cell-edge scene
    sg, sr, gain = calibrate(edge_base(rng), rng)
    print(f"EDGE sigma_global_px {sg:.5f} sigma_row_px {sr:.5f} gain {gain:.4f}")
    tp = tpls["HH"]
    r = {}
    for ex in ("full", "short"):
        for mode in ("global", "rolling"):
            s = sg if mode == "global" else sr
            t1, _ = M.threshold_px(tp, mode, M.EXPO[ex], s, 120, 1, rng); t10k, _ = M.threshold_px(tp, mode, M.EXPO[ex], s, 120, 10000, rng)
            r[(ex, mode)] = (t1, t10k)
            print(f"EDGE_THR HH {ex} {mode} single_nm {t1 * M.PIX_NM:.3f} avg10000_nm {t10k * M.PIX_NM:.4f}")
    # V3, V4: K+ speed from shape on the textured scene (L4029 calibration recomputed with seed 4045)
    sgt, srt, _, _, nrows = M.calibrate(rng)
    print(f"TEXTURE sigma_global_px {sgt:.5f} sigma_row_px {srt:.5f}")
    for ex, mode in (("short", "rolling"), ("full", "global")):
        ws = {}
        for name, tpl in tpls.items():
            s = sgt if mode == "global" else srt
            thr1k, _ = M.threshold_px(tpl, mode, M.EXPO[ex], s, nrows, 1000, rng)
            tt, prof = sta_profile(tpl, mode, M.EXPO[ex], s, nrows, 1000, 10 * thr1k, rng)
            ws[name] = eq_width(tt, prof)
            print(f"SHAPE {name} {ex} {mode} eqw_ms {ws[name] * 1e3:.4f}")
        print(f"SHAPE_RATIO {ex} {mode} slowK/HH {ws['slowK'] / ws['HH']:.4f}")
    # V5: tremor self-test of the phone tool
    vm.MODE = "tremorselftest"
    ok = vm.run_tremor("tremorselftest")
    print(f"TREMOR_SELFTEST {'PASS' if ok else 'FAIL'}")
    print("END")

if __name__ == '__main__':
    main()
