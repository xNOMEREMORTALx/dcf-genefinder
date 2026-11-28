#!/usr/bin/env python3
# DCF L4029 - can the visual microphone see nerve impulses? (Gary NOMEREMORTAL Mason). MIT. Not for medical use.
# Pre-registration: preregistrations/PREREGISTERED_L4029_neuron_visualmic.txt
import os, sys, math
import numpy as np
from math import erf, sqrt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dcf_visualmic as vm

SEED = 4029
F = 240.0                 # frames per second
PIX_NM = 100.0            # nm of sample per pixel
NOISE_DN = 2.0
H, W = 120, 160
ROWFRAC = 0.9
EXPO = {"full": 1.0 / F, "short": 0.25e-3}
Z99 = 2.3263478740408408  # one-sided 1% false alarm
NPHASE = 200

def Phi(x): return 0.5 * (1 + erf(x / sqrt(2)))

# ---------------- neuron models (Hodgkin-Huxley, mV, ms) ----------------
def hh_spike(kfac, T=20.0, dt=0.002, tmax=60.0, istim=40.0, tstim=(5.0, 5.5)):
    phi = 3 ** ((T - 6.3) / 10)
    gNa, gK, gL, ENa, EK, EL = 120.0, 36.0, 0.3, 50.0, -77.0, -54.387
    def an(V): return 0.01 * (V + 55) / (1 - math.exp(-(V + 55) / 10)) if abs(V + 55) > 1e-7 else 0.1
    def bn(V): return 0.125 * math.exp(-(V + 65) / 80)
    def am(V): return 0.1 * (V + 40) / (1 - math.exp(-(V + 40) / 10)) if abs(V + 40) > 1e-7 else 1.0
    def bm(V): return 4 * math.exp(-(V + 65) / 18)
    def ah(V): return 0.07 * math.exp(-(V + 65) / 20)
    def bh(V): return 1 / (1 + math.exp(-(V + 35) / 10))
    V = -65.0
    m = am(V) / (am(V) + bm(V)); h = ah(V) / (ah(V) + bh(V)); n = an(V) / (an(V) + bn(V))
    # settle to rest with this variant
    for _ in range(int(50 / dt)):
        INa = gNa * m**3 * h * (V - ENa); IK = gK * n**4 * (V - EK); IL = gL * (V - EL)
        V += dt * (-(INa + IK + IL))
        m += dt * phi * (am(V) * (1 - m) - bm(V) * m)
        h += dt * phi * (ah(V) * (1 - h) - bh(V) * h)
        n += dt * phi * kfac * (an(V) * (1 - n) - bn(V) * n)
    vrest = V
    ts = np.arange(0, tmax, dt); vs = np.empty_like(ts)
    for i, t in enumerate(ts):
        I = istim if tstim[0] <= t < tstim[1] else 0.0
        INa = gNa * m**3 * h * (V - ENa); IK = gK * n**4 * (V - EK); IL = gL * (V - EL)
        V += dt * (I - (INa + IK + IL))
        m += dt * phi * (am(V) * (1 - m) - bm(V) * m)
        h += dt * phi * (ah(V) * (1 - h) - bh(V) * h)
        n += dt * phi * kfac * (an(V) * (1 - n) - bn(V) * n)
        vs[i] = V
    return ts / 1000.0, vs, vrest     # seconds

def fwhm(t, y):
    p = np.argmax(y); half = y[p] / 2
    i = p
    while i > 0 and y[i] > half: i -= 1
    j = p
    while j < len(y) - 1 and y[j] > half: j += 1
    return t[j] - t[i]

# ---------------- camera noise calibration with the real MotionMeter ----------------
def texture(rng):
    base = rng.random((H, W * 2)).astype(np.float32)
    if vm.cv2 is not None:
        base = vm.cv2.GaussianBlur(base, (0, 0), 2.0)
    return base

def frame(base, d_px, rng, noise=NOISE_DN):
    xs = np.arange(W * 2, dtype=np.float32)
    img = np.empty((H, W), np.float32)
    for r in range(H):
        img[r] = np.interp(np.arange(W) + 40 - d_px, xs, base[r])
    img = img * 255 + rng.normal(0, noise, img.shape)
    return np.clip(np.round(img), 0, 255).astype(np.uint8)

def calibrate(rng):
    base = texture(rng)
    ref = vm.to_gray_small(frame(base, 0.0, rng))
    meter = vm.MotionMeter(ref)
    gx, rows = [], []
    for _ in range(300):
        g = vm.to_gray_small(frame(base, 0.0, rng))
        gx.append(meter.shift(g)[0]); rows.append(meter.row_shift(g))
    sig_g = float(np.std(gx)); sig_r = float(np.median(np.std(np.array(rows), axis=0)))
    # G1: known 0.3 px, 7 Hz motion, 400 frames
    k = np.arange(400); truth = 0.3 * np.sin(2 * np.pi * 7 * k / F)
    est = np.array([meter.shift(vm.to_gray_small(frame(base, d, rng)))[0] for d in truth])
    est = est - est.mean()
    gain = float(np.dot(est, truth) / np.dot(truth, truth))
    resid = float(np.std(est - gain * truth))
    return sig_g, sig_r, gain, resid, vm.to_gray_small(ref).shape[0]

# ---------------- sampling of the deformation by the camera ----------------
class Template:
    def __init__(self, t, y):
        self.t = t; self.y = y
        self.c = np.concatenate([[0], np.cumsum((y[1:] + y[:-1]) / 2 * np.diff(t))])
    def boxavg(self, a, b):   # mean of y over [a, b] (zero outside)
        ca = np.interp(a, self.t, self.c, left=0, right=self.c[-1])
        cb = np.interp(b, self.t, self.c, left=0, right=self.c[-1])
        return (cb - ca) / (b - a)

def sample_times(mode, nrows, phase, t0, t1):
    ks = np.arange(math.floor((t0 - phase) * F) - 1, math.ceil((t1 - phase) * F) + 1)
    if mode == "global":
        return (ks / F + phase)[:, None]
    r = np.arange(nrows)[None, :]
    return ks[:, None] / F + phase + ROWFRAC * r / (F * nrows)

def norm2(tpl, mode, expo, nrows, phase):
    st = sample_times(mode, nrows, phase, tpl.t[0], tpl.t[-1])
    s = tpl.boxavg(st, st + expo)
    return float(np.sum(s ** 2))

def threshold_px(tpl, mode, expo, sigma, nrows, nspk, rng):
    phases = rng.random(NPHASE) / F
    n2 = np.array([norm2(tpl, mode, expo, nrows, p) for p in phases])
    if nspk == 1:
        def pdet(A): return np.mean([Phi(A * sqrt(v) / sigma - Z99) for v in n2])
        lo, hi = 1e-6, 1e4
        for _ in range(100):
            mid = sqrt(lo * hi)
            lo, hi = (mid, hi) if pdet(mid) < 0.5 else (lo, mid)
        return sqrt(lo * hi), pdet
    return Z99 * sigma / sqrt(nspk * n2.mean()), None

def sta_halfwidth(tpl, mode, expo, sigma, nrows, nspk, A, rng, binw=0.05e-3, smooth_bins=5):
    edges = np.arange(tpl.t[0], tpl.t[-1] + binw, binw)
    acc = np.zeros(len(edges) - 1); cnt = np.zeros(len(edges) - 1)
    for _ in range(nspk):
        ph = rng.random() / F
        st = sample_times(mode, nrows, ph, tpl.t[0], tpl.t[-1])
        y = A * tpl.boxavg(st, st + expo) + rng.normal(0, sigma, st.shape)
        tc = (st + expo / 2).ravel()
        idx = np.digitize(tc, edges) - 1
        ok = (idx >= 0) & (idx < len(acc))
        np.add.at(acc, idx[ok], y.ravel()[ok]); np.add.at(cnt, idx[ok], 1)
    prof = np.where(cnt > 0, acc / np.maximum(cnt, 1), 0.0)
    prof = np.convolve(prof, np.ones(smooth_bins) / smooth_bins, "same")
    return fwhm((edges[:-1] + edges[1:]) / 2, prof)

def main():
    rng = np.random.default_rng(SEED)
    print("DCF L4029 tool output - neuron visual microphone feasibility")
    sig_g, sig_r, gain, resid, nrows = calibrate(rng)
    print(f"CALIB sigma_global_px {sig_g:.5f} sigma_row_px {sig_r:.5f} rows {nrows}")
    print(f"G1 gain {gain:.4f} resid_px {resid:.5f} ratio_to_sigma {resid / sig_g:.4f}")
    tpls = {}
    for name, kf in (("HH", 1.0), ("fastK", 4.0)):
        t, v, vr = hh_spike(kf)
        vp = v.max(); y = (v - vr) / (vp - vr)
        hw = fwhm(t, v - vr)
        print(f"NEURON {name} rest {vr:.2f} peak {vp:.2f} fired {int(vp > 0)} halfwidth_ms {hw * 1e3:.4f}")
        tpls[name] = Template(t, y)
    res = {}
    for name, tpl in tpls.items():
        for expo_name, expo in EXPO.items():
            for mode in ("global", "rolling"):
                sigma = sig_g if mode == "global" else sig_r
                thr1, pdet = threshold_px(tpl, mode, expo, sigma, nrows, 1, rng)
                thr1k, _ = threshold_px(tpl, mode, expo, sigma, nrows, 1000, rng)
                thr10k, _ = threshold_px(tpl, mode, expo, sigma, nrows, 10000, rng)
                p1nm = pdet(1.0 / PIX_NM)
                res[(name, expo_name, mode)] = (thr1, thr1k, thr10k)
                print(f"THR {name} {expo_name} {mode} single_nm {thr1 * PIX_NM:.3f} avg1000_nm {thr1k * PIX_NM:.4f} "
                      f"avg10000_nm {thr10k * PIX_NM:.4f} pdet_1nm {p1nm:.5f} pdet_3nm {pdet(3.0 / PIX_NM):.5f}")
    for expo_name, expo in EXPO.items():
        hws = {}
        for name, tpl in tpls.items():
            thr1k = res[(name, expo_name, "global")][1]
            hws[name] = sta_halfwidth(tpl, "global", expo, sig_g, nrows, 1000, 10 * thr1k, rng)
            print(f"STA {name} {expo_name} global amp_x_thr 10 halfwidth_ms {hws[name] * 1e3:.4f}")
        print(f"STA_RATIO {expo_name} HH/fastK {hws['HH'] / hws['fastK']:.4f}")
    # descriptive: noise needed for a single 3 nm spike at 50% (global, short exposure, fastest variant HH)
    thr1 = res[("HH", "short", "global")][0]
    need = sig_g * (3.0 / PIX_NM) / thr1
    print(f"NEED sigma_global_px_for_single_3nm {need:.6f} (now {sig_g:.5f}; factor {sig_g / need:.1f})")
    print("END")

if __name__ == "__main__":
    main()
