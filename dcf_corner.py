#!/usr/bin/env python3
"""dcf_corner.py - a 1-D "corner camera": see what is hidden around a corner from one photo of a wall strip.
Gary NOMEREMORTAL Mason, Digital Cosmos Framework. MIT. numpy only (opencv or Pillow optional, to read the photo). Experimental.

Idea (Bouman et al. 2017; Saunders, Murray-Bruce & Goyal, Nature 2019): an opaque VERTICAL edge (door frame, chair leg) between a
hidden scene and a white wall casts a penumbra on the wall; the wall's brightness along a horizontal strip is a linear function of
the hidden scene, so a single photo can be inverted (flatland model, 1-D hidden scene).

Geometry (metres, top view): wall along y = 0 from WALL_X0 to WALL_X1 (the strip you photograph, left image edge = WALL_X0);
occluder = segment at depth OCC_Y from OCC_X0 to OCC_X1; hidden scene = line at depth HID_Y from HID_X0 to HID_X1 (N patches).
Usage: python3 -I dcf_corner.py photo.jpg --wall 0 2 --occ 0.30 0.35 1.0 --hidden -1 1 2.0 [--rows 0.45 0.55] [--n 48]
       python3 -I dcf_corner.py --selftest
Output: brightness of each hidden patch (left to right) and the brightest positions. Not a measuring instrument.
"""
import sys, argparse
import numpy as np

def transport(wx, hx, occ, hid_y):
    A = np.zeros((len(wx), len(hx)))
    ox0, ox1, oy = occ
    for i, w in enumerate(wx):
        dx = w - hx; r2 = dx * dx + hid_y * hid_y; cosf = hid_y / np.sqrt(r2)
        xc = hx + dx * (1 - oy / hid_y)              # where the ray from patch to wall point crosses depth oy
        vis = ~((xc >= ox0) & (xc <= ox1))
        A[i] = vis * cosf * cosf / r2
    return A / A.max()

def reconstruct(profile, wx, hx, occ, hid_y, rel_noise=0.01):
    A = transport(wx, hx, occ, hid_y)
    A1 = np.hstack([A, np.ones((len(wx), 1))])       # extra column: unknown ambient room light
    y = profile / profile.max(); sig = rel_noise * y.mean(); lam = sig ** 2
    reg = lam * np.eye(A1.shape[1]); reg[-1, -1] = 0
    sol = np.linalg.solve(A1.T @ A1 + reg, A1.T @ y)
    return sol[:-1], sol[-1]

def read_profile(path, rows):
    try:
        import cv2; img = cv2.imread(path, cv2.IMREAD_GRAYSCALE).astype(float)
    except Exception:
        from PIL import Image; img = np.asarray(Image.open(path).convert('L'), float)
    h = img.shape[0]; band = img[int(rows[0] * h):max(int(rows[1] * h), int(rows[0] * h) + 1)]
    return band.mean(axis=0)

def peaks(xh, hx, k=2):
    idx = [i for i in range(1, len(xh) - 1) if xh[i] >= xh[i - 1] and xh[i] >= xh[i + 1]]
    idx = sorted(idx, key=lambda i: -xh[i])[:k]
    return sorted(hx[i] for i in idx)

def selftest(seed=4109, trials=10):
    """L4109: render a synthetic wall photo (two lamps, ambient light, vignetting, noise), reconstruct with exact and with
    slightly wrong geometry; report how often both lamps are found within 1 and 2 patches."""
    rng = np.random.default_rng(seed); N = 48; hx = np.linspace(-1, 1, N); W = 400; wx = np.linspace(0, 2, W)
    occ = (0.30, 0.35, 1.0); hid_y = 2.0; ok_exact = ok_wrong = 0
    for _ in range(trials):
        i = int(rng.integers(4, 20)); j = int(rng.integers(28, 44))
        x = np.zeros(N); x[i] = 1.0; x[j] = 0.7
        y = transport(wx, hx, occ, hid_y) @ x + 0.5            # ambient light
        y = y * (1 - 0.15 * ((wx - 1) / 1) ** 2)                 # vignetting
        y = y + rng.normal(0, 0.01 * y.mean(), W)
        step = hx[1] - hx[0]
        for geom, tol, tag in ((occ, 1, 'exact'), ((0.32, 0.37, 1.0), 2, 'wrong')):
            xh, amb = reconstruct(y, wx, hx, geom, hid_y)
            p = peaks(xh, hx); good = len(p) == 2 and abs(p[0] - hx[i]) <= tol * step + 1e-9 and abs(p[1] - hx[j]) <= tol * step + 1e-9
            if tag == 'exact': ok_exact += good
            else: ok_wrong += good
    print("SELFTEST corner lamps_found_exact_geometry %d/%d lamps_found_occluder_off_2cm %d/%d" % (ok_exact, trials, ok_wrong, trials))
    return ok_exact, ok_wrong

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('photo', nargs='?'); ap.add_argument('--selftest', action='store_true')
    ap.add_argument('--wall', nargs=2, type=float, default=[0, 2]); ap.add_argument('--occ', nargs=3, type=float, default=[0.30, 0.35, 1.0])
    ap.add_argument('--hidden', nargs=3, type=float, default=[-1, 1, 2.0]); ap.add_argument('--rows', nargs=2, type=float, default=[0.45, 0.55])
    ap.add_argument('--n', type=int, default=48)
    a = ap.parse_args()
    if a.selftest or not a.photo:
        selftest(); return
    prof = read_profile(a.photo, a.rows); wx = np.linspace(a.wall[0], a.wall[1], len(prof)); hx = np.linspace(a.hidden[0], a.hidden[1], a.n)
    xh, amb = reconstruct(prof, wx, hx, tuple(a.occ), a.hidden[2])
    print("hidden scene brightness, left to right (relative):"); print(" ".join("%.2f" % v for v in xh / max(xh.max(), 1e-12)))
    print("brightest positions (m):", ", ".join("%.2f" % v for v in peaks(xh, hx)), "| ambient %.3f" % amb)

if __name__ == '__main__':
    main()
