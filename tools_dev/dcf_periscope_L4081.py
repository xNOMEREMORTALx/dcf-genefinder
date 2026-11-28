#!/usr/bin/env python3
# DCF L4081 - penumbra ("corner camera") periscopy in 2-D flatland: a diffuse wall seen by an ordinary camera, with and
# without an occluder between the hidden scene and the wall. Gary NOMEREMORTAL Mason. MIT.
import numpy as np
N, M, SEED = 48, 200, 4081
def transport(occluder):
    xs = np.linspace(-1.0, 1.0, N); xm = np.linspace(0.0, 2.0, M)       # hidden scene on y = 2, wall on y = 0
    A = np.zeros((M, N))
    for i, w in enumerate(xm):
        for j, s in enumerate(xs):
            dx, dy = w - s, 2.0; r2 = dx * dx + dy * dy
            cos_w = dy / np.sqrt(r2); cos_s = dy / np.sqrt(r2)
            v = 1.0
            if occluder:
                xc = s + dx * 0.5                                         # where the ray crosses y = 1
                if 0.2 <= xc <= 0.35: v = 0.0
            A[i, j] = v * cos_w * cos_s / r2
    return A / A.max()
def scene(rng):
    x = np.zeros(N); edges = np.sort(rng.choice(np.arange(1, N), 3, replace=False))
    for a, b in zip(np.r_[0, edges], np.r_[edges, N]): x[a:b] = rng.random()
    return x
def recon(A, y, sigma, prior_var):
    lam = sigma ** 2 / prior_var
    return np.linalg.solve(A.T @ A + lam * np.eye(N), A.T @ y + lam * 0.5 * np.ones(N))   # Gaussian prior: mean 0.5, variance 1/12
def main():
    rng = np.random.default_rng(SEED); print("DCF L4081 tool output")
    for occ in (False, True):
        A = transport(occ); sv = np.linalg.svd(A, compute_uv=False)
        print(f"COND occluder {int(occ)} cond {sv[0] / sv[-1]:.4e} effective_rank_1pct {int(np.sum(sv > 0.01 * sv[0]))}")
    for noise in (0.01, 0.10):
        errs = {}
        for occ in (False, True):
            A = transport(occ); r = np.random.default_rng(SEED + int(noise * 1000)); e = []
            for _ in range(50):
                x = scene(r); y0 = A @ x; sig = noise * y0.mean(); y = y0 + r.normal(0, sig, M)
                xh = recon(A, y, sig, 1 / 12)
                e.append(np.linalg.norm(xh - x) / np.linalg.norm(x))
            errs[occ] = float(np.mean(e))
        print(f"ERR noise {noise} no_occluder {errs[False]:.4f} occluder {errs[True]:.4f} ratio {errs[False] / errs[True]:.3f}")
    print("END")
if __name__ == '__main__':
    main()
