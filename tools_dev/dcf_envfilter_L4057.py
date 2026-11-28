#!/usr/bin/env python3
# DCF L4057 - simulated environments (small room, pipe) as filters: what can filtering pull out of audio? Gary NOMEREMORTAL Mason. MIT.
import numpy as np
FS = 8000; C = 343.0; SEED = 4057

def room_rir(room=(3.0, 4.0, 2.5), src=(1.0, 1.0, 1.2), mic=(2.0, 3.0, 1.5), beta=0.7, order=8, length=0.5):
    """Image-source shoebox room impulse response (nearest-sample delays, 1/r spreading, wall reflection beta per bounce)."""
    h = np.zeros(int(length * FS)); L = np.array(room); s = np.array(src); m = np.array(mic)
    d0 = np.linalg.norm(s - m)
    rng = range(-order, order + 1)
    for nx in rng:
        for ny in rng:
            for nz in rng:
                for px in (0, 1):
                    for py in (0, 1):
                        for pz in (0, 1):
                            p = np.array([px, py, pz]); n = np.array([nx, ny, nz])
                            img = (1 - 2 * p) * s + 2 * n * L
                            refl = abs(2 * nx - px) + abs(2 * ny - py) + abs(2 * nz - pz)
                            d = np.linalg.norm(img - m); k = int(round(d / C * FS))
                            if k < len(h): h[k] += beta ** refl * d0 / d
    return h, int(round(d0 / C * FS))

def pipe_H(f, L=1.0, Q=30.0, nmodes=5):
    """Closed-open pipe as a modal filter: odd quarter-wave resonances f_n = (2n-1)c/4L, pressure gain Q/(2n-1) at resonance."""
    H = np.zeros_like(f, dtype=complex); ff = np.maximum(f, 1e-6)
    for n in range(1, nmodes + 1):
        fn = (2 * n - 1) * C / (4 * L)
        H += (Q / (2 * n - 1)) / (1 + 1j * Q * (ff / fn - fn / ff))
    return H

def conv_fft(x, H):
    return np.fft.irfft(np.fft.rfft(x) * H, len(x))

def main():
    rng = np.random.default_rng(SEED)
    N = 2 * FS; t = np.arange(N) / FS; f = np.fft.rfftfreq(N, 1 / FS)
    print("DCF L4057 tool output")
    h, k0 = room_rir(); Hr = np.fft.rfft(h, N)
    E_dir = h[k0] ** 2; E_tot = np.sum(h ** 2)
    print(f"ROOM rir_taps {np.count_nonzero(h)} direct_delay_ms {1000 * k0 / FS:.2f} energy_total_over_direct {E_tot / E_dir:.3f}")
    # A1: dereverberation of speech-like bursts, 20 dB SNR
    s = np.zeros(N)
    for b in range(4):
        a0 = int((0.2 + 0.4 * b) * FS); seg = rng.normal(0, 1, int(0.2 * FS))
        S = np.fft.rfft(seg); fb = np.fft.rfftfreq(len(seg), 1 / FS); S[(fb < 100) | (fb > 3000)] = 0
        s[a0:a0 + len(seg)] = np.fft.irfft(S, len(seg))
    wet = conv_fft(s, Hr); sig_n = np.sqrt(np.mean(wet ** 2) / 10 ** (20 / 10)); y = wet + rng.normal(0, sig_n, N)
    Ps = np.mean(np.abs(np.fft.rfft(s)) ** 2); Pn = sig_n ** 2 * N
    W = np.conj(Hr) / (np.abs(Hr) ** 2 + Pn / Ps)
    s_hat = np.fft.irfft(np.fft.rfft(y) * W, N)
    r_raw = np.corrcoef(np.roll(y, -k0), s)[0, 1] / h[k0] * h[k0]
    r_raw = abs(np.corrcoef(np.roll(y, -k0), s)[0, 1]); r_w = abs(np.corrcoef(s_hat, s)[0, 1])
    print(f"A1 r_raw {r_raw:.4f} r_wiener {r_w:.4f}")
    # A2 + gate: detection of a known chirp in the room (white sensor noise sigma = 1); deflection d = <g, h*s>/(sigma ||g||)
    sc = np.zeros(N); n_c = int(0.5 * FS); tc = np.arange(n_c) / FS
    sc[:n_c] = np.sin(2 * np.pi * (200 * tc + 0.5 * (1000 - 200) / 0.5 * tc ** 2))
    hs = conv_fft(sc, Hr)
    d_aware = np.linalg.norm(hs)
    d_dry = np.dot(sc, hs) / np.linalg.norm(sc)
    d_dry_best = max(np.dot(np.roll(sc, k), hs) / np.linalg.norm(sc) for k in range(0, int(0.05 * FS)))
    Wd = np.conj(Hr) / (np.abs(Hr) ** 2 + 1e-2 * np.mean(np.abs(Hr) ** 2))
    g_w = np.fft.irfft(np.fft.rfft(sc) * np.conj(Wd), N)          # Wiener-then-dry-MF expressed as one template on y
    d_wiener = np.dot(g_w, hs) / np.linalg.norm(g_w)
    db = lambda r: 20 * np.log10(r)
    print(f"A2 d_aware {d_aware:.4f} d_dry_bestlag {d_dry_best:.4f} d_wiener_then_dry {d_wiener:.4f} aware_minus_dry_dB {db(d_aware / d_dry_best):.3f} wiener_minus_aware_dB {db(d_wiener / d_aware):.4f}")
    # A3: pipe as pre-amplifier: source amplitude needed for d = 5 (1 s tone, sensor noise sigma = 1 after the pipe)
    Hp = pipe_H(f); f1 = C / 4.0; f_off = 2 * f1
    def need(freq, H):
        tone = np.sin(2 * np.pi * freq * t[:FS]); x = np.zeros(N); x[:FS] = tone
        return 5.0 / np.linalg.norm(conv_fft(x, H))
    one = np.ones_like(f)
    a_on, a_off, a_free = need(f1, Hp), need(f_off, Hp), need(f1, one)
    print(f"A3 f_on {f1:.2f} f_off {f_off:.2f} need_on {a_on:.6f} need_off {a_off:.6f} need_free {a_free:.6f} "
          f"on_vs_off_dB {db(a_off / a_on):.3f} on_vs_free_dB {db(a_free / a_on):.3f}")
    print("END")

if __name__ == '__main__':
    main()
