"""
dcf_visualmic.py  -  DCF Visual Microphone for Pydroid 3 (Android) or any Python 3
Gary NOMEREMORTAL Mason, Digital Cosmos Framework, 24 November 2026

Turns a silent video of a vibrating object (crisp packet, leaf, paper, foil) into
a sound file and a spectrum. A simplified version of the "visual microphone"
(Davis, Rubinstein et al., SIGGRAPH 2014): it measures the object's sub-pixel
motion in every frame and treats that motion as the audio signal.

HOW TO USE (Pydroid 3)
  1. Install once (Pydroid menu > Pip):  opencv-python  numpy   (matplotlib optional)
  2. Film a light object close to a speaker playing a steady tone, phone on a
     stand (not hand-held), good light, object filling much of the frame.
  3. Edit the SETTINGS block below (VIDEO_PATH, CAPTURE_FPS, TEST_TONE_HZ).
  4. Press Run. Output goes to OUT_DIR: a .wav file, a spectrum .png (if
     matplotlib is installed) and a text report.
  5. With MODE = "selftest" no camera is needed: the script makes a synthetic
     vibrating video and checks it recovers the known frequency.

THE ONE RULE THAT DECIDES SUCCESS (Nyquist)
  A camera taking F frames per second can only record vibrations below F/2 Hz.
     30 fps  -> below 15 Hz      240 fps -> below 120 Hz      960 fps -> below 480 Hz
  So with 240 fps slow motion, play a LOW tone (e.g. 60-100 Hz), not 440 Hz.
  Slow-motion files are often saved re-timed to 30 fps playback, so the file
  header can say 30 even though the camera shot 240. Set CAPTURE_FPS to the real
  shooting rate.
  ROLLING = True tries the rolling-shutter trick (each image row is exposed at a
  slightly different moment), which can reach higher frequencies but needs a
  guess of the row timing (ROW_FRACTION); treat it as experimental.
"""

# ============================== SETTINGS ===================================
MODE = "video"            # "video" | "selftest" | "tremor" (3-15 Hz body/hand motion) | "tremorselftest"
VIDEO_PATH = ""           # e.g. "/storage/emulated/0/DCIM/Camera/VID_20261124_221500.mp4"
                          # leave "" to use the newest video in VIDEO_DIR
VIDEO_DIR = "/storage/emulated/0/DCIM/Camera"
OUT_DIR = "/storage/emulated/0/Download"
CAPTURE_FPS = 240         # real shooting rate of the camera (frames per second)
TEST_TONE_HZ = 80         # tone you played (0 = unknown); used only for the report
ROI = None                # crop (x, y, width, height) in pixels, or None for whole frame
MAX_WIDTH = 320           # frames are shrunk to this width for speed
MAX_FRAMES = 4000         # stop after this many frames
ROLLING = False           # experimental rolling-shutter mode
ROW_FRACTION = 0.9        # fraction of the frame period spent reading rows (rolling mode)
REF_FRAMES = 31           # L4065: rolling mode uses the median of the first 31 frames as reference (waveform r 0.71 -> 0.99)
CALIBRATE_TONE_HZ = 0     # L4065: rolling mode only; a KNOWN tone in this video -> set ROW_FRACTION automatically (0 = off)
# ===========================================================================

TREMOR_BAND = (3.0, 15.0)  # Hz searched in tremor mode
import os, sys, wave, math
import numpy as np

try:
    import cv2
except ImportError:
    cv2 = None


def to_gray_small(frame):
    """BGR or gray frame -> float32 gray, cropped to ROI and shrunk to MAX_WIDTH."""
    if frame.ndim == 3:
        g = 0.114 * frame[..., 0] + 0.587 * frame[..., 1] + 0.299 * frame[..., 2]
    else:
        g = frame.astype(np.float32)
    if ROI:
        x, y, w, h = ROI
        g = g[y:y + h, x:x + w]
    g = g.astype(np.float32)
    if g.shape[1] > MAX_WIDTH:
        f = MAX_WIDTH / g.shape[1]
        if cv2 is not None:
            g = cv2.resize(g, (MAX_WIDTH, max(8, int(g.shape[0] * f))), interpolation=cv2.INTER_AREA)
        else:
            step = int(round(1 / f))
            g = g[::step, ::step]
    return g


def smooth(img):
    """Light blur so the motion estimate uses edges, not pixel noise."""
    if cv2 is not None:
        return cv2.GaussianBlur(img, (0, 0), 1.5)
    k = np.array([1, 4, 6, 4, 1], np.float32) / 16
    img = np.apply_along_axis(lambda r: np.convolve(r, k, 'same'), 1, img)
    return np.apply_along_axis(lambda c: np.convolve(c, k, 'same'), 0, img)


class MotionMeter:
    """Global sub-pixel shift of each frame relative to a reference frame
    (one-step Lucas-Kanade on the whole picture, using edges where the gradient
    is strong). Returns (dx, dy) in pixels, or per-row dx in rolling mode."""

    def __init__(self, ref):
        self.ref = smooth(ref)
        gy, gx = np.gradient(self.ref)
        self.gx, self.gy = gx, gy
        w = np.hypot(gx, gy)
        self.w = (w > np.percentile(w, 75)).astype(np.float32)   # keep strongest edges
        A = np.array([[np.sum(self.w * gx * gx), np.sum(self.w * gx * gy)],
                      [np.sum(self.w * gx * gy), np.sum(self.w * gy * gy)]], np.float64)
        self.Ainv = np.linalg.pinv(A)
        self.rowden = np.sum(self.w * gx * gx, axis=1) + 1e-6

    def shift(self, img):
        it = smooth(img) - self.ref
        b = -np.array([np.sum(self.w * self.gx * it), np.sum(self.w * self.gy * it)])
        return self.Ainv @ b

    def row_shift(self, img):
        it = smooth(img) - self.ref
        return -np.sum(self.w * self.gx * it, axis=1) / self.rowden


def frames_from_video(path):
    if cv2 is None:
        sys.exit("OpenCV is missing: in Pydroid open Pip and install opencv-python")
    cap = cv2.VideoCapture(path)
    if not cap.isOpened():
        sys.exit("Could not open video: " + path)
    header_fps = cap.get(cv2.CAP_PROP_FPS) or 0
    print("video:", path)
    print("file header says %.1f fps; using CAPTURE_FPS = %s" % (header_fps, CAPTURE_FPS))
    n = 0
    while n < MAX_FRAMES:
        ok, fr = cap.read()
        if not ok:
            break
        n += 1
        yield fr
    cap.release()


def synthetic_frames(freq, fps, seconds=2.0, amp_px=0.15, size=(120, 160), seed=1):
    """A textured scene shifted sideways by amp_px * sin(2 pi f t): a known answer."""
    rng = np.random.default_rng(seed)
    h, w = size
    base = rng.random((h, w * 2)).astype(np.float32)
    if cv2 is not None:
        base = cv2.GaussianBlur(base, (0, 0), 2.0)
    xs = np.arange(w * 2, dtype=np.float32)
    for k in range(int(seconds * fps)):
        d = amp_px * math.sin(2 * math.pi * freq * k / fps)
        img = np.empty((h, w), np.float32)
        for r in range(h):
            img[r] = np.interp(np.arange(w) + 40 - d, xs, base[r])
        img += rng.normal(0, 0.002, img.shape).astype(np.float32)   # sensor noise
        yield (img * 255).astype(np.uint8)


def rolling_rows(frames):
    """Per-row sideways shift of every frame (frames x rows), against a median reference when REF_FRAMES > 1."""
    buf, meter, rows = [], None, []
    for fr in frames:
        g = to_gray_small(fr)
        if meter is None:
            buf.append(g)
            if len(buf) < REF_FRAMES: continue
            meter = MotionMeter(np.median(np.array(buf), axis=0) if REF_FRAMES > 1 else buf[0])
            for b in buf: rows.append(meter.row_shift(b))
            continue
        rows.append(meter.row_shift(g))
    if meter is None:
        sys.exit("No frames read.")
    return np.array(rows)


def resample_rows(R, fps, row_fraction):
    """Place every row sample at its own capture time and resample evenly. Returns signal, rate, grid, in-readout mask."""
    n_f, n_r = R.shape
    t = (np.arange(n_f)[:, None] + row_fraction * np.arange(n_r)[None, :] / n_r) / fps
    order = np.argsort(t.ravel()); tt, vv = t.ravel()[order], R.ravel()[order]
    rate = min(fps * n_r, 48000)
    grid = np.arange(tt[0], tt[-1], 1.0 / rate)
    mask = ((grid * fps) % 1.0) <= row_fraction          # False inside the dead time between frames (interpolated, not measured)
    return np.interp(grid, tt, vv), rate, grid, mask


def calibrate_row_fraction(R, fps, tone_hz, lo=0.5, hi=1.0, step=0.005):
    """L4065: film a KNOWN tone once; the row timing that puts the most energy exactly at the tone is this phone's ROW_FRACTION."""
    best, best_rf = -1.0, None
    for rf in np.arange(lo, hi + 1e-9, step):
        x, rate, _, _ = resample_rows(R, fps, rf); x = clean(x, rate)
        n = len(x); X = np.abs(np.fft.rfft(x * np.hanning(n))); k = int(round(tone_hz * n / rate))
        v = X[max(k - 1, 0):k + 2].max()
        if v > best: best, best_rf = v, float(rf)
    return best_rf


LAST_ROLL = None


def measure(frames):
    global LAST_ROLL
    if ROLLING:
        global ROW_FRACTION
        R = rolling_rows(frames)
        if CALIBRATE_TONE_HZ:
            ROW_FRACTION = calibrate_row_fraction(R, CAPTURE_FPS, CALIBRATE_TONE_HZ)
            print("row timing calibrated from the %.1f Hz tone: ROW_FRACTION = %.3f (reuse this value for later videos from this phone)"
                  % (CALIBRATE_TONE_HZ, ROW_FRACTION))
        x, rate, grid, mask = resample_rows(R, CAPTURE_FPS, ROW_FRACTION)
        LAST_ROLL = (grid, mask)
        return x, rate
    meter, sig = None, []
    for fr in frames:
        g = to_gray_small(fr)
        if meter is None:
            meter = MotionMeter(g)
        sig.append(meter.shift(g))
    if meter is None:
        sys.exit("No frames read.")
    S = np.array(sig)                              # frames x 2 (dx, dy)
    S = S - S.mean(axis=0)
    u, s, vt = np.linalg.svd(S, full_matrices=False)   # main direction of vibration
    return S @ vt[0], CAPTURE_FPS


def clean(x, rate):
    x = x - np.mean(x)
    n = len(x)
    if n > 8:                                      # remove slow drift (< 2 Hz) in frequency domain
        X = np.fft.rfft(x)
        f = np.fft.rfftfreq(n, 1.0 / rate)
        X[f < 2.0] = 0
        x = np.fft.irfft(X, n)
    return x


def spectrum(x, rate):
    n = len(x)
    win = np.hanning(n)
    X = np.abs(np.fft.rfft(x * win))
    f = np.fft.rfftfreq(n, 1.0 / rate)
    return f, X


def write_wav(path, x, rate, out_rate=8000):
    """Write the motion signal as audio. Low sample rates are upsampled so phones can play it."""
    peak = np.max(np.abs(x)) or 1.0
    y = x / peak * 0.9
    if rate < out_rate:
        t_old = np.arange(len(y)) / rate
        t_new = np.arange(0, t_old[-1], 1.0 / out_rate)
        y = np.interp(t_new, t_old, y)
        rate = out_rate
    with wave.open(path, 'wb') as w:
        w.setnchannels(1); w.setsampwidth(2); w.setframerate(int(rate))
        w.writeframes((y * 32767).astype('<i2').tobytes())


def report(f, X, rate, tone, lines):
    nyq = rate / 2.0
    k = np.argsort(X[1:])[::-1][:5] + 1
    lines.append("Nyquist limit: %.1f Hz (nothing above this can be recovered)" % nyq)
    lines.append("Top 5 spectral peaks:")
    for i in k:
        lines.append("   %8.2f Hz   strength %.3g" % (f[i], X[i]))
    best = f[k[0]]
    if tone:
        if tone >= nyq:
            alias = abs(tone - round(tone / rate) * rate)
            lines.append("WARNING: the tone %.1f Hz is above the Nyquist limit; it can only appear "
                         "folded down near %.1f Hz (aliasing)." % (tone, alias))
            target = alias
        else:
            target = tone
        ok = abs(best - target) <= max(1.0, 2 * rate / len(X) / 2)
        lines.append("Expected peak: %.2f Hz  |  strongest peak: %.2f Hz  ->  %s"
                     % (target, best, "MATCH" if ok else "no match"))
        return ok
    return None


def plot(f, X, path, title):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return False
    plt.figure(figsize=(8, 4))
    plt.plot(f[1:], X[1:], lw=1)
    plt.xlabel("frequency (Hz)"); plt.ylabel("strength"); plt.title(title)
    plt.tight_layout(); plt.savefig(path, dpi=120); plt.close()
    return True


def tremor_peak(x, rate, band=None):
    """Dominant frequency inside the tremor band, refined by parabolic interpolation of the spectrum peak."""
    lo, hi = band or TREMOR_BAND
    n = len(x); X = np.abs(np.fft.rfft((x - x.mean()) * np.hanning(n))); f = np.fft.rfftfreq(n, 1.0 / rate)
    sel = np.where((f >= lo) & (f <= min(hi, rate / 2)))[0]
    if len(sel) == 0:
        return None
    k = sel[np.argmax(X[sel])]
    if 0 < k < len(X) - 1:
        a, b, c = X[k - 1], X[k], X[k + 1]; d = 0.5 * (a - c) / (a - 2 * b + c) if (a - 2 * b + c) != 0 else 0.0
        return float(f[k] + d * (f[1] - f[0]))
    return float(f[k])


def tremor_selftest_frames(freq, fps, seconds=10.0, amp_px=0.05, noise_dn=2.0, seed=7):
    rng = np.random.default_rng(seed); h, w = 120, 160
    base = rng.random((h, w * 2)).astype(np.float32)
    if cv2 is not None:
        base = cv2.GaussianBlur(base, (0, 0), 2.0)
    xs = np.arange(w * 2, dtype=np.float32)
    for k in range(int(seconds * fps)):
        d = amp_px * math.sin(2 * math.pi * freq * k / fps)
        img = np.empty((h, w), np.float32)
        for r in range(h):
            img[r] = np.interp(np.arange(w) + 40 - d, xs, base[r])
        yield np.clip(np.round(img * 255 + rng.normal(0, noise_dn, img.shape)), 0, 255).astype(np.uint8)


def newest_video(folder):
    vids = [os.path.join(folder, n) for n in os.listdir(folder)
            if n.lower().endswith((".mp4", ".mov", ".3gp", ".mkv", ".avi"))]
    if not vids:
        sys.exit("No videos found in " + folder)
    return max(vids, key=os.path.getmtime)


def run(mode=MODE):
    global CAPTURE_FPS
    if mode in ("tremor", "tremorselftest"):
        return run_tremor(mode)
    out_dir = OUT_DIR if os.path.isdir(OUT_DIR) else os.getcwd()
    lines = ["DCF Visual Microphone report", "mode: %s%s" % (mode, " (rolling)" if ROLLING else "")]
    if mode == "selftest":
        tone, CAPTURE_FPS = 37.0, 240
        frames = synthetic_frames(tone, CAPTURE_FPS)
        name = "visualmic_selftest"
        lines.append("synthetic scene: 0.15-pixel sideways vibration at %.1f Hz, filmed at %d fps" % (tone, CAPTURE_FPS))
    else:
        path = VIDEO_PATH or newest_video(VIDEO_DIR)
        frames = frames_from_video(path)
        tone = TEST_TONE_HZ
        name = "visualmic_" + os.path.splitext(os.path.basename(path))[0]
        lines.append("video: " + path)
    x, rate = measure(frames)
    x = clean(x, rate)
    lines.append("samples: %d at %.1f samples/s (%.2f s)" % (len(x), rate, len(x) / rate))
    lines.append("motion size (rms): %.4f pixels" % float(np.sqrt(np.mean(x ** 2))))
    f, X = spectrum(x, rate)
    ok = report(f, X, rate, tone, lines)
    wav = os.path.join(out_dir, name + ".wav"); write_wav(wav, x, rate)
    png = os.path.join(out_dir, name + "_spectrum.png")
    lines.append("sound file: " + wav)
    if plot(f, X, png, name):
        lines.append("spectrum picture: " + png)
    txt = os.path.join(out_dir, name + "_report.txt")
    lines.append("report: " + txt)
    with open(txt, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    return ok


def run_tremor(mode):
    """Neuron-driven motion a phone CAN see: tremor of a hand or finger resting on a card (motor-unit firing shows up as
    3-15 Hz movement). Film the card held still, phone on a stand. NOT FOR MEDICAL USE: this is a physics demo, not a test."""
    global CAPTURE_FPS
    lines = ["DCF Visual Microphone - tremor mode (not for medical use)"]
    if mode == "tremorselftest":
        ok_all = True
        for fps in (30, 60, 240):
            CAPTURE_FPS = fps
            x, rate = measure(tremor_selftest_frames(9.3, fps))
            pk = tremor_peak(clean(x, rate), rate)
            ok = pk is not None and abs(pk - 9.3) <= 0.2; ok_all &= ok
            lines.append("self-test %3d fps: 0.05-pixel 9.30 Hz motion -> peak %.2f Hz  %s" % (fps, pk or -1, "MATCH" if ok else "no match"))
        print("\n".join(lines)); return ok_all
    path = VIDEO_PATH or newest_video(VIDEO_DIR)
    x, rate = measure(frames_from_video(path)); x = clean(x, rate)
    pk = tremor_peak(x, rate)
    lines.append("video: " + path)
    lines.append("dominant motion frequency in %.0f-%.0f Hz: %s" % (TREMOR_BAND[0], TREMOR_BAND[1], "%.2f Hz" % pk if pk else "none (band above Nyquist?)"))
    lines.append("motion size (rms): %.4f pixels" % float(np.sqrt(np.mean(x ** 2))))
    print("\n".join(lines)); return pk


if __name__ == "__main__":
    run()
