# DCF L4097: a sideband criterion for row-timing calibration (0 of 3) + the 37-family literature check
### Gary NOMEREMORTAL Mason · Digital Cosmos Framework · 28 November 2026

**Protocol.** Pre-registered and pushed before computing (commit 50f120f); output committed before grading; separate grader.

**Record:** 664 + 0 = **664 of 930 passed, plus 3 split.**

## L4097 (0 of 3)
| true row timing | energy criterion error | sideband criterion error |
|---|---|---|
| 0.65 | +0.006 | −0.006 |
| 0.70 | −0.003 | −0.005 |
| 0.85 | −0.002 | −0.008 |
| 0.95 | +0.017 | −0.001 |

**Graded:**
- **S1 FAIL:** worst sideband error 0.008 (needed ≤ 0.005).
- **S2 FAIL:** mean −0.005 (needed ≤ 0.003).
- **S3 FAIL:** sideband better on only 2 of 4.

**Reading.**
- **The sideband criterion has a tighter spread** (worst 0.008 vs 0.017) but leans slightly low.
- **The energy criterion is unbiased on average here** (mean +0.005) but occasionally misses by 0.017.
- **Neither reaches ±0.005.** With one tone in 2 s of video, about ±0.01–0.02 looks like the practical floor; L4065 showed that is enough to recover unknown tones.
- **No default change.** Per the frozen interpretation, the energy criterion stays.

## The 37-family theorem against the literature
- **Already known:** OEIS A074992 already lists 1, 37, 3367, 333667, … and notes (Fiorini, 2016) that the reciprocals have period 3n.
- **Not on that page:**
  - the centred-hexagonal identity, a(n) = A003215((10ⁿ−1)/3);
  - a proof;
  - the general base b ≡ 1 (mod 3) (verified in 72 cases).
- **Written up:** `notes/Note_37_family_centred_hexagonal_v0.1.md` in dcf-math, with a ready OEIS comment for Gary to submit under his own account.

## Summary
- **Neither calibration criterion beats about ±0.01;** a single short tone sets the floor. The default is unchanged.
- **The 37-family period is already in OEIS,** but its centred-hexagonal identity, the proof and the base generalisation are not; an OEIS comment is drafted.

## Potential leads for future research
- **[Two-tone calibration]:** two known tones, or a 10-second clip; pre-register ±0.005.
- **[OEIS comment]:** Gary submits the A074992 comment.
- **[Other cyclotomic families]:** which Φ_m(bᵏ)/g are figurate?
