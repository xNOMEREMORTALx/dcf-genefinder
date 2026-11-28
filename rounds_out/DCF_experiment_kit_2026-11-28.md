# DCF experiment kit: three phone experiments for today
### Gary NOMEREMORTAL Mason · Digital Cosmos Framework · 28 November 2026

**Before you start:**
- Pydroid 3 with numpy and opencv-python installed.
- Phone bundle **6.23** (menu 38), plus `dcf_corner.py` in the same folder.
- A tripod or stand. **None of this is for medical use.**

## 1. Visual microphone above the frame rate (rolling shutter)
1. **Set up:** a crisp packet, foil or a paper card about 10–20 cm from a speaker, in bright light. Phone on a stand, about 30 cm away, filling the frame with the object.
2. **Play the tones:** `dcf_two_tone_530_1530Hz_60s.wav` at moderate volume. If the speaker sounds weak, use `dcf_two_tone_400_900Hz_60s.wav`.
3. **Film:** 10–20 s in slow motion at **240 fps**. Don't touch the phone.
4. **Analyse:** menu 38 → `video` → the video path → frame rate 240 → tone 1530 (or 900) → rolling-shutter mode **y** → calibration tone **1530** (or 900).
5. **Report back:** the calibrated row timing, and whether the 530 Hz (or 400 Hz) peak shows. Or send me the video.

## 2. Tremor (neuron-driven motion)
1. **Set up:** rest a hand lightly on a stiff card, held out in the air or resting on a table edge.
2. **Film:** the card edge for 10 s at 240 fps (or 60 fps), phone on a stand.
3. **Analyse:** menu 38 → `tremor` → the video path → the shooting rate.
4. **Report back:** the dominant frequency. It is a physics demo, **not a health test**.

## 3. Corner camera (one photo)
1. **Set up:** a white or pale wall in a dim room. Two lamps (or phone torches) on the hidden side of a **vertical edge** (door frame or chair leg) that stands between them and the wall.
2. **Measure with a tape, in metres, from the left end of the wall strip you will photograph** (to about 1 cm; L4109 shows 2 cm errors hurt):
   - the length of the wall strip in the photo;
   - the edge's distance from the wall, and its left and right positions along the wall;
   - the lamps' distance from the wall, and the left–right span they could be in.
3. **Photograph:** the wall strip, straight on, phone level, no flash.
4. **Analyse:**

       python3 -I dcf_corner.py photo.jpg --wall 0 WALL_LEN --occ EDGE_LEFT EDGE_RIGHT EDGE_DEPTH --hidden SPAN_LEFT SPAN_RIGHT LAMP_DEPTH

5. **Report back:** the "brightest positions" against where the lamps really were. Or send me the photo and the measurements.

## Summary
- **Three experiments, about 30 minutes in total.** Each gives a pre-registerable real-world result: rolling-shutter calibration, a tremor frequency, and a corner-camera reconstruction.

## Potential leads for future research
- **[Real-phone rounds]:** I will pre-register predictions for each before you send data. Tell me which you are doing.
