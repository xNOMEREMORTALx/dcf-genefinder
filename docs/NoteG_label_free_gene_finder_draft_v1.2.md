# Note G: a label-free, phone-sized gene finder from the period-3 rhythm, with genetic-code detection and honest blind testing
### Draft v1.2 · Gary NOMEREMORTAL Mason · Digital Cosmos Framework · 9 November 2026

**Status.**
- Tool: `dcf_genefinder_mobile.py` **v2.32** (SHA-256 fa30258b…297c), numpy only; it runs in Pydroid 3 on Android.
- Code and grader: public repository `xNOMEREMORTALx/dcf-genefinder` (MIT).
- **New in v1.2:** v2.32, start-codon re-pick fix, blind-tested twice (§3k); code record 134 of 137.
- **New in v1.1:** v2.31 (genome-size prior for the code decision); the code record through L3185 (**114 of 117**); blind sets 9–15 (§3i); the metal profile and selenocysteine finder with their limits (§5c); two new corrections (§6).
- **New in v0.8:** blind sets 4 and 5 plus *S. epidermidis*; start accuracy as a third measure; v2.28 (start regression fixed), v2.29 (learned start signal), v2.30 (30-nt window).
- **New in v0.7:**
  - four blind genome sets (L2953, L2957, L2965, L2977/L2981);
  - the GC-rich fix (v2.26) and its extension to all genomes (v2.27);
  - one regression and two labelling errors, all corrected;
  - the short-gene limit.
- Every figure is from the cited round. *Development* means the result was tuned on the same genomes. *Blind* means pre-registered and frozen before the file was opened.

---

## 1. Idea
**The signal.** Protein-coding DNA carries a period-3 rhythm: the three codon positions use the four bases differently. The three "phase sums" of each base are 5–30 times more unequal in genes than in shuffled DNA, in all 31 genomes tested (L2866).

**The tool.**
1. **Self-training.** It trains a period-3 model on the genome itself, with no reference genes.
2. **Calling.** It scores every open reading frame **from its best start codon** (v2.26/v2.27) and keeps the frames that fit.
3. **Clean-up.** It removes antisense shadows and duplicate frames.
4. **Starts.** It anchors starts with the genome's own 16S rRNA anti-Shine–Dalgarno sequence, then **learns the genome's own start signal** (bases −30..−1 and start-codon preferences) from its confident genes and re-picks every start (v2.29–v2.30).
5. **Annotation is used only for grading,** never for training.

## 2. Genetic-code detection (L2720–L2831, L2957)
| test | what it decides | principle |
|---|---|---|
| Read-through | Is TGA a stop? | Does coding signal continue past TGA? |
| **Length ratio** (L2827) | TGA = Trp? | Median gene length with TGA as sense ÷ as stop: standard 0.98–1.38, TGA = Trp 1.64–1.98 |
| Disagreement rule | read-through vs ratio | Outside the 1.4–1.6 band, the ratio decides |
| **Tiny AT-rich rule** (v2.25) | genomes < 300 kb, GC < 0.25 | The ratio decides at a 1.45 cut, with no uncertain band (fixes *Nasuia*) |
| AGR fragmentation (L2728) | Are AGA/AGG stops (mitochondria)? | Top-gene length ratio: stops 0.86–0.97, sense 0.31–0.40 |
| Reliability self-check (L2818) | Is the model failing? | Median predicted gene < 400 nt raises a warning |

| **Genome-size prior** (v2.31, L3149) | read-through says 11, ratio says 4 | Keep 11 if the genome is > 2.5 Mb (known table-4 bacteria are ≤ about 2.2 Mb) |

**Blind record: 134 correct of 137 genomes** (L2953–L3229, after the L2973 correction).
- **The third error:** *Gluconobacter oxydans* (L3145). The length ratio (1.69) inflated in a medium-GC proteobacterium with about 70% TGA stops and overrode a correct read-through call. v2.31 fixes this for genomes > 2.5 Mb; blind test L3153, 4 of 4 passed. Genomes < 2.5 Mb with ratio > 1.6 remain unprotected.
- **The first two errors:** *Nasuia* (L2953, flagged UNCERTAIN; now fixed by v2.25, post hoc) and *Hodgkinia* (an early round; fixed in v2.20).
- **The read-through test is unreliable in AT-rich genomes.** In *Sulcia*, *Nasuia* and *Buchnera* it was wrong each time it disagreed with the length ratio, and the ratio was right.

**Known limits.**
- Sense-codon reassignments (AGR = Ser, AAA = Asn) are invisible to a rhythm method.
- The tiny-genome branch of v2.25 has not yet fired in a blind genome.

## 3. Blind accuracy (all pre-registered; one grader for both tools)
**Grading:**
- **found** = annotated gene covered ≥ 50% by a same-strand prediction;
- **correct** = prediction covered ≥ 50% by an annotated gene.

### 3a. Set 1 (L2953, v2.24)
| genome | GC | DCF found / correct | Prodigal found / correct |
|---|---|---|---|
| *Micrococcus luteus* | 0.73 | 74.3 / 91.1 | 98.5 / 99.2 |
| *Corynebacterium glutamicum* | 0.54 | 85.8 / 97.6 | 98.3 / 96.3 |
| *Sulcia muelleri* GWSS | 0.22 | 88.3 / 96.6 | 98.2 / 96.5 |
| *Nasuia* (code miscalled) | 0.17 | 83.9 / 94.7 | 77.4 / 92.9\* |
| *Blochmannia floridanus* | 0.27 | 91.9 / 99.8 | 99.5 / 98.2 |
| *Prochlorococcus* MED4 | 0.31 | 77.0 / 99.5 | 97.9 / 95.8 |
| *Pelagibacter ubique* | 0.30 | 92.3 / 99.4 | 99.6 / 99.2 |

\* Prodigal was run with the DCF tool's wrong code for *Nasuia*.

### 3b. *Arabidopsis* mitochondrion (L2957, v2.24)
- 22 of 24 single-piece genes found; 10.9% of the 202 calls correct. Prodigal: 100% / 9.6%.
- The genome is about 90% non-coding. The tool's own reliability warning fired.

### 3c. Set 2 (L2965, v2.25)
| genome | GC | DCF found / correct | Prodigal found / correct |
|---|---|---|---|
| *Carsonella* HC | 0.14 | 96.4 / 98.4 | 98.4 / 99.0 |
| *Carsonella* CE | 0.14 | 96.8 / 98.9 | 97.9 / 99.5 |
| *Sulcia* SMDSEM | 0.23 | 88.2 / 93.0 | 99.2 / 93.4 |
| *Uzinura diaspidicola* | 0.30 | 91.2 / 98.6 | 99.6 / 89.7 |
| *Buchnera aphidicola* Cc | 0.20 | 91.0 / 99.7 | 99.2 / 99.2 |
| *Streptomyces coelicolor* | 0.72 | 74.5 / 91.8 | 98.1 / 98.8 |
| *M. genitalium*† | 0.32 | 92.2 / 98.2 | 98.9 / 97.9 |
| *Deinococcus* chr 1† | 0.67 | 82.9 / 94.6 | 98.2 / 98.9 |

† Not truly unseen: both had been used in earlier development (corrected in L2973). They are reported but excluded from the blind counts.

### 3d. The GC-rich fix, held-out test (L2969, v2.26)
**Cause.** In GC-rich DNA, stop codons are rare, so reading frames run far upstream of the real gene. The non-coding part dragged whole-frame scores below zero; this was the reason for 362 of 543 *Micrococcus* misses.

**Fix.** Score each frame from its best start codon, and lower the minimum to 180 nt.

| genome (never used for tuning) | v2.25 | **v2.26** | Prodigal |
|---|---|---|---|
| *Streptomyces coelicolor* | 74.5 / 91.8 | **94.9 / 98.6** | 98.1 / 98.8 |
| *Deinococcus* chr 1 | 82.9 / 94.6 | **93.3 / 97.3** | 98.2 / 98.9 |

### 3e. Set 3 (L2977/L2981, v2.27 vs v2.26)
| genome | GC | v2.26 | **v2.27** | Prodigal |
|---|---|---|---|---|
| *Vibrio cholerae* chr I | 0.48 | 88.2 / 99.4 | **89.7 / 98.8** | 98.2 / 98.6 |
| *Listeria monocytogenes* | 0.38 | 91.3 / 99.7 | **93.5 / 99.1** | 99.3 / 99.0 |
| *Lactococcus lactis* | 0.35 | 86.0 / 98.6 | **89.3 / 98.0** | 98.3 / 96.8 |
| *Clostridium acetobutylicum* | 0.31 | 89.3 / 99.4 | **92.7 / 98.8** | 98.6 / 98.4 |
| *Methanosarcina acetivorans* | 0.43 | 75.4 / 96.8 | **79.5 / 93.9** | 93.3 / 93.8 |
| *Synechocystis* PCC 6803 | 0.47 | 80.6 / 98.6 | **84.3 / 97.0** | 97.2 / 97.9 |
| *Pseudomonas aeruginosa* | 0.66 | 96.6 / 97.6 | **96.6 / 97.6** | 99.4 / 97.7 |
| *Caulobacter crescentus* | 0.67 | 89.9 / 98.0 | **89.9 / 98.0** | 92.9 / 98.9 |

### 3f. Set 4 (L2997): the start fix, v2.28 vs v2.27
**Measure.** Exact start = the predicted start equals the annotated start, among genes whose stop is predicted on the right strand.

| genome | found | correct | exact start v2.27 → **v2.28** | Prodigal start |
|---|---|---|---|---|
| *Chlamydia trachomatis* | 95.0 | 96.8 | 59.5 → **74.5** | 90.3 |
| *Halobacterium salinarum* | 90.1 | 98.6 | 70.8 → **76.8** | 93.4 |
| *Thermococcus kodakarensis* | 89.5 | 97.9 | 75.1 → **83.0** | 85.9 |
| *Streptococcus pneumoniae* | 81.4 | 97.8 | 57.6 → **82.3** | 91.3 |
| *M. tuberculosis* | 94.6 | 90.9 | 50.7 → **64.1** | 77.3 |
| *Rhodopseudomonas palustris* | 94.3 | 97.7 | 62.8 → **75.2** | 89.5 |
| *Bacteroides thetaiotaomicron* | 89.1 | 97.4 | 39.2 → **52.4** | 88.7 |

### 3g. Set 5 (L3005) and *S. epidermidis* (L3013): learned starts, v2.29/v2.30 vs v2.28
| genome | found | correct | exact start v2.28 → **v2.29** | Prodigal start |
|---|---|---|---|---|
| *Porphyromonas gingivalis* (no SD) | 82.4 | 95.0 | 58.1 → **69.3** | 86.0 |
| *Cytophaga hutchinsonii* (no SD) | 93.1 | 97.8 | 61.5 → **71.1** | 89.3 |
| *Sulfolobus solfataricus* (leaderless) | 89.2 | 97.8 | 70.3 → **72.3** | 85.9 |
| *Francisella tularensis* | 90.4 | 93.4 | 75.2 → **76.9** | 86.6 |
| *Corynebacterium diphtheriae* | 88.6 | 94.9 | 73.2 → **74.9** | 84.5 |
| *Lactobacillus plantarum* | 91.5 | 97.8 | 81.8 → **87.8** | 91.7 |
| *C. difficile* | 91.1 | 99.2 | 89.3 → **92.0** | 93.5 |
| *S. epidermidis* RP62A (v2.30) | 89.1 | 99.1 | 87.4 → 91.2 → **92.0** | 93.1 |

**Two pre-registered start predictions failed in set 5.**
- Gains of 3 or more points came in only 3 of 7 genomes.
- The mean gain in the non-SD group was 7.6, against a target of 8.
- The learned signal helps most where the fixed 16S anchor was wrong for the genome. It helps little where there is no upstream signal at all (leaderless genes).

### 3h. Set 6 (L3053/L3057, v2.30) and the promoter self-check
| genome | found | correct | exact start | promoter E17 |
|---|---|---|---|---|
| *S. pyogenes* | 88.8 | 98.9 | 87.0 | **3.78** clear |
| *E. faecalis* V583 | 88.4 | 95.9 | 78.1 | **2.53** clear |
| *C. perfringens* | 94.5 | 99.5 | 94.1 | **2.17** clear |
| *B. cereus* | 87.8 | 98.3 | 86.8 | **1.59** clear |
| *Xanthomonas* | 90.1 | 95.2 | 74.5 | 1.22 |
| *Legionella* | 93.3 | 97.6 | 80.0 | 1.13 |
| *Acinetobacter* | 93.2 | 97.9 | 82.9 | 1.13 |
| *Shewanella* | 91.3 | 96.2 | 81.6 | 0.89 |

**Promoter self-check (`dcf_promoters.py`).**
- **What it does.** It learns −35/−10 boxes from the 250 nt upstream of predicted operon-first genes, and reports a clear signal when the 17-nt spacer share is above 1.5× its shuffled control.
- **Record: 22 of 22 genomes sorted correctly** (14 development, 8 blind). Low-GC Firmicutes are clear (1.59–3.78); proteobacteria, actinobacteria, Bacteroidetes and archaea are not (0.76–1.22).
- **In clear genomes,** strong sites sit before likely operon starts at 6–13× the shuffled rate.

**Overall blind picture.**
- **Precision** is 91–99.8%, equal to or above Prodigal's in about half of the genomes.
- **Sensitivity** is now 89–97% in GC-rich genomes (after v2.26) and 79–96% elsewhere, against Prodigal's 93–99.6%.
- **Exact starts** are 52–92% (v2.28–v2.30), against Prodigal's 77–93.5%. The gap ranges from about 1 point (*S. epidermidis*) to 17 points (*Porphyromonas*).

### 3i. Blind sets 9–15 (L3137–L3185, v2.30/v2.31)
These rounds tested mainly the metal and selenocysteine add-ons; the gene finder ran blind on every genome.
| round | genomes | code correct | genes found | correct | exact starts |
|---|---|---|---|---|---|
| L3153 (set 11) | 6 | 6 | 89.2–94.6% | 91.0–97.7% | 72.1–91.5% |
| L3165 (set 12) | 10 | 10 | 86–94% | 94–99% | 69–94% |
| L3181/L3185 (sets 14–15) | 10 | 10 | 86.0–95.3% | 95.5–99.2% | 79.0–91.0% |

(Sets 9, 10 and 13: code 6/6, 8/9 and 8/8; full grades are in their reports.)

### 3j. Head-to-head with Prodigal on sets 14–15 (L3205)
Same 10 genomes, same grader (`dcf_grade.py`), Prodigal via pyrodigal 3.7.1 in single-genome mode with table 11.

| measure | DCF v2.31 mean (range) | Prodigal mean (range) | DCF ≥ Prodigal |
|---|---|---|---|
| genes found | 91.8% (86.0–95.3) | 97.9% (96.3–98.9) | 0 of 10 |
| predictions correct | **98.3%** (95.5–99.2) | 98.6% (96.3–99.3) | 3 of 10 |
| exact starts | 84.8% (79.0–91.0) | 91.9% (88.6–94.4) | 0 of 10 |

**Reading.**
- **Precision is level with Prodigal** (0.3 points apart on average).
- **The gaps are sensitivity (about 6 points) and starts (about 7 points).** These are the limits already listed in §7.
- **Caveat:** the reference annotations come from NCBI's pipeline, which uses its own gene models, so neither tool is graded against experimental truth.

### 3k. Gene starts, v2.32 (L3209–L3229)
**Diagnosis (development, blind_14/15).** Start errors in GC-rich genomes were mostly too short (2:1). The commonest case was a downstream TTG/GTG chosen over the annotated ATG, because the learned codon weights were circular, i.e. learned from the tool's own picks.

**Fix.**
- Fixed penalties GTG −4 and TTG −5 in the learned re-pick, and coding weight 0.75.
- From v2.32b: a re-picked start must leave the minimum gene length. Without this guard, a few genes shrank to 9–120 nt fragments.

| blind test | genomes | exact starts v2.31 → v2.32 | genes found | verdict |
|---|---|---|---|---|
| L3213 (v2.32 candidate, no guard) | 13 | 82.3 → 84.5 (+2.2) | −0.6 in 1 genome | T1, T2 pass; **T3 fail** (not released) |
| **L3229 (v2.32b)** | 6 | +2.3 (largest +5.5, +5.3 in GC-rich) | +0.04 to +0.48 | **D4–D6 pass → released** |

Prodigal on the L3213 set: exact starts 88.4%. The gap to Prodigal narrows from about 6 points to about 4.

## 4. Development accuracy and drivers (unchanged from v0.6)
- Precision falls with GC (ρ = −0.63 across 29 genomes, L2871). At a given GC, a stronger codon rhythm raises it.
- Animal mitochondria: 85–100% found. Reduced symbionts: 83–98%.

## 5. Extra outputs (all label-free)
| output | validation |
|---|---|
| Per-gene fit and class (possible foreign / strong codon bias) | phage genes ×5.0 in the lowest 5% (*B. subtilis*), translation genes ×8.6 in the top 5% (*S. aureus*) (L2731) |
| Genomic-island map | *B. subtilis* SPβ 92% phage-type genes (L2735) |
| Frameshift check | 100% of synthetic 1-nt indels caught, 2% false alarms (L2731) |
| Stem-loop scan (< 100 kb) | human mtDNA 5 vs 1.6 expected, including O_L (L2678) |
| Codon-phase strength (fixed in v2.24) | — |
| Helical rhythm | median 11.05 bp; 26/31 above null (L2913). It is **not** explained by GC, size or lifestyle (L2953, symbiont hypothesis refuted) |
| Dispersed-repeat flag | *B. pertussis* IS481 157/160 caught (L2921) |

## 5b. Read-through detection: pyrrolysine and selenocysteine (`dcf_pylflag.py`, L3069–L3089)
**Method.** Two consecutive same-strand predictions joined **in frame across a single stop codon** (each part ≥ 150 nt, no other in-frame stop between them) are listed as read-through candidates:
- **TAG:** pyrrolysine?
- **TGA:** selenocysteine?

**Gene-level results:**
- **Pyrrolysine.** In 5 pyrrolysine users, all **27 of 27** methylamine methyltransferases carrying an in-frame TAG are caught (*M. acetivorans* 7, *M. barkeri* 6, *M. mazei* 7, *M. burtonii* 6, *D. hafniense* 1).
  - The TAG sits at the same codon in each family: **MtmB codon 202, MttB codon 334, MtbB codon 356**. I recall these as the known pyrrolysine positions but did not check a source this round.
  - The automated RefSeq annotation marks **all 27 as pseudogenes** ("internal stop").
- **Selenocysteine.** *E. coli* *fdhF*, *fdnG* and *fdoG* are caught as TGA joins matching the annotated coordinates exactly (3 of 3).

**The genome-level score is only a hint.** The TAG excess over the genome's own TAA/TGA joins:
- flags pyrrolysine methanogens (*M. mazei* 5.06, *M. burtonii* 3.72, *M. barkeri* 2.06, *M. acetivorans* 2.90);
- misses *D. hafniense* (1.28; many ordinary TAG stops);
- non-users sit around 1.35, with *S. pneumoniae* at 2.04 (TAG pseudogenes).

Blind set 7: 3 of 4 predictions passed.

**Blind set 8 (L3093, pre-registered, 3 of 3 passed).** In *Methanohalophilus*, *Methanolobus*, *Methanosalsum* and *Methanomethylovorans*, **21 of 21** pyrrolysine methyltransferases are listed (9–20 TAG candidates per genome). Their TAG codons are again conserved (MtbB 355/356, MttB 334/324, MtmB 202 in 7 of 10). The genome-level excess would have flagged only 1 of the 4.

## 5c. Metal profile and selenocysteine candidates (`dcf_metals.py`, `dcf_secflag.py`)
- **Fe–S profile:** C-x-x-C motifs per 1,000 amino acids of the predicted proteins, with a frozen cut of 0.50 (L3141).
  - It separates Fe–S-rich anaerobes (methanogens, acetogens, sulfate reducers) and aerobic H2 oxidisers from aerobes and lactic fermenters.
  - It is **not** an oxygen meter (L3173).
  - The anaerobe/aerobe contrast itself is published (Andreini *et al.*, PLoS ONE 2017); see Note H in dcf-math.
  - R (observed/expected CxxC) is unstable below 0.4% cysteine, and the tool now warns.
- **Selenocysteine finder:** recall 65% overall, 57% blind in AT-rich genomes, but too many candidates (> 100 per genome).
  - An optional coding-potential filter (`--coding 0.05`, L3197) kept 28/28 caught sites with AUC 0.87 in development.
  - Its blind test is pre-registered.

## 6. Corrections log
1. **v2.22–v2.23:** a strand-coordinate bug in the codon-phase strength line. Fixed in v2.24; gene calls were unaffected.
2. **L2965:** *M. genitalium* and *Deinococcus* were wrongly labelled unseen. Their blind status was withdrawn (L2973), and a list of every genome used now gates new requests.
3. **v2.26:** best-start scoring also altered code detection, and *Hodgkinia* (code 4) was called standard. The L2969 side-effect check covered only genomes with GC ≤ 0.55. In **v2.27**, code detection is restored to v2.25, and a 12-genome code check (codes 2, 4, 11) now runs before release.
5. **Reader (L3069):** `dcf_grade.py` and `dcf_promoters.py` dropped IUPAC ambiguity letters, shifting coordinates in some genomes. Grades were unaffected (coordinates only); the promoter verdicts were re-checked and unchanged. Fixed.
6. **Selenocysteine counts (L3197):** counts in L3181/L3185 used `grep 'aa:Sec'`, which also matches tRNA-Sec anticodon lines. They are recounted from `transl_except`; no verdict changed.
7. **Seen list (L3177):** *L. plantarum* (set 5) was missing from `seen_accessions.txt`; a name list (`seen_names.txt`) now backs it up.
4. **v2.26–v2.27:** start selection searched only downstream of the best-start trim, and exact starts fell (*E. coli* 81.4% → 63.0%). The 50%-overlap grader could not see this. **v2.28** fixed it, and exact start is now reported by `dcf_grade.py` for every release.

## 7. Honest limits
- **Short genes are the remaining gap.** In *E. coli* and *Methanosarcina*, about half of the missed genes are under 300 nt, and most sit where nothing is called on either strand.
- **A rhythm-only rescue pass was tried and rejected** (L2981). The loosest setting added 243 *E. coli* calls but found only about 43 real genes, at −3.6 points of precision.
- Closing the last 5–14 points to Prodigal needs non-rhythm evidence (for example, conservation).
- **Archaea:** *Methanosarcina* has the lowest sensitivity (79.5%). The leaderless genes of *Sulfolobus* give the tool little upstream signal to learn starts from.
- **Plant mitochondria** give low precision (about 11%), as does Prodigal.
- **Reference annotations** are themselves produced by gene-finding pipelines, which may favour Prodigal-like calls.

## Summary
- **What it is.** A gene finder that teaches itself each genome's coding rhythm, needs only numpy, and runs on a phone.
- **Genetic code.** It works out each genome's genetic code by itself, correct in 114 of 117 blind genomes.
- **Promoters.** A companion script tells whether a genome's promoters can be read from sequence alone, correct in 22 of 22 genomes, and learns them where they can.
- **Rare amino acids.** A second script lists genes that read through a stop codon. It caught all 48 pyrrolysine genes in nine microbes (21 of them in a blind test), which the official annotation had wrongly marked as broken, and *E. coli*'s three selenocysteine genes.
- **Blind accuracy.** On 57 genomes never seen during development, its calls are right 91–99.8% of the time, about as often as the standard tool, Prodigal.
- **Gene starts.** It now learns where each genome's genes start and gets 52–92% of starts exactly right, against Prodigal's 77–93.5%.
- **Sensitivity.** It finds fewer genes than Prodigal, mainly the very short ones. Its biggest earlier weakness, GC-rich genomes, was traced to one cause and fixed: held-out *Streptomyces* rose from 74.5% to 94.9% of genes found.
- **Honesty.** Every error found along the way is logged and corrected in §6.

## Potential leads for future research
- **[G10] Short genes with conservation.** A small built-in set of conserved short proteins (for example, ribosomal proteins) as non-rhythm evidence.
- **[G9] Archaeal starts.** Test whether leaderless mRNAs and weak Shine–Dalgarno signals explain the *Methanosarcina* gap. Develop on *Pyrococcus* and *M. jannaschii*.
- **[G6b] A tiny TGA = Trp symbiont,** to test the v2.25 branch blind.
- **[G2] tRNA-anticodon module** for sense-codon reassignments.
- **[Submission]** With four blind sets this is close to a short software note (for example, *Bioinformatics Advances* or the *Journal of Open Source Software*). Needed first: a public code repository and a frozen version number.
