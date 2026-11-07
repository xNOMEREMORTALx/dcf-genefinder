# DCF gene finder

A label-free, phone-sized prokaryotic gene finder built on the period-3 rhythm of protein-coding DNA. It needs only Python 3 and numpy, and runs on Android (Pydroid 3).

**Author:** Gary NOMEREMORTAL Mason (Digital Cosmos Framework).
**Version:** 2.30.

## What it does
- **Self-training.** Learns a period-3 coding model from the genome itself. No reference genes or training data are used.
- **Genetic code.** Detects the genetic code automatically (standard, TGA = Trp "code 4", and mitochondrial variants). Blind record: 66 correct of 68 unseen genomes.
- **Gene calls.** Calls genes, scoring each reading frame from its best start codon.
- **Starts.** Learns each genome's own start signal: the upstream motif and start-codon preferences.
- **Extra outputs:**
  - per-gene fit / possible-foreign tag;
  - genomic-island map;
  - frameshift candidates;
  - repeat flag;
  - helical-rhythm line;
  - stem-loops for small genomes.

## Usage
```
python3 dcf_genefinder_mobile.py genome.gb --table auto --out pred.gff3
python3 dcf_grade.py genome.gb pred.gff3      # found / correct / exact start against the file's own annotation
python3 dcf_promoters.py genome.gb pred.gff3  # learned -35/-10 promoter boxes with a self-check (clear in low-GC Firmicutes)
python3 dcf_pylflag.py genome.gb pred.gff3    # TAG read-through (pyrrolysine) flag: catches all 7 known Methanosarcina pyrrolysine genes
```
Promoter self-check, blind set 6 (pre-registered): Firmicutes 1.59–3.78 (3/3 clear), proteobacteria 0.89–1.22 (0/4 clear); exact starts 74.5–94.1%.
```
```
Input is a GenBank file with sequence (NCBI efetch `rettype=gbwithparts`). The annotation is used only for grading, never for training.

## Blind accuracy (pre-registered; predictions frozen before the files were opened)
| set | genomes | genes found | predictions correct | Prodigal found / correct |
|---|---|---|---|---|
| 1 (v2.24) | 7 | 74–92% | 91–99.8% | 98–99.6 / 93–99 |
| 2 (v2.25) | 6 new | 88–97% | 93–99.7% | 97.9–99.6 / 89.7–99.5 |
| GC-rich held-out (v2.26) | *Streptomyces*, *Deinococcus* | 94.9%, 93.3% | 98.6%, 97.3% | 98.1–98.2 / 98.8–98.9 |
| 3 (v2.27) | 8 | 79.5–96.6% | 93.9–99.1% | 92.9–99.4 / 93.8–99.0 |
| 4 (v2.28) | 7 | 81.4–95.0% | 90.9–98.6% | 93.0–99.3 / 93.1–98.9 |
| 5 (v2.29) | 7 | 82.4–93.1% | 93.4–99.2% | 91.9–99.4 / 92.5–99.1 |
| 6 (v2.30) | 7 | 87.8–94.5% | 95.2–99.5% | not run |

**Grading.**
- **found** = annotated CDS covered ≥ 50% by a same-strand prediction.
- **correct** = prediction covered ≥ 50% by an annotated CDS.

**Exact starts** (blind): set 4 (v2.28) 52–83%; set 5 (v2.29) 69–92%, higher than v2.28 in all 7 genomes (+1.7 to +11.2 points) but short of two pre-registered targets (3 of 5 predictions passed; see `PREREGISTERED_blind5.txt`). Prodigal: 84.5–93.5% on set 5. *S. epidermidis* RP62A (blind, v2.30): found 89.1%, correct 99.1%, exact start 92.0%.

## Honest limits
- **Sensitivity.** The tool finds fewer genes than Prodigal (gap 3–14 points). The missed genes are mostly short (< 300 nt) with weak coding rhythm, and a rhythm-only rescue was tested and rejected.
- **Starts.** Exact starts trail Prodigal, most where Shine–Dalgarno sites are rare.
- **Codon reassignments.** Sense-codon reassignments (AGR = Ser, AAA = Asn) are not detected.
- **Plant mitochondria.** Precision is low (≈ 11%, as for Prodigal).

## Corrections
Every error found so far is logged in `CHANGELOG.md`, including:
- a code-call regression in v2.26;
- a start-codon regression in v2.26–v2.27;
- two genomes wrongly labelled unseen.

## Citation
Mason, G. N. *DCF gene finder: a label-free period-3 gene finder with genetic-code detection*, v2.29 (2026).

## License
MIT. See `LICENSE`.
