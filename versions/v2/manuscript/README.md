# Modular recurrences and congruences for 14325 packings

The note proves both observations in Section 7, item 11, page 49 of
Blitvić, Kammoun and Steingrímsson, *Permutations, Pattern Covers, Posets and
Positivity: Part I: Patterns and Posets*, arXiv:2610.29846v1 (2026).
For all positive integers j the packing count is one modulo five; for j at
least two it is three times the j-th Catalan number modulo four. These
congruences imply both the modulo-ten rule and the exact exceptions to
four-divisibility, including equal powers and exponent zero.

These laws follow from the authors' earlier Jacobi fraction, Part II,
arXiv:2611.12829v1, Theorem 9.2. The alternative proof uses Michael Wallner's
recurrence recorded in OEIS A274644 on 13 February 2024. The block product
was proved by the posing authors in their 2023 Proposition 5.

Theorem 2 gives a finite matrix and linear recurrence modulo every m>=2
coprime to six. Corollary 3 gives A_j=1 (mod 35) for odd j and 6 for even j,
and A_j=36*6^(j-1)+35*C_j (mod 140) for j>=2. Corollaries 4–7 give mixed
overlap congruences, an exact separator criterion, Catalan-family valuation
consequences, and the transfer to the equivalent patterns 52341, 25314 and
41352. The mixed-overlap product is asserted for 14325, not for 25314.

## Contents

- `note.tex`: self-contained proof and bibliography.
- `note.pdf`: rendered note.
- `metadata.json`: creators, license, scientific summary and pinned Lean-source link.
- `coverage.json`: source mapping, prior results and mathematical verification scope.
- `verify.py`: independent finite checks using Python's standard library.
- `verification.json`: deterministic output of `verify.py`.

## Build

Use Tectonic 0.17.0:

```sh
tectonic --untrusted note.tex
```

Alternatively, use a standard LaTeX distribution:

```sh
pdflatex -interaction=nonstopmode -halt-on-error note.tex
pdflatex -interaction=nonstopmode -halt-on-error note.tex
```

The bibliography is embedded in the TeX source; no BibTeX step is required.

## Check

```sh
python3 verify.py > verification.json
```

All arithmetic is exact. Literal permutations are enumerated for j=1,2;
comparison-poset linear extensions for j=1 through 8; transfer counts and
the congruences for j=1 through 256; the period-five functions for k=0
through 99; and the formal-series identities through degree 64 for
k=0 through 16. All 15 overlap words of sizes one and two with at most
four windows are also counted independently, and the Jacobi continued
fraction is checked modulo twenty through degree 64. All modular models
for m<=100 coprime to six are checked through j=256, as are the modulo-35
and modulo-140 laws. The modulo-25 rational series is checked through
degree 255. Separator existence is checked for all 127 overlap words with
at most seven windows and for pure overlap two through sixteen windows.
The seven Catalan pairs and four equivalent patterns are independently
counted through eight windows.

The written proofs establish the stated all-parameter results. The pinned Lean declarations and their precise mathematical
scope are listed in the note and `coverage.json`.
The original single-block modulo-ten and four-divisibility laws have the
literal packing count as their formal conclusion. The transfer-row modulo-five
and Catalan modulo-four identities are linked to this count by the grid
equivalence. The finite modular models, stronger congruences, mixed products,
separator criterion and additional pattern-family claims have written proofs.

## Creators and license

Haobo Ma (ChronoAI Pte Ltd, Singapore), followed by Wenlin Zhang
(National University of Singapore; ChronoAI Pte Ltd, Singapore).
Licensed under CC BY 4.0 (`cc-by-4.0`).
