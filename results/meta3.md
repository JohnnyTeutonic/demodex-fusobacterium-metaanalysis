# Cross-cohort genus meta-analysis (ocular Demodex vs control)

Cohorts: discovery, rep2020, shotgun2022. Prevalence filter: genus present in >=50% of samples in EVERY cohort.

Genera passing prevalence filter in all cohorts: 38.

Meta-analysis: Stouffer combined Z of per-cohort Mann-Whitney signed z, weighted by sqrt(n). Concordant = same sign in all cohorts. Combined p is two-sided; Bonferroni over the tested genera.


## All prevalence-filtered genera ranked by meta-analytic |Z|

| genus | dir | combined Z | combined p | Bonferroni p | BH q | concordant | AUC[discovery] | AUC[rep2020] | AUC[shotgun2022] |
|---|---|---|---|---|---|---|---|---|---|
| Fusobacterium | down | -2.93 | 0.0034 | 0.131 | 0.088 | yes | 0.332 | 0.305 | 0.360 |
| Bradyrhizobium | up | +2.83 | 0.0046 | 0.176 | 0.088 | yes | 0.668 | 0.520 | 0.844 |
| Brevundimonas | down | -2.11 | 0.0346 | 1.000 | 0.325 | yes | 0.366 | 0.452 | 0.291 |
| Delftia | up | +2.00 | 0.0453 | 1.000 | 0.325 | no | 0.420 | 0.468 | 0.989 |
| Novosphingobium | up | +1.86 | 0.0631 | 1.000 | 0.325 | yes | 0.735 | 0.539 | 0.578 |
| Pseudomonas | up | +1.85 | 0.0641 | 1.000 | 0.325 | yes | 0.672 | 0.586 | 0.578 |
| Neisseria | down | -1.79 | 0.0729 | 1.000 | 0.325 | no | 0.689 | 0.321 | 0.222 |
| Veillonella | down | -1.79 | 0.0737 | 1.000 | 0.325 | yes | 0.496 | 0.451 | 0.245 |
| Acinetobacter | up | +1.77 | 0.0770 | 1.000 | 0.325 | no | 0.651 | 0.395 | 0.818 |
| Aeromonas | down | -1.57 | 0.1166 | 1.000 | 0.402 | no | 0.500 | 0.345 | 0.400 |
| Bacillus | up | +1.49 | 0.1364 | 1.000 | 0.402 | yes | 0.538 | 0.626 | 0.585 |
| Halomonas | up | +1.48 | 0.1380 | 1.000 | 0.402 | no | 0.723 | 0.393 | 0.705 |
| Corynebacterium | up | +1.40 | 0.1629 | 1.000 | 0.402 | yes | 0.685 | 0.555 | 0.520 |
| Bacteroides | up | +1.38 | 0.1672 | 1.000 | 0.402 | no | 0.483 | 0.633 | 0.607 |
| Stenotrophomonas | down | -1.38 | 0.1684 | 1.000 | 0.402 | yes | 0.382 | 0.436 | 0.433 |
| Cutibacterium | down | -1.37 | 0.1693 | 1.000 | 0.402 | no | 0.538 | 0.367 | 0.375 |
| Helicobacter | up | +1.29 | 0.1985 | 1.000 | 0.427 | yes | 0.506 | 0.579 | 0.633 |
| Prevotella | up | +1.28 | 0.2022 | 1.000 | 0.427 | no | 0.502 | 0.750 | 0.425 |
| Sphingomonas | up | +1.12 | 0.2633 | 1.000 | 0.509 | no | 0.282 | 0.557 | 0.840 |
| Paracoccus | down | -1.10 | 0.2693 | 1.000 | 0.509 | yes | 0.441 | 0.386 | 0.491 |
| Lysinibacillus | up | +1.08 | 0.2815 | 1.000 | 0.509 | yes | 0.601 | 0.505 | 0.600 |
| Paenarthrobacter | down | -1.04 | 0.2990 | 1.000 | 0.516 | no | 0.483 | 0.602 | 0.198 |
| Rothia | up | +0.94 | 0.3472 | 1.000 | 0.527 | no | 0.729 | 0.590 | 0.347 |
| Lactobacillus | up | +0.94 | 0.3491 | 1.000 | 0.527 | no | 0.408 | 0.752 | 0.447 |
| Streptococcus | up | +0.93 | 0.3521 | 1.000 | 0.527 | no | 0.664 | 0.617 | 0.376 |
| Flavobacterium | down | -0.90 | 0.3699 | 1.000 | 0.527 | no | 0.395 | 0.343 | 0.625 |
| Enterococcus | down | -0.89 | 0.3744 | 1.000 | 0.527 | no | 0.437 | 0.383 | 0.538 |
| Mycobacterium | up | +0.84 | 0.3992 | 1.000 | 0.542 | no | 0.412 | 0.550 | 0.676 |
| Sphingobacterium | up | +0.81 | 0.4186 | 1.000 | 0.549 | no | 0.374 | 0.383 | 0.920 |
| Streptomyces | up | +0.78 | 0.4359 | 1.000 | 0.552 | no | 0.605 | 0.467 | 0.589 |
| Chryseobacterium | up | +0.58 | 0.5597 | 1.000 | 0.686 | no | 0.424 | 0.379 | 0.829 |
| Nocardioides | up | +0.52 | 0.6042 | 1.000 | 0.709 | no | 0.450 | 0.550 | 0.582 |
| Haemophilus | up | +0.50 | 0.6159 | 1.000 | 0.709 | yes | 0.508 | 0.536 | 0.538 |
| Micrococcus | down | -0.42 | 0.6769 | 1.000 | 0.757 | no | 0.504 | 0.486 | 0.436 |
| Actinomyces | down | -0.33 | 0.7377 | 1.000 | 0.801 | no | 0.603 | 0.436 | 0.422 |
| Massilia | down | -0.27 | 0.7874 | 1.000 | 0.831 | no | 0.349 | 0.652 | 0.404 |
| Lactococcus | down | -0.22 | 0.8225 | 1.000 | 0.845 | no | 0.643 | 0.335 | 0.525 |
| Staphylococcus | up | +0.15 | 0.8791 | 1.000 | 0.879 | no | 0.571 | 0.424 | 0.556 |

## Concordant hits surviving multiple-testing correction

_None survive Bonferroni/BH at the current cohort count._

## Leave-one-cohort-out robustness (top concordant leads)

Combined Z (and two-sided p) recomputed after dropping each cohort. A robust signal keeps its sign and comparable magnitude throughout.

| genus | full Z | drop discovery | drop rep2020 | drop shotgun2022 |
|---|---|---|---|---|
| Fusobacterium | -2.93 | -2.45 (p=0.014) | -2.06 (p=0.039) | -2.64 (p=0.008) |
| Bradyrhizobium | +2.83 | +2.35 (p=0.019) | +3.46 (p=0.001) | +1.20 (p=0.231) |
| Brevundimonas | -2.11 | -1.70 (p=0.089) | -2.31 (p=0.021) | -1.20 (p=0.229) |
| Novosphingobium | +1.86 | +0.81 (p=0.420) | +2.05 (p=0.040) | +1.75 (p=0.080) |
| Pseudomonas | +1.85 | +1.17 (p=0.243) | +1.65 (p=0.099) | +1.74 (p=0.082) |
| Veillonella | -1.79 | -2.08 (p=0.037) | -1.87 (p=0.061) | -0.43 (p=0.666) |
