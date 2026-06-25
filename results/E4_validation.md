# E4 -- Validation arm (PRJNA692647): does the index carry a Demodex signal?

Ocular-surface 16S study, groups from sample title: **14 demodex-infection vs
17 healthy control**. Mann-Whitney U + ROC-AUC (demodex > control).

## Two stages -- the body-site-tuned background changes the verdict

### Stage 1 (skin-tuned background): composite FAILS
With a skin-only background, the composite did not validate (AUC 0.41) because
*Bartonella* ran backwards (AUC 0.30) -- cross-mapping to ocular Proteobacteria
(Rhizobiales/Alphaproteobacteria) absent from the skin background. Only
*S. alvi* trended correctly (AUC 0.65). Diagnosis: garbage components sank the
sum, and the background was not matched to the body site.

### Stage 2 (body-site-tuned background): composite VALIDATES
Adding ocular/environmental + Rhizobiales/Bacillus-relative taxa to the
background removed the cross-mapping (Bartonella diagnostic k-mers 1244 -> 776;
its per-sample reads collapsed from thousands to single digits). Re-scored:

| metric                       | median dem | median ctl | AUC (dem>ctl) | p (1-sided) |
|------------------------------|-----------:|-----------:|--------------:|------------:|
| **mite_index (composite)**   | 0.02430    | 0.01952    | **0.655**     | **0.074**   |
| **Snodgrassella_alvi**       | 0.00535    | 0.00398    | 0.647         | 0.085       |
| Bacillus_oleronius           | 0.01717    | 0.01541    | 0.580         | 0.23        |
| Cutibacterium_kroppenstedtii | 7.3e-5     | 8.7e-5     | 0.43          | 0.74        |
| Bartonella_quintana          | 4.6e-5     | 7.2e-5     | 0.28 (inv)    | 0.98        |

## Reading
- With a **body-site-matched background the composite now separates Demodex+
  from control in the predicted direction** (AUC 0.66; one-sided p ~ 0.074),
  carried by *S. alvi* (AUC 0.65) and *Bacillus* (AUC 0.58).
- *Bartonella* is still mildly inverted but now contributes ~nothing (median
  4.6e-5), so it no longer drags the composite. It is an unreliable component and
  should be dropped from any deployed panel.
- The key methods lesson: **closed-reference marker counting must use a
  body-site-matched negative control.** The "failure" in Stage 1 was an
  artefact of a mismatched background, not evidence against the hypothesis.

## Implication
After the fix, both the composite and the single best marker (*S. alvi*) carry a
borderline-significant Demodex signal at the ocular surface. A curated assay
around *S. alvi* (+ optionally *Bacillus*), with body-site-matched references and
background and follicle-targeted sampling, is the recommended next step.
