# E5 -- Method characterisation (sensitivity + specificity of the k-mer panel)

Run: `scripts/validate_panel.py`. Addresses the "approximate, untested" flag by
*measuring* the closed-reference k-mer method instead of asserting it.

- **Sensitivity** = fraction of simulated 250 bp windows over HELD-OUT same-taxon
  16S references (accessions beyond the build set) detected under the >=2-k-mer
  rule.
- **Specificity** = same, applied to congeneric/confamilial CHALLENGE organisms
  that must NOT be called as the target.

| taxon                        | sensitivity | challenge detections (should be ~0) |
|------------------------------|------------:|-------------------------------------|
| Snodgrassella_alvi           | 0.84        | N. gonorrhoeae 0.00; **Kingella kingae 0.47** |
| Bacillus_oleronius           | 0.90        | B. subtilis 0.00; **Heyndrickxia coagulans 0.66** |
| Cutibacterium_kroppenstedtii | 0.98        | C. acnes 0.00; **C. diphtheriae 0.43** |
| Bartonella_quintana          | 0.74        | **B. henselae 0.67**; Brucella abortus 0.00 |

## Interpretation
- **Good sensitivity** (0.74-0.98): the panels reliably detect their own taxon.
- **Species-specific vs distant relatives** (0.00 against N. gonorrhoeae,
  B. subtilis, C. acnes, Brucella -- several of which are in the background).
- **Genus/family-level cross-reactivity vs close congeners** (Kingella,
  Heyndrickxia, C. diphtheriae, B. henselae score 0.43-0.67). So the index
  resolves at roughly **genus level, not strict species** -- now quantified.

## Consequence for the conclusions
- The lead marker *S. alvi* reads clean against *Neisseria* (0.00); its main
  leak is to *Kingella* (oral/rare, not a dominant skin or ocular taxon), so the
  practical impact on the Demodex signal is limited but must be stated.
- *Bartonella quintana* cannot be separated from *B. henselae* by 16S k-mers
  (0.67) -- another reason to drop Bartonella from a deployed panel.
- Full species resolution would require ASV/oligotyping (DADA2) or a longer
  marker; that remains the stated upgrade path. The method is now *characterised*
  rather than merely *approximate*.
