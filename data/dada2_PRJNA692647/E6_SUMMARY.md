## E6 -- DADA2/ASV re-analysis of PRJNA692647 (ocular Demodex+/-): conclusions

Faithful DADA2 1.30 pipeline (cutadapt primer removal -> filterAndTrim ->
learnErrors -> dada -> mergePairs -> removeBimeraDenovo -> SILVA v138.1
assignTaxonomy + addSpecies). 31 samples (14 demodex / 17 control), 2x250 V4
515F/806R. 8,551 ASVs, 97.4% reads retained, ~94% pairs merged, merged length
peak 253 bp.

## Headline: the k-mer S. alvi lead does NOT survive ASV resolution

- Zero ASVs are classified as Snodgrassella. SILVA v138.1 contains the genus
  (25 refs) + Neisseriaceae (1,795), so this is genuine absence, not a
  reference gap.
- Applying our own S. alvi diagnostic 31-mer panel directly to the ASVs
  (E6b bridge) flags only non-Snodgrassella Neisseriaceae -- Neisseria,
  Eikenella, Uruburuella, Bergeriella -- plus Methylotenera/Psychrobacter.
- Sequence identity of the top panel-hit ASVs to the S. alvi 16S reference is
  92.5-94.5% over V4 (E6c) -- inter-genus distance, not the ~99% a true
  S. alvi would show.
- Conclusion: the preliminary S. alvi signal (k-mer AUC 0.65) was a
  conserved-region cross-mapping artefact onto abundant commensal
  Neisseriaceae. This is the E5 genus-level limitation, now demonstrated
  quantitatively. The closed-reference k-mer panel is not species-specific for
  these taxa at 16S.

## The defensible signal that survives: a Neisseria ASV

- ASV43 (SILVA Family Neisseriaceae, Genus Neisseria; 4,367 reads, the most
  abundant Neisseriaceae variant) is enriched in Demodex+ vs control:
  AUC 0.697, one-sided p = 0.032 -- the strongest, cleanest, species-resolved
  separation in the study.
- This relocates the ocular-surface Demodex-associated signal from S. alvi to
  Neisseriaceae / Neisseria. It does NOT corroborate the Homey et al. (JID 2025)
  facial S. alvi observation; body site differs (ocular here vs facial there),
  and at the ocular surface the signal is Neisseria, not S. alvi.

## Other taxa

- Bacillus panel is also partly non-specific (k-mers hit Anoxybacillus,
  Enterococcus, Geobacillus, Oscillibacter, Ligilactobacillus). Best genuine
  Bacillaceae ASV40 (Bacillus): AUC 0.643, p = 0.092 -- weak/borderline.
- Bartonella panel hits Rhizobiaceae (Allorhizobium-...-Rhizobium) and other
  Alphaproteobacteria -- confirms the previously-identified cross-mapping;
  Bartonella stays dropped.

## What this means for the manuscript

The result is stronger and more honest as a methods cautionary tale + a new lead:
1. Closed-reference k-mer panels for Demodex-associated taxa are not
   species-specific at 16S; they cross-map within family/order and can
   manufacture a spurious genus-level "signal" (here, S. alvi).
2. ASV/DADA2 resolution is required and overturns the k-mer attribution.
3. The real, exact-sequence Demodex+ association at the ocular surface is a
   Neisseria ASV (AUC 0.70, p = 0.03) -- the recommended lead going forward.
