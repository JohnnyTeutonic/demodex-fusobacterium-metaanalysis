# E3 -- Ivermectin treatment arm (PRJDB18292)

Design: 6 subjects x {lesional, non-lesional} x {pre, post} topical ivermectin,
24 paired skin-swab 16S samples (UC San Diego, 2019). Index from `build_index.py`
(closed-reference diagnostic k-mers, background-subtracted, >=2 k-mers/read).

## Result: NULL / mixed -- the surface proxy does not show a clean mite drop
Paired Wilcoxon (n=6 pairs/site) and exact sign test, mite_index and each
Demodex-associated taxon: **no stratum reaches significance** (all p >= 0.25).

| stratum                     | n | down | up | median delta | wilcoxon p |
|-----------------------------|---|------|----|--------------|------------|
| lesional / mite_index       | 6 | 2    | 4  | +0.0006      | 0.84       |
| non-lesional / mite_index   | 6 | 2    | 4  | +0.0008      | 0.69       |

Per-subject lesional mite_index pre->post (body-site-matched panel): subjects
1,2 fall (-0.013, -0.016); 3,4,6 ~flat-up; **subject 5 rises steeply (+0.029)**
and dominates. Non-lesional is similarly split. (Verdict is unchanged from the
skin-only panel; only the magnitudes shifted.)

## Honest interpretation (this matters for the protocol)
1. **n=6 is severely underpowered**; one outlier (subj 5) swings the mean.
2. **Swab vs follicle.** Demodex lives deep in the pilosebaceous follicle;
   surface swabs under-sample it (stated limitation in the source study and in
   Woo 2020). A *surface* bacterial proxy is therefore a weak readout of a
   *follicular* mite load.
3. **Mite lysis confound.** Killing mites can release their bacterial cargo
   (B. oleronius, Snodgrassella, etc.) onto the surface, so swab-detectable
   Demodex-associated bacteria could transiently *rise* as live mites fall --
   exactly the kind of mixed signal seen here.
4. Net: the surface 16S proxy is **NOT a validated substitute for a direct mite
   count**. This actually *reinforces* the clinical protocol point -- grading
   treatment needs a proper standardized skin-surface biopsy / follicular mite
   count, not a surface swab or its molecular shadow.

## Next
Validation arm (PRJNA692647, Demodex+ vs control): the decisive test of whether
the index carries any Demodex signal at all when sampling is designed around the
mite. If it separates there but not here, the swab/follicle gap (point 2) is the
explanation; if it fails there too, the bacterial proxy is simply too weak.
