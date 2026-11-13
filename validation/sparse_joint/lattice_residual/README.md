# Few frequent karyotypes: is a lattice-correlated residual (change C) worth it?

The ALFA-K paper benchmark (`alfak2-paper-benchmark`) found that on the most
rugged landscapes (wavelength 0.2, three or four frequent karyotypes) a
frequent karyotype held out by `as_alfakR()`'s cross-validation is predicted
three to four standard deviations too low, with a held-out correlation near
zero. Change C would replace the independent residual of
`fit_sparse_joint()` (log birth rate minus trend) by one correlated along the
karyotype lattice, so that unobserved karyotypes borrow the deviations of
their neighbours from the trend. This study first checks whether the data
hold anything for such a residual to borrow (C0, C1), before implementing it.

## Data (C0)

The development simulator of `benchmark/sparse_two_timepoint` (sourced, not
modified) at the population sizes of the paper's ABM: each passage grows from
2e5 to 2e7 cells. At the development benchmark's 5e4 to 2e6, wavelength 0.2
gives a median of two frequent karyotypes over 20 pilot trajectories; at the
paper's sizes 3.5, with a largest clone of 0.79 (paper: three to four, about
0.8).

- Wavelengths 0.2 (12 landscapes), 0.4 and 0.8 (8 each); seed base 95000,
  ids `lr_L..`. The four last wavelength-0.2 landscapes (L25-L28) were added
  after the first eight, half of which had one or two frequent karyotypes;
  the twelve have 1, 1, 4, 4, 4, 2, 3, 1, 1, 1, 4 and 2 (at least 10 cells
  over both samples).
- First sample at day 120; intervals of 6, 18 and 42 days; 1000 cells at
  each sample. 84 datasets.
- Fitted with the package as committed (af362af) at its defaults (`base`).

A held-out set (seed base 96000, ids `lh_L..`) is defined for a final
confirmation and has not been generated.

## The weakness reproduces

Cross-validation of `as_alfakR()` (minobs 10) against the truth, medians:

| | this set, w 0.2 | paper ABM, w 0.2 |
| --- | --- | --- |
| folds | 3 | 3-4 |
| held-out slope | 0.13 | 0.06-0.17 |
| held-out correlation | 0.08 | 0.06-0.11 |
| bias of held-out frequent karyotypes (sd) | -4.3 | -2.8 to -3.8 |
| share of error common to a fold | 0.89 | 0.80-0.83 |
| landscape r | 0.23 | 0.20 |

In use, where the frequent karyotypes are in the data, the predictions for
unobserved karyotypes at wavelength 0.2 are flat and conservative (scorer of
the development benchmark; RMSE relative to predicting no differences):

| wavelength | distance | RMSE ratio | Spearman | 95% coverage | bias |
| --- | --- | --- | --- | --- | --- |
| 0.2 | 1 | 0.39 | 0.15 | 1.00 | -0.04 |
| 0.2 | 2 | 0.33 | 0.07 | 1.00 | -0.07 |
| 0.4 | 1 | 0.55 | 0.34 | 0.87 | +0.06 |
| 0.4 | 2 | 0.56 | 0.28 | 0.82 | +0.11 |
| 0.8 | 1 | 0.54 | 0.82 | 0.96 | 0.00 |
| 0.8 | 2 | 0.53 | 0.78 | 0.96 | +0.01 |

A frequent karyotype held out by the cross-validation is treated as never
observed, and a never-observed karyotype next to the remaining majors is, on
average, much less fit than they are: at wavelength 0.2 unobserved one-step
neighbours are 0.47 per day below the count-weighted observed karyotypes
(median over datasets). The model predicts that level, which is
right for an unobserved neighbour (bias -0.04) but not for a karyotype that
selection had made frequent. Most of the cross-validation's bias is this
selection, not an error the model makes in use.

## Screen (C1)

`screen.R` refits each dataset at the history its `base` fit chose and takes
residuals (fitness minus the fitted trend, relative to the majors). A residual
correlated along the lattice moves an unobserved karyotype towards the
residuals of its neighbours; the screen applies such shifts, rho times the
neighbours' mean, to the base predictions and scores them. The sources of the
shift are the fitted residuals of the neighbouring majors, of all observed
neighbours, or of all domain neighbours (pooling the zero counts of unobserved
ones), and, as bounds, the true residuals of the observed neighbours or the
true mean residual of the unobserved neighbours of each major. Geometric mean
RMSE relative to the base fit, at rho = 1 (smaller rho lies in between):

| source | w 0.2, d1 | w 0.4, d1 | w 0.8, d1 | d2, all wavelengths |
| --- | --- | --- | --- | --- |
| majors, fitted | 1.001 | 1.015 | 0.990 | 0.999-1.000 |
| observed neighbours, fitted | 1.000 | 1.014 | 0.991 | 0.999-1.000 |
| all domain neighbours, fitted | 1.000 | 1.008 | 0.995 | 0.999-1.000 |
| observed neighbours, true (bound) | 0.878 | 0.879 | 0.851 | 0.995-0.999 |
| offset around each major, true (bound) | 0.954 | 0.922 | 0.937 | 0.987-0.992 |

The true residual is correlated along the lattice (lag-1 autocorrelation 0.49
at wavelength 0.2, 0.67-0.74 at 0.4, 0.77-0.81 at 0.8), but what the counts
say about it is not. At the majors the trend already fits, so their residuals
are small (sd 0.03, against an error of 0.12-0.19 at distance one). The
residuals the fit gives unobserved domain karyotypes have a standard
deviation of 0.001-0.003 against a true 0.035-0.16, and correlate with the
truth at 0.02-0.08: a zero count does not tell which unobserved neighbour is
fitter. Observed non-major karyotypes, mostly single cells, are where the
bound's gain comes from, and their own fitness is what one or two cells
cannot pin down.

## Conclusion

Change C, in any form that draws on the fitted model, changes the error of
unobserved karyotypes by less than 2% in either direction; it was not
implemented. The second weakness has two parts, neither of which C
addresses: the cross-validation's selection of frequent karyotypes, and, at
wavelength 0.2, the absence of information about unobserved karyotypes in
two samples, which the model reports as flat predictions with full coverage.

A finding outside C: at wavelength 0.4 in these young populations the model
overestimates unobserved karyotypes (bias +0.06 at distance one, +0.11 at
two) and its intervals undercover (0.87 and 0.82; 0.83 and 0.73 at 6-day
intervals). On the development benchmark (first sample at day 1440) coverage
at wavelength 0.4 was 94-98%.

## Running

```sh
R=validation/sparse_joint/lattice_residual/run.R
Rscript $R --stage=truth --workers=6
Rscript $R --stage=datasets --workers=6
Rscript $R --stage=fit --methods=zero,base --workers=6
Rscript $R --stage=score --workers=6
Rscript validation/sparse_joint/lattice_residual/screen.R --workers=6
```

The stage code is that of `../short_interval` (`setup.R`, `stages.R`),
sourced unchanged. Results go to `results/` (not tracked); the library of
af362af is built under `results/lib/`, never into the source tree.
