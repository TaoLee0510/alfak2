# A corrected rule for initial_state = "auto": step 0

Approved by the user on 2026-11-16 after the paper benchmark of alfak2
827c4b5 showed that `"auto"` switched to free initial abundances at time
points one passage apart even at the true missegregation rate, and lost
accuracy there (`../README.md`, step 4). The aim: a rule that switches where
the assumed rate is wrong or calls are noisy, and not where the cloud is large
for other reasons. Step 0 was to look, without new rules, for what tells the
two apart, in two ways: (1) the paper benchmark's results directly, and (2)
simulations of the paper's population size and age generated here.

## Data

- **(1) Paper benchmark** (`paper_features.R`, read only): every
  two-time-point ABM and missegregation input of round 05a8855, the level of
  the `"auto"` check from round 827c4b5, and the benchmark's paired effect of
  the switches against the truth (`vs_truth_fits.tsv`). The benchmark fits
  at the true rate, so these inputs only show where switching must not
  happen.
- **(2) Simulations**: the study's six earlier sets (dev, short, diverse and
  their confirmation sets) and a new paper-sized set, `paper3` (seed base
  104000): populations growing from 2e5 to 2e7 cells per passage, first time
  point at day 120, intervals of 6, 18 and 42 days, 1000 cells per time point,
  wavelengths 0.2–1.6, five landscapes each, under three conditions: the
  assumed rate is the true one (`rate1`), the true rate is 8 times the assumed
  one (`rate8`), or 2% of chromosomes are called a copy off (`call2`). Fits
  `base`, `free`, `level_dt` and, new, `level_age`: the level measured against
  the balance over the population's age (history = the first time point's
  day). 1416 fits (array 21984600).
- Switching "helped" where free lowered RMSE at distance one by more than 5%
  or raised the landscape correlation by more than 0.05, "hurt" in the
  opposite case.

Tables: `results/rule_step0/` and `results/logs/rule_step0.log`
(`step0.R`, `interval_rule.R`).

## Findings

**The cloud's diversity does not mark a wrong rate.** The paper's
missegregation simulations, fitted at their true rate, are as diverse as the
study's wrong-rate sets (domain coverage 0.84, a fifth of the cells on
karyotypes seen once), and switching hurt 110 of 119 of them.

**The level against the balance over the population's age does not separate
them either** (AUC 0.52 between helped and hurt inputs). In old populations
most neighbours are less fit than the majors, so the cloud stays far below
the neutral balance over the whole age, whatever the rate (level −5 to −8 on
the diverse sets at day 1440). Within the paper-sized set (one age) it
separates the conditions only partly.

**The level at the interval separates best but depends on the interval**
(AUC 0.75). At the true rate its median falls from 1.5 at 6 days to 0.0 at 42
days; with an 8 times higher rate from 2.9 to 1.0. A threshold that falls
with the interval, `level > 1.75 − 0.25 log(dt / 6)`, chosen on these same data
(exploratory), switches:

| Inputs | Switched | Of the switched, switching hurt |
| --- | --- | --- |
| Paper benchmark ABM, one passage apart | 2% | (all switches there hurt) |
| Paper benchmark missegregation, one passage apart | 6% | (all switches there hurt) |
| Development scenarios | 0–3% | |
| Short intervals, 6 days | 20% | 50–88% |
| Paper-sized, true rate, 6 days | 40% | 63% |
| Paper-sized, 8 times the rate | 45–75% | 0–20% |
| Paper-sized, 2% miscalled | 15–100% | 0–20% |
| Diverse, wrong rate or miscalls | 100% | 4–13% |
| Diverse, true rate, old population on a smooth landscape | 92% | 55–91% |

It catches half of the inputs where switching helps and switches 14% of those
where it hurts. The failures are populations whose cloud is large for
reasons the counts share with a wrong rate: an old population on a smooth
landscape, and young paper-sized populations sampled six days apart.

## Conclusion

Two time points do not tell a wrong missegregation rate or miscalled cells
from a population whose cloud is large for other reasons. This is the
limit already met when choosing the rate itself (`../../beta_choice/`).
The rule would fail the plan's criteria (no loss on short intervals).

On the real data the question matters less than it seemed: on the Salehi
passages, forecasts with free initial abundances are no better than with the
balance prior (`../../salehi_horizon/`). The default stays `"balance"` and
`"auto"` stays an option with its limitation documented. A rate measured
outside the counts, or a range of rates (`fit_sparse_joint_rates()`), remain
the remedies for a misspecified rate.

## Decision

The user stopped the corrected rule after this step on 2026-11-16, on this
recommendation: steps 1 (fix the rule) and 2 (confirm on new seeds) were not
run, since the best rule already fails on the data it was chosen on.
