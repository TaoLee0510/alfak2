# Arms of the cloud-level study on the Salehi inputs, as methods alfak2_<arm>
# of run_candidate.R (a copy of the paper benchmark's scripts/run_task.R with
# a hook that adds an arm's arguments to the fit; the hook is the
# alfak2-salehi-variants project's). Sourced inside that script, after the
# benchmark's R/ files, with the input `yi` in the global environment.
#   base      the build at its defaults: must reproduce the round's alfak2 fit
#   level     balance prior with a free level, history among the candidates
#   level_dt  balance prior with a free level, history fixed at the interval
apply_variant <- function(variant) {
  switch(variant,
    base = list(),
    level = list(balance_level = "free"),
    level_dt = list(balance_level = "free", history = diff(input_times(get("yi", envir = globalenv())))),
    stop("unknown variant ", variant))
}
