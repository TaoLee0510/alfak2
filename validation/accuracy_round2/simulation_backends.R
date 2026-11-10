# Independent exact continuous-time finite-domain branching simulator.
simulate_exact_division_abm <- function(karyotypes, fitness, beta, dt, initial_probability,
                                         initial_population, seed, max_events = 10000000L) {
  stopifnot(nrow(karyotypes) == length(fitness), length(fitness) == length(initial_probability),
            all(fitness >= 0), beta >= 0, beta <= 1, dt > 0, initial_population > 0)
  set.seed(seed)
  labels <- apply(karyotypes, 1L, paste, collapse = ".")
  population <- as.numeric(stats::rmultinom(1, initial_population, initial_probability))
  initial <- population
  elapsed <- 0; events <- 0L
  repeat {
    rates <- population * fitness
    total_rate <- sum(rates)
    if (total_rate <= 0) break
    next_time <- elapsed + stats::rexp(1, rate = total_rate)
    if (next_time > dt) break
    elapsed <- next_time
    parent <- sample.int(length(rates), 1L, prob = rates)
    x <- karyotypes[parent, ]
    population[parent] <- population[parent] - 1
    errors <- stats::rbinom(1, sum(x), beta)
    if (errors == 0L) population[parent] <- population[parent] + 2 else {
      copies <- rep(seq_along(x), x)
      picked <- copies[sample.int(length(copies), errors)]
      signs <- sample(c(-1L, 1L), errors, replace = TRUE)
      delta <- tabulate(picked[signs > 0], length(x)) - tabulate(picked[signs < 0], length(x))
      for (daughter in list(x + delta, x - delta)) {
        index <- match(paste(daughter, collapse = "."), labels)
        if (!is.na(index)) population[index] <- population[index] + 1
      }
    }
    events <- events + 1L
    if (events > max_events) stop("Exact ABM exceeded max_events")
    if (sum(population) == 0) stop("Exact ABM extinction")
  }
  list(pop0 = initial, pop1 = population, n_steps = NA_integer_, events = events)
}
