# Reusable nested MLE fits on veteran (init = 0, fixed degree for speed)
nested_ph_mle <- function(degree = 5L) {
  f0 <- bpph(
    Surv(time, status) ~ 1,
    data = veteran,
    approach = "mle",
    degree = degree,
    init = 0
  )
  f1 <- bpph(
    Surv(time, status) ~ karno,
    data = veteran,
    approach = "mle",
    degree = degree,
    init = 0
  )
  list(null = f0, full = f1, degree = degree)
}

nested_po_mle <- function(degree = 5L) {
  f0 <- bppo(
    Surv(time, status) ~ 1,
    data = veteran,
    approach = "mle",
    degree = degree,
    init = 0
  )
  f1 <- bppo(
    Surv(time, status) ~ karno,
    data = veteran,
    approach = "mle",
    degree = degree,
    init = 0
  )
  list(null = f0, full = f1, degree = degree)
}

nested_aft_mle <- function(degree = 5L) {
  f0 <- bpaft(
    Surv(time, status) ~ 1,
    data = veteran,
    approach = "mle",
    degree = degree,
    init = 0
  )
  f1 <- bpaft(
    Surv(time, status) ~ karno,
    data = veteran,
    approach = "mle",
    degree = degree,
    init = 0
  )
  list(null = f0, full = f1, degree = degree)
}

# Shared MLE fit for information-criterion consistency checks
mle_ph_fit <- function(degree = 5L) {
  bpph(
    Surv(time, status) ~ karno + factor(celltype),
    data = veteran,
    approach = "mle",
    degree = degree,
    init = 0
  )
}
