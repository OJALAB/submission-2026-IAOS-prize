library(stringr)
library(doFuture)
library(progressr)
library(doRNG)

source("code/functions.R")
load("data-raw/data.RData")

data <- data[str_sub(poz_kodZawodu, 1, 1) != "0"]
data <- data[, c("poz_kodZawodu", "poz_lWolnychMiejsc")]
data[, poz_kodZawodu := as.factor(str_sub(poz_kodZawodu, 1, 1))]
setnames(data, old = c("poz_kodZawodu", "poz_lWolnychMiejsc"), new = c("code", "vac"))

Pi_values <- c(
  82.8, 11.2,  1.7,  1.7,  2.2,  0.0,  0.4,  0.0,  0.0,
  3.5, 88.2,  4.8,  1.6,  1.6,  0.0,  0.3,  0.0,  0.0,
  2.0, 13.9, 68.8,  2.6,  4.4,  0.0,  6.0,  1.2,  1.2,
  1.1, 10.2,  6.0, 71.3,  3.0,  0.0,  0.0,  1.1,  7.2,
  1.5, 4.3, 6.3,  1.0, 85.0,  0.0,  0.8,  0.0,  1.3,
  0.0, 14.3, 0.0,  0.0,  0.0, 71.4,  0.0,  0.0, 14.3,
  0.6, 1.4, 1.7,  0.2,  0.4,  0.0, 90.7,  2.7,  2.3,
  0.0,  0.0,  1.9,  1.2,  0.4,  0.0, 18.2, 77.9,  0.4,
  0.5, 0.0, 2.7,  5.0,  2.7,  2.3,  9.1,  0.9, 76.8
) / 100
Pi <- matrix(Pi_values, nrow = 9, ncol = 9, byrow = TRUE)
rownames(Pi) <- 1:9
colnames(Pi) <- 1:9
Pi <- Pi / rowSums(Pi)

registerDoFuture()
plan(multisession, workers = 10)
on.exit(plan(sequential), add = TRUE)

handlers("txtprogressbar")

iterations <- 100

estimates <- with_progress({
  p <- progressor(steps = iterations)
  set.seed(123)
  foreach(i = 1:iterations, 
          .combine = rbind, 
          .packages = c("data.table", "brglm2")) %dorng% {
            p(sprintf("Iteration: %g", i))
            sim_vac(df = data,
                    p_easy = c(0.8, 0.2, 0.1),
                    p_hard = c(0.1, 0.8, 0.6),
                    prob_hard_vec = c(0.4, 0.3, 0.2, 0.3, 0.4, 0.3, 0.2, 0.3, 0.4),
                    cens_frac = 0.25,
                    missing_frac = 0.25,
                    Pi = Pi,
                    val_sample_size = 1000,
                    max_iter_em = 1000,
                    tol_em = 1e-4,
                    trace_em = FALSE,
                    use_estimated_Pi = TRUE)
          }
})

results <- eval(estimates)
names(results) <- c("total", "by_code")

results_total <- results[["total"]]
results_by_code <- results[["by_code"]]

save(results_total, file = "results/results_total.RData")
save(results_by_code, file = "results/results_by_code.RData")