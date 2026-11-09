trial_correlation_fixture <- function() {
  traces <- array(c(
    0, 1, 2, 4, 0, 3, 1, 0, 0, 2, 3, 3,
    1, 0, 3, 4, 2, 1, 1, 2, 0, 2, 0, 1
  ), dim = c(6L, 4L, 1L))
  metadata <- data.frame(
    Experiment = factor(rep("trial-correlation", 4L)),
    Step = ordered(rep(1, 4L)),
    RecordingID = rep(7L, 4L),
    RunUID = seq_len(4L)
  )
  EPhysData::newEPhysContinuous(
    traces, 0:5, metadata, "A", data.frame(row.names = "A")
  )
}

test_that("the public trial correlation test defaults to CPU", {
  x <- trial_correlation_fixture()
  result <- TrialCorrelationTest(x, n_perm = 20L, seed = 7L, parallel = FALSE, progress = FALSE)
  expect_identical(result, TrialCorrelationTest(
    x, n_perm = 20L, seed = 7L, parallel = FALSE, progress = FALSE, GPU = FALSE
  ))
  expect_equal(result$n_trials, 4L)
  expect_equal(result$n_perm_used, 20L)
  expect_error(TrialCorrelationTest(x, GPU = NA), "requires GPU to be TRUE or FALSE")

  if (!requireNamespace("torch", quietly = TRUE) ||
      !torch::torch_is_installed() || !torch::cuda_is_available()) {
    expect_error(TrialCorrelationTest(x, GPU = TRUE), "GPU = TRUE requires")
  }
})

test_that("the torch batches agree with CPU permutations and preserve the RNG", {
  skip_if_not_installed("torch")
  skip_if_not(torch::torch_is_installed(), "R torch runtime is unavailable")
  x <- t(trial_correlation_fixture()@Data[, , 1L])
  for (stop_upper in c(1, 0)) {
    cpu <- trial_cor_pvalue_perm(
      x, n_perm = 71L, seed = 27L, min_perm = 3L, stop_lower = 0,
      stop_upper = stop_upper, return_perm = TRUE
    )
    cpu_seed <- .Random.seed
    torch_result <- trial_cor_pvalue_perm(
      x, n_perm = 71L, seed = 27L, min_perm = 3L, stop_lower = 0,
      stop_upper = stop_upper, return_perm = TRUE, gpu_device = "cpu"
    )
    expect_equal(torch_result, cpu, tolerance = 1e-10)
    expect_identical(torch_result$p_value, cpu$p_value)
    expect_identical(.Random.seed, cpu_seed)
    expect_identical(torch_result$n_perm_used, if (stop_upper == 0) 3L else 71L)
  }

  tied <- rbind(c(0, 1), c(2, 0), c(0, 3))
  cpu <- trial_cor_pvalue_perm(tied, n_perm = 4L, seed = 27L, return_perm = TRUE)
  torch_result <- trial_cor_pvalue_perm(
    tied, n_perm = 4L, seed = 27L, return_perm = TRUE, gpu_device = "cpu"
  )
  expect_identical(torch_result, cpu)
})

test_that("the public CUDA path matches serial CPU with the same seed", {
  skip_if_not_installed("torch")
  skip_if_not(torch::torch_is_installed(), "R torch runtime is unavailable")
  skip_if_not(torch::cuda_is_available(), "CUDA is unavailable")
  x <- trial_correlation_fixture()
  args <- list(X = x, n_perm = 71L, seed = 27L, min_perm = 3L,
               stop_lower = 0, stop_upper = 0, parallel = FALSE, progress = FALSE)
  cpu <- do.call(TrialCorrelationTest, args)
  cpu_seed <- .Random.seed
  cuda <- do.call(TrialCorrelationTest, c(args, list(GPU = TRUE)))
  expect_identical(cuda, cpu)
  expect_identical(.Random.seed, cpu_seed)
})
