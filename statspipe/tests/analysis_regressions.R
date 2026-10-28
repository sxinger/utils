# Run from the repository root: Rscript statspipe\tests\analysis_regressions.R
source(file.path("statspipe", "scripts", "statspipe.R"))

check_analysis_regressions <- function() {
  set.seed(87)
  data <- data.frame(
    age = rnorm(200, 50, 10),
    sex = factor(rep(c("Female", "Male"), 100)),
    group = rep(seq_len(40), each = 5),
    response = rnorm(200),
    time = rexp(200),
    status = rep(c(0, 1, 2, 1), 50)
  )
  args <- list(
    data = data, covariates = "age", selected = "age",
    outcomes = list(Survival = list(type = "cox", time = "time", event = "event"))
  )
  args$data$event <- as.integer(data$status == 1)

  survival_result <- do.call(run_statspipe, c(args, list(conf_level = 0.9)))
  stopifnot(
    nrow(survival_result$multivariable$combined) == 1L,
    survival_result$multivariable$combined$n == nrow(data),
    survival_result$univariate$combined$n == nrow(data),
    all(grepl("90% CI", survival_result$forest$data$outcome_label, fixed = TRUE)),
    grepl("90% confidence", survival_result$forest$plot$labels$subtitle, fixed = TRUE),
    survival_result$settings$conf_level == 0.9
  )
  args$outcomes$Survival$method <- "cause_specific"
  args$outcomes$Survival$event <- "status"
  args$outcomes$Survival$event_of_interest <- 1
  cause_specific <- do.call(run_statspipe, args)
  stopifnot(cause_specific$multivariable$combined$n == nrow(data))

  args$outcomes$Survival$method <- "finegray"
  finegray_result <- do.call(run_statspipe, args)
  reference_data <- data
  reference_data$id <- seq_len(nrow(data))
  reference_data$state <- factor(data$status, levels = c(0, 1, 2))
  expanded <- survival::finegray(
    survival::Surv(time, state) ~ age + id, data = reference_data, etype = "1"
  )
  reference <- survival::coxph(
    survival::Surv(fgstart, fgstop, fgstatus) ~ age,
    data = expanded, weights = fgwt, cluster = id, robust = TRUE
  )
  fitted <- finegray_result$multivariable$models$Survival
  stopifnot(
    !is.null(fitted),
    finegray_result$multivariable$combined$n == nrow(data),
    isTRUE(all.equal(unname(stats::coef(fitted)), unname(stats::coef(reference)))),
    isTRUE(all.equal(unname(stats::vcov(fitted)), unname(stats::vcov(reference))))
  )
  args$data$status <- factor(data$status, levels = c(0, 1, 2),
                             labels = c("censored", "target", "competing"))
  args$outcomes$Survival$event_of_interest <- "target"
  args$outcomes$Survival$censor_value <- "censored"
  factor_status <- do.call(run_statspipe, args)
  stopifnot(isTRUE(all.equal(
    unname(stats::coef(factor_status$multivariable$models$Survival)),
    unname(stats::coef(reference))
  )))

  args$outcomes <- list(Response = list(type = "linear", variable = "response"))
  args$cluster <- "group"
  args$data$age[c(2, 19, 88)] <- NA_real_
  previous_options <- options(na.action = "na.omit")
  on.exit(options(previous_options), add = TRUE)
  omitted <- do.call(run_statspipe, args)
  options(na.action = "na.exclude")
  excluded <- do.call(run_statspipe, args)
  stopifnot(
    !is.null(excluded$multivariable$combined),
    isTRUE(all.equal(omitted$multivariable$combined, excluded$multivariable$combined))
  )
  options(previous_options)

  output_dir <- tempfile("statspipe-regression-")
  on.exit(unlink(output_dir, recursive = TRUE), add = TRUE)
  args$cluster <- NULL
  args$strata <- "sex"
  args$save <- TRUE
  args$output_dir <- output_dir
  saved <- do.call(run_statspipe, args)
  for (ext in c("rds", "csv")) {
    path <- file.path(output_dir, "table1",
                      paste0("statistical_analysis_table1_by_strata_sex.", ext))
    stopifnot(file.exists(path), file.info(path)$size > 0)
  }
  reloaded <- readRDS(saved$saved_files$complete_results)
  stopifnot("sex" %in% names(reloaded$table1$by_strata))

  sparse_args <- args
  sparse_args$save <- FALSE
  sparse_args$strata <- "stratum"
  sparse_args$table1_vars <- "category"
  sparse_args$data$stratum <- factor(rep(c("A", "B", "C"), length.out = nrow(data)))
  sparse_args$data$category <- factor(rep(seq_len(17), length.out = nrow(data)))
  set.seed(123)
  sparse_result <- do.call(run_statspipe, sparse_args)
  sparse_table <- sparse_result$table1$by_strata$stratum$table_body
  p_rows <- sparse_table$variable == "category" & sparse_table$row_type == "label"
  stopifnot(
    sum(p_rows) == 1L,
    sparse_table$test_name[p_rows] == "fisher.test",
    is.finite(sparse_table$p.value[p_rows]),
    sparse_table$p.value[p_rows] > 0,
    sparse_table$p.value[p_rows] <= 1
  )

  args$save <- FALSE
  expect_error <- function(input, text) {
    error <- tryCatch(do.call(run_statspipe, input), error = identity)
    stopifnot(inherits(error, "error"),
              grepl(text, conditionMessage(error), fixed = TRUE))
  }
  for (invalid in list(0, 1, NA_real_, Inf, c(0.9, 0.95))) {
    expect_error(c(args, list(conf_level = invalid)), "'conf_level'")
  }
  invalid_args <- args
  invalid_args$outcomes <- setNames(rep(args$outcomes, 2), c("same", "same"))
  expect_error(invalid_args, "'outcomes'")
  invalid_args$outcomes <- list()
  expect_error(invalid_args, "'outcomes'")
  invalid_args <- args
  invalid_args$forest_annotate <- NA
  expect_error(invalid_args, "'forest_annotate'")
  invalid_args <- args
  invalid_args$save <- "yes"
  expect_error(invalid_args, "'save'")
  message("Analysis regression checks passed.")
}

check_analysis_regressions()
