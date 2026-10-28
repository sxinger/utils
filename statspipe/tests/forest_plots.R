# Run from the repository root: Rscript statspipe\tests\forest_plots.R
source(file.path("statspipe", "scripts", "statspipe.R"))

check_forest_plots <- function() {
  set.seed(42)
  data <- data.frame(
    AGE = rnorm(300, 50, 10),
    SEX = factor(rep(c("Female", "Male"), 150)),
    event = rbinom(300, 1, 0.4)
  )
  outcomes <- setNames(
    rep(list(list(type = "logistic", variable = "event")), 5),
    paste("Outcome", seq_len(5))
  )
  args <- list(
    data = data, outcomes = outcomes,
    covariates = c("AGE", "SEX"), selected = c("AGE", "SEX"),
    forest_ncol = 3, forest_annotate = FALSE
  )
  baseline <- do.call(run_statspipe, args)
  labels <- c(AGE = "Age (years)", SEX = "Sex")
  labeled <- do.call(run_statspipe, c(args, list(forest_labels = labels)))
  stopifnot(
    inherits(labeled$forest$plot, "patchwork"),
    length(labeled$forest$plots) == 5L,
    identical(baseline$multivariable$combined, labeled$multivariable$combined),
    identical(baseline$forest$data$term, labeled$forest$data$term),
    all(labeled$forest$data$term_label[labeled$forest$data$term == "AGE"] ==
          "Age (years)"),
    all(labeled$forest$data$term_label[labeled$forest$data$term == "SEXMale"] ==
          "Sex: Male")
  )

  check_rows <- function(plot, counts) {
    layout <- patchwork::patchworkGrob(plot)$layout
    bounds <- lapply(seq_along(counts), function(i) {
      panels <- layout[
        grepl(paste0("^panel-[0-9]+-[0-9]+-", i, "$"), layout$name), ,
        drop = FALSE
      ]
      stopifnot(nrow(panels) == counts[i])
      c(left = min(panels$l), right = max(panels$r))
    })
    stopifnot(all(vapply(bounds, identical, logical(1), bounds[[1L]])))
  }
  check_rows(labeled$forest$plot, c(3L, 2L))

  args$forest_ncol <- 2
  args$forest_source <- "univariate"
  args$forest_annotate <- TRUE
  override <- do.call(run_statspipe, c(
    args, list(forest_labels = c(labels, SEXMale = "Male vs female"))
  ))
  check_rows(override$forest$plot, c(2L, 2L, 1L))
  stopifnot(all(
    override$forest$data$term_label[override$forest$data$term == "SEXMale"] ==
      "Male vs female"
  ))
  for (plot in override$forest$plots) {
    axis_labels <- ggplot2::ggplot_build(plot)$layout$panel_params[[1L]]$y$get_labels()
    stopifnot(all(c("Age (years)", "Male vs female") %in% axis_labels))
  }

  args$forest_ncol <- NULL
  args$forest_nrow <- 2
  by_rows <- do.call(run_statspipe, args)
  check_rows(by_rows$forest$plot, c(3L, 2L))

  args$forest_nrow <- NULL
  default_layout <- do.call(run_statspipe, args)
  check_rows(default_layout$forest$plot, c(3L, 2L))
  args$outcomes <- outcomes[1:4]
  complete_grid <- do.call(run_statspipe, args)
  stopifnot(
    !inherits(complete_grid$forest$plot, "patchwork"),
    nrow(ggplot2::ggplot_build(complete_grid$forest$plot)$layout$layout) == 4L
  )
  args$outcomes <- outcomes[1:3]
  complete <- do.call(run_statspipe, args)
  stopifnot(!inherits(complete$forest$plot, "patchwork"))
  args$outcomes <- outcomes[1]
  single <- do.call(run_statspipe, args)
  stopifnot(length(single$forest$plots) == 1L)

  for (invalid in list(c("Age"), c(AGE = ""), c(AGE = NA_character_),
                       c(AGE = "Age", AGE = "Duplicate"))) {
    error <- tryCatch(
      do.call(run_statspipe, c(args, list(forest_labels = invalid))),
      error = identity
    )
    stopifnot(inherits(error, "error"),
              grepl("'forest_labels'", conditionMessage(error), fixed = TRUE))
  }

  files <- c(tempfile(fileext = ".pdf"), tempfile(fileext = ".png"))
  on.exit(unlink(files), add = TRUE)
  for (file in files) {
    ggplot2::ggsave(file, labeled$forest$plot, width = 12, height = 8)
  }
  stopifnot(all(file.info(files)$size > 0))
  message("Forest plot checks passed.")
}

check_forest_plots()
