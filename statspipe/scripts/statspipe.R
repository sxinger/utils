run_statspipe <- function(
    data,
    outcomes,
    covariates,
    selected,
    cohort_filters = NULL,
    strata = NULL,
    cluster = NULL,
    table1_vars = NULL,
    conf_level = 0.95,
    
    # Forest plot options
    forest_source = c("multivariable", "univariate"),
    forest_include = NULL,
    forest_exclude_levels = NULL,
    forest_annotate = TRUE,
    forest_ncol = NULL,
    forest_nrow = NULL,
    ncol = NULL,  
    nrow = NULL,   
    verbose = FALSE,
    
    # Saving options
    save = FALSE,
    output_dir = NULL,
    prefix = "statistical_analysis",
    forest_labels = NULL
) {
  
  # ============================================================
  # 0. PACKAGE CHECKS
  # ============================================================
  
  required_packages <- c(
    "dplyr",
    "broom",
    "survival",
    "gtsummary",
    "ggplot2",
    "rlang"
  )
  
  missing_packages <- required_packages[
    !vapply(
      required_packages,
      requireNamespace,
      quietly = TRUE,
      FUN.VALUE = logical(1)
    )
  ]
  
  if (length(missing_packages) > 0) {
    stop(
      paste0(
        "Please install the following packages first: ",
        paste(missing_packages, collapse = ", ")
      )
    )
  }
  
  forest_source <- match.arg(forest_source)
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || !is.finite(conf_level) ||
      conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be a single number strictly between 0 and 1.")
  }
  for (arg_name in c("forest_annotate", "save")) {
    value <- get(arg_name)
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop("'", arg_name, "' must be TRUE or FALSE.")
    }
  }
  ci_label <- paste0(format(100 * conf_level, trim = TRUE), "%")
  forest_subtitle <- paste0(
    "Points represent estimates and lines represent ", ci_label,
    " confidence intervals"
  )

  if (!is.null(ncol)) {
    if (!is.null(forest_ncol) && !identical(as.numeric(ncol), as.numeric(forest_ncol))) {
      stop("Specify only one value for 'ncol' and 'forest_ncol'.")
    }
    forest_ncol <- ncol
  }
  if (!is.null(nrow)) {
    if (!is.null(forest_nrow) && !identical(as.numeric(nrow), as.numeric(forest_nrow))) {
      stop("Specify only one value for 'nrow' and 'forest_nrow'.")
    }
    forest_nrow <- nrow
  }

  for (arg_name in c("forest_ncol", "forest_nrow")) {
    value <- get(arg_name)
    if (!is.null(value) &&
        (!is.numeric(value) || length(value) != 1L ||
         is.na(value) || !is.finite(value) || value < 1 || value != floor(value))) {
      stop("'", arg_name, "' must be NULL or a positive integer.")
    }
  }
  if (!is.null(forest_labels) &&
      (!is.character(forest_labels) || !is.null(dim(forest_labels)) ||
       is.null(names(forest_labels)) || anyNA(names(forest_labels)) ||
       any(!nzchar(trimws(names(forest_labels)))) ||
       anyDuplicated(names(forest_labels)) || anyNA(forest_labels) ||
       any(!nzchar(trimws(forest_labels))))) {
    stop("'forest_labels' must be NULL or a uniquely named character vector of non-empty labels.")
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("'verbose' must be TRUE or FALSE.")
  }
  log_progress <- function(...) {
    if (isTRUE(verbose)) message(sprintf("[run_statspipe] %s", paste0(..., collapse = "")))
  }
  log_progress("Starting analysis")
  
  
  # ============================================================
  # 1. INPUT CHECKS
  # ============================================================
  
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame.")
  }

  # ------------------------------------------------------------
  # Cohort eligibility is evaluated separately for each outcome.
  # A named list maps outcome names to character expression vectors
  # (AND-combined), or functions returning row masks/data.frames.
  # NULL entries and unspecified outcomes receive the full cohort.
  # A bare vector/function is accepted as a shared legacy filter.
  # ------------------------------------------------------------

  filter_one_cohort <- function(source, filters, outcome_name) {
    if (is.null(filters)) return(source)
    if (is.character(filters)) {
      if (anyNA(filters) || any(!nzchar(trimws(filters)))) {
        stop("'cohort_filters' contains an empty or missing expression.")
      }
      expressions <- tryCatch(
        rlang::parse_exprs(filters),
        error = function(e) stop("Invalid 'cohort_filters': ", conditionMessage(e))
      )
      filtered <- tryCatch(
        dplyr::filter(source, !!!expressions),
        error = function(e) stop("Unable to filter outcome '", outcome_name, "': ", conditionMessage(e))
      )
    } else if (is.function(filters)) {
      filtered <- filters(source)
      if (is.logical(filtered) && length(filtered) == nrow(source)) {
        filtered <- source[!is.na(filtered) & filtered, , drop = FALSE]
      }
    } else {
      stop("Each outcome filter must be NULL, a character vector, or a function.")
    }
    if (!is.data.frame(filtered) || !identical(names(filtered), names(source))) {
      stop("Filter for outcome '", outcome_name, "' must preserve data columns.")
    }
    if (nrow(filtered) == 0L || nrow(filtered) > nrow(source)) {
      stop("Invalid cohort size after filtering outcome '", outcome_name, "'.")
    }
    filtered
  }
  
  if (!is.list(outcomes) || length(outcomes) == 0L) {
    stop(
      "'outcomes' must be a named list."
    )
  }
  
  if (
    is.null(names(outcomes)) ||
    anyNA(names(outcomes)) ||
    any(!nzchar(trimws(names(outcomes)))) ||
    anyDuplicated(names(outcomes))
  ) {
    stop(
      "'outcomes' must be a non-empty, uniquely named list."
    )
  }

  if (is.list(cohort_filters)) {
    if (is.null(names(cohort_filters)) ||
        anyNA(names(cohort_filters)) ||
        any(!nzchar(names(cohort_filters))) ||
        anyDuplicated(names(cohort_filters))) {
      stop("'cohort_filters' must be a uniquely named list keyed by outcome.")
    }
    unknown <- setdiff(names(cohort_filters), names(outcomes))
    if (length(unknown)) {
      stop("Unknown cohort filter outcome(s): ", paste(unknown, collapse = ", "))
    }
  }

  cohort_data <- setNames(vector("list", length(outcomes)), names(outcomes))
  cohort_counts <- data.frame(
    outcome = names(outcomes),
    n_before = nrow(data),
    n_after = NA_integer_,
    n_excluded = NA_integer_,
    stringsAsFactors = FALSE
  )
  for (i in seq_along(outcomes)) {
    nm <- names(outcomes)[i]
    rule <- if (is.list(cohort_filters)) cohort_filters[[nm]] else cohort_filters
    cohort_data[[nm]] <- filter_one_cohort(data, rule, nm)
    cohort_counts$n_after[i] <- nrow(cohort_data[[nm]])
    cohort_counts$n_excluded[i] <- nrow(data) - cohort_counts$n_after[i]
    log_progress("Cohort ", nm, ": ", cohort_counts$n_after[i],
                 "/", nrow(data), " records eligible")
  }
  
  
  # ------------------------------------------------------------
  # Covariates
  # ------------------------------------------------------------
  
  missing_covariates <- setdiff(
    covariates,
    names(data)
  )
  
  if (length(missing_covariates) > 0) {
    warning(
      paste0(
        "Covariates not found and skipped: ",
        paste(
          missing_covariates,
          collapse = ", "
        )
      )
    )
  }
  
  covariates <- intersect(
    covariates,
    names(data)
  )
  
  
  # ------------------------------------------------------------
  # selected may be:
  #
  # selected = c("AGE", "SEX", ...)
  #
  # or
  #
  # selected = list(
  #   Outcome1 = c(...),
  #   Outcome2 = c(...)
  # )
  # ------------------------------------------------------------
  
  if (
    !is.character(selected) &&
    !is.list(selected)
  ) {
    stop(
      "'selected' must be a character vector or named list."
    )
  }
  
  
  # ------------------------------------------------------------
  # Strata
  # ------------------------------------------------------------
  
  if (!is.null(strata)) {
    
    missing_strata <- setdiff(
      strata,
      names(data)
    )
    
    if (length(missing_strata) > 0) {
      
      warning(
        paste0(
          "Strata variables not found and skipped: ",
          paste(
            missing_strata,
            collapse = ", "
          )
        )
      )
    }
    
    strata <- intersect(
      strata,
      names(data)
    )
  }
  
  
  # ------------------------------------------------------------
  # Cluster
  # ------------------------------------------------------------

  if (!is.null(cluster)) {

    if (!is.character(cluster) || length(cluster) != 1) {
      stop("'cluster' must be NULL or a single column name.")
    }

    if (!cluster %in% names(data)) {
      stop(
        paste0(
          "Cluster variable '",
          cluster,
          "' not found in data."
        )
      )
    }
  }


  # ------------------------------------------------------------
  # Table 1 variables
  # ------------------------------------------------------------
  
  if (is.null(table1_vars)) {
    
    if (is.list(selected)) {
      
      selected_all <- unique(
        unlist(
          selected,
          use.names = FALSE
        )
      )
      
    } else {
      
      selected_all <- selected
    }
    
    table1_vars <- unique(
      c(
        covariates,
        selected_all
      )
    )
  }
  
  missing_table1 <- setdiff(
    table1_vars,
    names(data)
  )
  
  if (length(missing_table1) > 0) {
    warning(
      paste0(
        "Table 1 variables not found and skipped: ",
        paste(
          missing_table1,
          collapse = ", "
        )
      )
    )
  }
  
  table1_vars <- intersect(
    table1_vars,
    names(data)
  )
  
  
  # ============================================================
  # HELPER: clean model names
  # ============================================================
  
  clean_name <- function(x) {
    gsub(
      "`",
      "",
      x,
      fixed = TRUE
    )
  }
  
  
  # ============================================================
  # HELPER: outcome specifications
  # ============================================================
  
  get_outcome_info <- function(
    outcome_name,
    spec
  ) {
    
    if (is.null(spec$type)) {
      stop(
        paste0(
          "Outcome '",
          outcome_name,
          "' must specify type."
        )
      )
    }
    
    type <- tolower(
      spec$type
    )
    
    
    # Logistic
    if (
      type %in%
      c(
        "logistic",
        "binary"
      )
    ) {
      
      variable <- spec$variable
      
      if (is.null(variable)) {
        variable <- outcome_name
      }
      
      if (!variable %in% names(data)) {
        stop(
          paste0(
            "Outcome variable '",
            variable,
            "' not found."
          )
        )
      }
      
      return(
        list(
          name = outcome_name,
          type = "logistic",
          variable = variable,
          time = NULL,
          event = NULL,
          measure = "OR"
        )
      )
    }
    
    
    # Poisson
    if (type == "poisson") {
      
      variable <- spec$variable
      
      if (is.null(variable)) {
        variable <- outcome_name
      }
      
      if (!variable %in% names(data)) {
        stop(
          paste0(
            "Outcome variable '",
            variable,
            "' not found."
          )
        )
      }
      
      return(
        list(
          name = outcome_name,
          type = "poisson",
          variable = variable,
          time = NULL,
          event = NULL,
          measure = "IRR"
        )
      )
    }
    
    
    # Linear
    if (
      type %in%
      c(
        "linear",
        "gaussian"
      )
    ) {
      
      variable <- spec$variable
      
      if (is.null(variable)) {
        variable <- outcome_name
      }
      
      if (!variable %in% names(data)) {
        stop(
          paste0(
            "Outcome variable '",
            variable,
            "' not found."
          )
        )
      }
      
      return(
        list(
          name = outcome_name,
          type = "linear",
          variable = variable,
          time = NULL,
          event = NULL,
          measure = "Beta"
        )
      )
    }
    
    
    # Cox survival
    if (
      type %in%
      c(
        "cox",
        "survival"
      )
    ) {
      
      if (
        is.null(spec$time) ||
        is.null(spec$event)
      ) {
        stop(
          paste0(
            "Survival outcome '",
            outcome_name,
            "' requires 'time' and 'event'."
          )
        )
      }
      
      missing_surv <- setdiff(
        c(
          spec$time,
          spec$event
        ),
        names(data)
      )
      
      if (length(missing_surv) > 0) {
        
        stop(
          paste0(
            "Missing survival variables: ",
            paste(
              missing_surv,
              collapse = ", "
            )
          )
        )
      }
      
      method <- if (is.null(spec$method)) {
        "cox"
      } else {
        tolower(spec$method)
      }

      if (!method %in% c("cox", "cause_specific", "finegray")) {
        stop(
          paste0(
            "Unsupported survival method for outcome '",
            outcome_name,
            "': ",
            method
          )
        )
      }

      event_of_interest <- spec$event_of_interest
      censor_value <- if (is.null(spec$censor_value)) 0 else spec$censor_value

      if (method %in% c("cause_specific", "finegray") &&
          is.null(event_of_interest)) {
        stop(
          paste0(
            "Survival outcome '",
            outcome_name,
            "' with method='",
            method,
            "' requires 'event_of_interest'."
          )
        )
      }

      return(
        list(
          name = outcome_name,
          type = "survival",
          variable = NULL,
          time = spec$time,
          event = spec$event,
          method = method,
          event_of_interest = event_of_interest,
          censor_value = censor_value,
          measure = if (method == "finegray") "sHR" else "HR"
        )
      )
    }
    
    
    stop(
      paste0(
        "Unsupported model type: ",
        type
      )
    )
  }
  
  
  # ============================================================
  # HELPER: selected variables by outcome
  # ============================================================
  
  get_selected_variables <- function(
    outcome_name
  ) {
    
    if (is.character(selected)) {
      
      vars <- selected
      
    } else {
      
      if (!outcome_name %in% names(selected)) {
        
        warning(
          paste0(
            "No selected variables specified for outcome '",
            outcome_name,
            "'."
          )
        )
        
        return(
          character(0)
        )
      }
      
      vars <- selected[[outcome_name]]
    }
    
    missing_vars <- setdiff(
      vars,
      names(data)
    )
    
    if (length(missing_vars) > 0) {
      
      warning(
        paste0(
          "Selected variables missing for ",
          outcome_name,
          ": ",
          paste(
            missing_vars,
            collapse = ", "
          )
        )
      )
    }
    
    intersect(
      vars,
      names(data)
    )
  }
  
  
  # ============================================================
  # HELPER: formula
  # ============================================================
  
  build_formula <- function(
    info,
    predictors
  ) {
    
    predictors <- unique(
      predictors
    )
    
    predictors <- intersect(
      predictors,
      names(data)
    )
    
    if (length(predictors) == 0) {
      stop(
        "No valid predictors supplied."
      )
    }
    
    rhs <- paste(
      sprintf(
        "`%s`",
        predictors
      ),
      collapse = " + "
    )
    
    
    if (info$type == "survival") {
      
      lhs <- paste0(
        "survival::Surv(`",
        info$time,
        "`, `",
        info$event,
        "`)"
      )
      
    } else {
      
      lhs <- paste0(
        "`",
        info$variable,
        "`"
      )
    }
    
    
    stats::as.formula(
      paste(
        lhs,
        "~",
        rhs
      )
    )
  }
  
  
  # ============================================================
  # HELPER:
  # map coefficient names to original variables
  #
  # Example:
  #
  # SEXMale      -> SEX
  # RACEBlack    -> RACE
  # ============================================================
  
  get_term_variable_map <- function(
    model
  ) {
    
    mm <- tryCatch(
      stats::model.matrix(model),
      error = function(e) NULL
    )
    
    if (is.null(mm)) {
      return(
        data.frame(
          term = character(0),
          variable = character(0)
        )
      )
    }
    
    assignment <- attr(
      mm,
      "assign"
    )
    
    mm_names <- colnames(
      mm
    )
    
    model_terms <- tryCatch(
      stats::terms(model),
      error = function(e) NULL
    )
    
    if (is.null(model_terms)) {
      return(
        data.frame(
          term = mm_names,
          variable = NA_character_
        )
      )
    }
    
    term_labels <- attr(
      model_terms,
      "term.labels"
    )
    
    variable <- rep(
      NA_character_,
      length(mm_names)
    )
    
    for (i in seq_along(mm_names)) {
      
      a <- assignment[i]
      
      if (
        !is.na(a) &&
        a == 0
      ) {
        
        variable[i] <- "(Intercept)"
        
      } else if (
        !is.na(a) &&
        a > 0 &&
        a <= length(term_labels)
      ) {
        
        variable[i] <- clean_name(
          term_labels[a]
        )
      }
    }
    
    data.frame(
      term = mm_names,
      variable = variable,
      stringsAsFactors = FALSE
    )
  }

  # ============================================================
  # HELPER: cluster-robust covariance without sandwich package
  # HC1 finite-sample correction; independent clusters assumed.
  # ============================================================

  cluster_vcov <- function(model, cluster_values) {
    X <- stats::model.matrix(model)
    beta <- stats::coef(model)
    keep <- !is.na(beta)
    X <- X[, names(beta)[keep], drop = FALSE]
    n <- nrow(X)
    p <- ncol(X)

    # Match the fitted rows after automatic NA omission.
    model_rows <- rownames(stats::model.frame(model))
    source_rows <- names(cluster_values)
    if (is.null(source_rows)) {
      stop("Cluster values must be named by source data row names.")
    }
    idx <- match(model_rows, source_rows)
    if (anyNA(idx)) stop("Could not align cluster IDs to model rows.")
    groups <- cluster_values[idx]
    if (anyNA(groups)) stop("Cluster IDs cannot be missing in fitted rows.")
    G <- length(unique(groups))
    if (G < 2L) stop("Cluster-robust inference requires at least two clusters.")
    if (n <= p) stop("Insufficient residual degrees of freedom for cluster covariance.")

    if (inherits(model, "glm")) {
      mu <- model$fitted.values
      family <- model$family
      prior_w <- model$prior.weights
      score_factor <- prior_w * (model$y - mu) *
        family$mu.eta(model$linear.predictors) / family$variance(mu)
      info_w <- model$weights  # final IRLS working weights
    } else if (inherits(model, "lm")) {
      prior_w <- model$weights
      if (is.null(prior_w)) prior_w <- rep(1, n)
      score_factor <- prior_w * model$residuals
      info_w <- prior_w
    } else {
      stop("cluster_vcov supports only lm and glm objects.")
    }

    if (any(!is.finite(score_factor)) || any(!is.finite(info_w))) {
      stop("Non-finite model scores or weights prevent robust covariance estimation.")
    }
    scores <- X * as.numeric(score_factor)
    cluster_scores <- rowsum(scores, group = groups, reorder = FALSE)
    bread <- solve(crossprod(X, X * as.numeric(info_w)))
    correction <- (G / (G - 1)) * ((n - 1) / (n - p))
    V <- correction * bread %*% crossprod(cluster_scores) %*% bread
    dimnames(V) <- list(colnames(X), colnames(X))
    attr(V, "n_clusters") <- G
    V
  }
  
  
  # ============================================================
  # HELPER: fit model
  # ============================================================
  
  fit_model <- function(
    info,
    predictors,
    analysis_data
  ) {

    data <- analysis_data
    
    predictors <- unique(
      predictors
    )
    
    formula <- build_formula(
      info,
      predictors
    )
    
    
    # Logistic
    if (info$type == "logistic") {
      
      model <- stats::glm(
        formula,
        data = data,
        family = stats::binomial(
          link = "logit"
        )
      )
      
      exponentiate <- TRUE
    }
    
    
    # Poisson
    else if (info$type == "poisson") {
      
      model <- stats::glm(
        formula,
        data = data,
        family = stats::poisson(
          link = "log"
        )
      )
      
      exponentiate <- TRUE
    }
    
    
    # Linear
    else if (info$type == "linear") {
      
      model <- stats::lm(
        formula,
        data = data
      )
      
      exponentiate <- FALSE
    }
    
    
    # Survival / competing risks
    else if (info$type == "survival") {

      if (info$method == "cox") {
        if (is.null(cluster)) {
          model <- survival::coxph(
            formula, data = data, x = TRUE, model = TRUE
          )
        } else {
          model <- survival::coxph(
            formula, data = data, cluster = data[[cluster]],
            robust = TRUE, x = TRUE, model = TRUE
          )
        }
        exponentiate <- TRUE
      }

      else if (info$method == "cause_specific") {
        analysis_data <- data
        analysis_data$.statspipe_cs_event <- as.integer(
          analysis_data[[info$event]] == info$event_of_interest
        )

        rhs <- paste(sprintf("`%s`", predictors), collapse = " + ")
        cs_formula <- stats::as.formula(
          paste0(
            "survival::Surv(`", info$time,
            "`, .statspipe_cs_event) ~ ", rhs
          )
        )

        if (is.null(cluster)) {
          model <- survival::coxph(
            cs_formula, data = analysis_data, x = TRUE, model = TRUE
          )
        } else {
          model <- survival::coxph(
            cs_formula, data = analysis_data,
            cluster = analysis_data[[cluster]], robust = TRUE,
            x = TRUE, model = TRUE
          )
        }
        exponentiate <- TRUE
      }

      else if (info$method == "finegray") {
        required_vars <- unique(c(info$time, info$event, predictors, cluster))
        fg_source <- data[
          stats::complete.cases(data[, required_vars, drop = FALSE]),
          required_vars, drop = FALSE
        ]
        original_n <- nrow(fg_source)

        # Preserve subject IDs through expansion for subject-level robust SEs.
        fg_source$.statspipe_fg_id <- seq_len(original_n)
        # Fine-Gray requires a multi-state status factor with censoring first.
        fg_events <- as.character(fg_source[[info$event]])
        fg_source$.statspipe_fg_status <- factor(
          fg_events,
          levels = unique(c(as.character(info$censor_value), fg_events))
        )
        target_level <- as.character(info$event_of_interest)

        if (!target_level %in% levels(fg_source$.statspipe_fg_status)) {
          stop(paste0(
            "event_of_interest=", info$event_of_interest,
            " is not present in outcome '", info$name, "'."
          ))
        }

        # Include cluster in finegray() input so it survives row expansion,
        # but do not include it as a modeled predictor.
        fg_rhs_terms <- c(sprintf("`%s`", predictors))
        fg_rhs_terms <- c(fg_rhs_terms, ".statspipe_fg_id")
        if (!is.null(cluster)) {
          fg_rhs_terms <- c(fg_rhs_terms, sprintf("`%s`", cluster))
        }

        fg_formula <- stats::as.formula(
          paste0(
            "survival::Surv(`", info$time,
            "`, .statspipe_fg_status) ~ ",
            paste(fg_rhs_terms, collapse = " + ")
          )
        )

        fg_data <- survival::finegray(
          formula = fg_formula,
          data = fg_source,
          etype = target_level
        )

        fg_cox_formula <- stats::as.formula(
          paste(
            "survival::Surv(fgstart, fgstop, fgstatus) ~",
            paste(sprintf("`%s`", predictors), collapse = " + ")
          )
        )

        if (is.null(cluster)) {
          model <- survival::coxph(
            fg_cox_formula, data = fg_data, weights = fgwt,
            cluster = fg_data$.statspipe_fg_id, robust = TRUE,
            x = TRUE, model = TRUE
          )
        } else {
          model <- survival::coxph(
            fg_cox_formula, data = fg_data, weights = fgwt,
            cluster = fg_data[[cluster]], robust = TRUE,
            x = TRUE, model = TRUE
          )
        }

        attr(model, "statspipe_original_n") <- original_n
        exponentiate <- TRUE
      }
    }
    
    
    # Tidy coefficients
    if (!is.null(cluster) && info$type %in% c("logistic", "poisson", "linear")) {
      cluster_values <- data[[cluster]]
      names(cluster_values) <- rownames(data)
      robust_vcov <- cluster_vcov(model, cluster_values)
      beta <- stats::coef(model)
      se <- rep(NA_real_, length(beta))
      names(se) <- names(beta)
      se[names(robust_vcov[, 1])] <- sqrt(diag(robust_vcov))
      cluster_df <- attr(robust_vcov, "n_clusters") - 1L
      zcrit <- stats::qt(1 - (1 - conf_level) / 2, df = cluster_df)
      est <- beta
      lo <- beta - zcrit * se
      hi <- beta + zcrit * se
      pval <- 2 * stats::pt(abs(beta / se), df = cluster_df, lower.tail = FALSE)

      if (isTRUE(exponentiate)) {
        est <- exp(est); lo <- exp(lo); hi <- exp(hi)
      }

      tidy_result <- data.frame(
        term = names(beta), estimate = unname(est),
        std.error = unname(se), statistic = unname(beta / se),
        p.value = unname(pval), conf.low = unname(lo),
        conf.high = unname(hi), stringsAsFactors = FALSE
      )
    } else {
      tidy_result <- broom::tidy(
        model, conf.int = TRUE, conf.level = conf_level,
        exponentiate = exponentiate
      )
    }
    
    
    # Map coefficients back to variables
    term_map <- get_term_variable_map(
      model
    )
    
    tidy_result <- dplyr::left_join(
      tidy_result,
      term_map,
      by = "term"
    )
    
    
    # Fallback mapping
    missing_mapping <- is.na(
      tidy_result$variable
    )
    
    if (any(missing_mapping)) {
      
      for (
        j in which(missing_mapping)
      ) {
        
        current_term <- clean_name(
          tidy_result$term[j]
        )
        
        matches <- predictors[
          startsWith(
            current_term,
            predictors
          )
        ]
        
        if (length(matches) > 0) {
          
          matches <- matches[
            order(
              nchar(matches),
              decreasing = TRUE
            )
          ]
          
          tidy_result$variable[j] <-
            matches[1]
        }
      }
    }
    
    
    tidy_result <- tidy_result |>
      dplyr::mutate(
        outcome = info$name,
        model_type = info$type,
        measure = info$measure,
        n = if (
          info$type == "survival" &&
          !is.null(info$method) &&
          info$method == "finegray" &&
          !is.null(attr(model, "statspipe_original_n"))
        ) {
          attr(model, "statspipe_original_n")
        } else {
          if (info$type == "survival") model$n else stats::nobs(model)
        }
      )
    
    
    list(
      model = model,
      tidy = tidy_result
    )
  }
  
  
  # ============================================================
  # HELPER: format regression table
  # ============================================================
  
  format_regression_table <- function(
    x
  ) {
    
    x |>
      dplyr::mutate(
        
        estimate_ci = dplyr::case_when(
          
          measure %in%
            c(
              "OR",
              "IRR",
              "HR",
              "sHR"
            ) ~
            sprintf(
              "%.2f (%.2f, %.2f)",
              estimate,
              conf.low,
              conf.high
            ),
          
          TRUE ~
            sprintf(
              "%.3f (%.3f, %.3f)",
              estimate,
              conf.low,
              conf.high
            )
        ),
        
        p_value = dplyr::case_when(
          
          is.na(p.value) ~
            NA_character_,
          
          p.value < 0.001 ~
            "<0.001",
          
          TRUE ~
            sprintf(
              "%.3f",
              p.value
            )
        )
      ) |>
      
      dplyr::select(
        outcome,
        variable,
        term,
        measure,
        estimate,
        conf.low,
        conf.high,
        estimate_ci,
        p.value,
        p_value,
        n,
        dplyr::everything()
      )
  }
  
  
  # ============================================================
  # 2. TABLE 1
  # ============================================================
  
  table1 <- list()
  log_progress("Generating descriptive tables")
  
  
  # Overall source cohort (for backward compatibility)
  table1$overall <-
    data |>
    dplyr::select(
      dplyr::all_of(
        table1_vars
      )
    ) |>
    gtsummary::tbl_summary(
      missing = "ifany"
    ) |>
    gtsummary::add_n()
  
  
  # ------------------------------------------------------------
  # Stratified Table 1; strata can either be outcome or exposure
  # ------------------------------------------------------------

  table1$by_strata <- list()
  if (!is.null(strata) && length(strata) > 0) {
    for (s in strata) {
      
      vars_for_table <- unique(
        c(
          s,
          table1_vars
        )
      )
      
      tbl_data <-
        data |>
        dplyr::select(
          dplyr::all_of(
            vars_for_table
          )
        )
      
      tbl <-
        gtsummary::tbl_summary(
          tbl_data,
          by = !!rlang::sym(s),
          missing = "ifany"
        ) |>
        gtsummary::add_overall() |>
        gtsummary::add_n()
      
      tbl <- tryCatch(
        gtsummary::add_p(
          tbl,
          test.args = gtsummary::all_tests("fisher.test") ~
            list(simulate.p.value = TRUE, B = 10000L)
        ),
        
        error = function(e) {
          warning(
            paste0(
              "Unable to calculate Table 1 p-values for ",
              s,
              ": ",
              e$message
            )
          )
          
          tbl
        }
      )
      
      table1$by_strata[[s]] <- tbl
    }
  }
  
  # ============================================================
  # 3. UNIVARIATE REGRESSION
  # ============================================================
  
  univariate_models <- list()
  univariate_results <- list()
  
  
  for (
    outcome_name in names(outcomes)
  ) {
    log_progress("Univariate models: ", outcome_name, " (", length(covariates), " predictors)")
    
    info <- get_outcome_info(
      outcome_name,
      outcomes[[outcome_name]]
    )
    
    outcome_models <- list()
    outcome_results <- list()
    
    
    for (var in covariates) {
      
      fit <- tryCatch(
        
        fit_model(
          info,
          var,
          cohort_data[[outcome_name]]
        ),
        
        error = function(e) {
          
          warning(
            paste0(
              "Unable to fit univariate model ",
              outcome_name,
              " ~ ",
              var,
              ": ",
              e$message
            )
          )
          
          NULL
        }
      )
      
      
      if (!is.null(fit)) {
        
        outcome_models[[var]] <-
          fit$model
        
        
        result_i <-
          fit$tidy |>
          dplyr::filter(
            term != "(Intercept)"
          )
        
        
        # Explicit assignment protects categorical predictors
        result_i$variable <- var
        
        result_i$analysis <-
          "Univariate"
        
        
        outcome_results[[var]] <-
          result_i
      }
    }
    
    
    univariate_models[[outcome_name]] <-
      outcome_models
    
    
    if (length(outcome_results) > 0) {
      
      univariate_results[[outcome_name]] <-
        dplyr::bind_rows(
          outcome_results
        ) |>
        format_regression_table()
    }
  }
  
  
  # ============================================================
  # 4. MULTIVARIABLE REGRESSION
  # ============================================================
  
  multivariable_models <- list()
  multivariable_results <- list()
  
  
  for (
    outcome_name in names(outcomes)
  ) {
    log_progress("Multivariable model: ", outcome_name)
    
    info <- get_outcome_info(
      outcome_name,
      outcomes[[outcome_name]]
    )
    
    
    selected_vars <-
      get_selected_variables(
        outcome_name
      )
    
    
    if (length(selected_vars) == 0) {
      
      warning(
        paste0(
          "No valid selected variables for ",
          outcome_name,
          ". Multivariable model skipped."
        )
      )
      
      next
    }
    
    
    fit <- tryCatch(
      
      fit_model(
        info,
        selected_vars,
        cohort_data[[outcome_name]]
      ),
      
      error = function(e) {
        
        warning(
          paste0(
            "Unable to fit multivariable model for ",
            outcome_name,
            ": ",
            e$message
          )
        )
        
        NULL
      }
    )
    
    
    if (!is.null(fit)) {
      
      multivariable_models[[outcome_name]] <-
        fit$model
      
      
      multivariable_results[[outcome_name]] <-
        fit$tidy |>
        dplyr::filter(
          term != "(Intercept)"
        ) |>
        dplyr::mutate(
          analysis = "Multivariable"
        ) |>
        format_regression_table()
    }
  }
  
  
  # ============================================================
  # 5. COMBINED RESULTS
  # ============================================================
  
  if (length(univariate_results) > 0) {
    
    univariate_combined <-
      dplyr::bind_rows(
        univariate_results
      )
    
  } else {
    
    univariate_combined <- NULL
  }
  
  
  if (length(multivariable_results) > 0) {
    
    multivariable_combined <-
      dplyr::bind_rows(
        multivariable_results
      )
    
  } else {
    
    multivariable_combined <- NULL
  }
  
  
  # ============================================================
  # 6. FOREST PLOT
  #
  # One outcome per column.
  #
  # OR / IRR / HR only.
  #
  # Annotation:
  #
  # 1.32 (1.15, 1.51)
  # p<0.001
  #
  # ============================================================
  
  if (forest_source == "univariate") {
    
    forest_data <- univariate_combined
    
  } else {
    
    forest_data <- multivariable_combined
  }
  
  
  forest_plot <- NULL
  forest_plots <- list()
  log_progress("Preparing forest plots from ", forest_source, " results")
  
  
  if (!is.null(forest_data)) {
    
    # ----------------------------------------------------------
    # Ratio-based models only
    # ----------------------------------------------------------
    
    forest_data <-
      forest_data |>
      dplyr::filter(
        measure %in%
          c(
            "OR",
            "IRR",
            "HR",
            "sHR"
          ),
        is.finite(estimate),
        is.finite(conf.low),
        is.finite(conf.high),
        estimate > 0,
        conf.low > 0,
        conf.high > 0
      )
    
    
    # ----------------------------------------------------------
    # Optional variable/term filtering
    # ----------------------------------------------------------
    
    if (
      !is.null(forest_include) &&
      nrow(forest_data) > 0
    ) {
      
      forest_data <-
        forest_data |>
        dplyr::filter(
          variable %in% forest_include |
            term %in% forest_include
        )
    }


    # ----------------------------------------------------------
    # Optional exclusion of categorical levels from final plot
    #
    # Example:
    # forest_exclude_levels = "NI"
    #
    # removes terms such as RACE_NI / RACENI from the forest
    # plot only. Model fitting and regression tables are unchanged.
    # ----------------------------------------------------------

    if (
      !is.null(forest_exclude_levels) &&
      length(forest_exclude_levels) > 0 &&
      nrow(forest_data) > 0
    ) {

      exclude_levels <- as.character(forest_exclude_levels)

      forest_data <-
        forest_data |>
        dplyr::mutate(
          .term_clean = clean_name(as.character(term)),
          .variable_clean = clean_name(as.character(variable))
        ) |>
        dplyr::rowwise() |>
        dplyr::filter(
          !any(vapply(
            exclude_levels,
            function(level) {
              term_i <- .term_clean
              var_i <- .variable_clean

              identical(term_i, paste0(var_i, level)) ||
                identical(term_i, paste0(var_i, "_", level)) ||
                identical(term_i, level)
            },
            logical(1)
          ))
        ) |>
        dplyr::ungroup() |>
        dplyr::select(-.term_clean, -.variable_clean)
    }
    
    
    if (nrow(forest_data) > 0) {

      n_outcomes_plot <- dplyr::n_distinct(forest_data$outcome)
      if (!is.null(forest_ncol) && !is.null(forest_nrow) &&
          forest_ncol * forest_nrow < n_outcomes_plot) {
        stop("forest_ncol * forest_nrow is smaller than the number of plotted outcomes.")
      }
      plot_cols <- if (!is.null(forest_ncol)) forest_ncol else if (
        !is.null(forest_nrow)) ceiling(n_outcomes_plot / forest_nrow) else
        grDevices::n2mfrow(n_outcomes_plot)[1L]
      plot_cols <- min(plot_cols, n_outcomes_plot)
      plot_rows <- ceiling(n_outcomes_plot / plot_cols)
      incomplete_row <- plot_rows > 1L && n_outcomes_plot %% plot_cols != 0L
      if (incomplete_row && !requireNamespace("patchwork", quietly = TRUE)) {
        stop("Please install 'patchwork' to align incomplete forest plot rows.")
      }
      
      # --------------------------------------------------------
      # Preserve outcome order
      # --------------------------------------------------------
      
      forest_data$outcome <-
        factor(
          forest_data$outcome,
          levels = names(outcomes)
        )
      
      
      # --------------------------------------------------------
      # Outcome headers
      #
      # Mortality
      # OR (95% CI)
      # --------------------------------------------------------
      
      outcome_labels <-
        forest_data |>
        dplyr::group_by(
          outcome
        ) |>
        dplyr::summarise(
          measure_label =
            dplyr::first(
              measure
            ),
          .groups = "drop"
        ) |>
        dplyr::mutate(
          outcome_label = paste0(
            as.character(
              outcome
            ),
            "\n",
            measure_label,
            " (", ci_label, " CI)"
          )
        )
      
      
      forest_data <-
        dplyr::left_join(
          forest_data,
          outcome_labels,
          by = "outcome"
        )
      
      
      label_levels <-
        outcome_labels$outcome_label[
          order(
            match(
              outcome_labels$outcome,
              names(outcomes)
            )
          )
        ]
      
      
      forest_data$outcome_label <-
        factor(
          forest_data$outcome_label,
          levels = label_levels
        )
      
      
      # --------------------------------------------------------
      # Determine variable ordering
      # --------------------------------------------------------
      
      if (forest_source == "univariate") {
        
        desired_variable_order <-
          covariates
        
      } else {
        
        if (is.character(selected)) {
          
          desired_variable_order <-
            selected
          
        } else {
          
          desired_variable_order <-
            unique(
              unlist(
                selected,
                use.names = FALSE
              )
            )
        }
      }
      
      
      # --------------------------------------------------------
      # Order factor levels together under parent variable
      # --------------------------------------------------------
      
      term_order <-
        forest_data |>
        dplyr::mutate(
          variable_order = match(
            variable,
            desired_variable_order
          )
        ) |>
        dplyr::arrange(
          variable_order,
          variable
        ) |>
        dplyr::distinct(
          variable,
          term
        ) |>
        dplyr::pull(
          term
        )
      
      
      forest_data$term <-
        factor(
          forest_data$term,
          levels = rev(
            unique(
              term_order
            )
          )
        )

      forest_term_labels <- setNames(
        levels(forest_data$term), levels(forest_data$term)
      )
      for (term_name in names(forest_term_labels)) {
        variable_name <- as.character(forest_data$variable[
          match(term_name, as.character(forest_data$term))
        ])
        if (term_name %in% names(forest_labels)) {
          forest_term_labels[term_name] <- forest_labels[[term_name]]
        } else if (!is.na(variable_name) && variable_name %in% names(forest_labels)) {
          level <- clean_name(term_name)
          variable_clean <- clean_name(variable_name)
          if (startsWith(level, variable_clean)) {
            level <- substring(level, nchar(variable_clean) + 1L)
            level <- sub("^_", "", level)
          }
          forest_term_labels[term_name] <- if (nzchar(level)) {
            paste0(forest_labels[[variable_name]], ": ", level)
          } else {
            forest_labels[[variable_name]]
          }
        }
      }
      forest_data$term_label <- unname(
        forest_term_labels[as.character(forest_data$term)]
      )
      
      
      # --------------------------------------------------------
      # Annotation text
      # --------------------------------------------------------
      
      forest_data <-
        forest_data |>
        dplyr::mutate(
          
          estimate_label =
            sprintf(
              "%.2f (%.2f, %.2f)",
              estimate,
              conf.low,
              conf.high
            ),
          
          p_label =
            dplyr::case_when(
              
              is.na(p.value) ~
                "",
              
              p.value < 0.001 ~
                "p<0.001",
              
              TRUE ~
                paste0(
                  "p=",
                  sprintf(
                    "%.3f",
                    p.value
                  )
                )
            ),
          
          annotation =
            paste0(
              estimate_label,
              "\n",
              p_label
            )
        )
      
      
      # ========================================================
      # FOREST PLOT
      # ========================================================
      
      forest_plot <-
        ggplot2::ggplot(
          forest_data,
          ggplot2::aes(
            x = estimate,
            y = term
          )
        ) +
        
        # Null line
        ggplot2::geom_vline(
          xintercept = 1,
          linetype = 2,
          linewidth = 0.5
        ) +
        
        # CI
        ggplot2::geom_errorbar(
          ggplot2::aes(
            xmin = conf.low,
            xmax = conf.high
          ),
          orientation = "y",
          width = 0.15,
          linewidth = 0.6
        ) +
        
        # Estimate
        ggplot2::geom_point(
          size = 2.5
        )
      
      
      # --------------------------------------------------------
      # Annotation
      # --------------------------------------------------------
      
      if (isTRUE(forest_annotate)) {
        
        forest_plot <-
          forest_plot +
          
          ggplot2::geom_text(
            ggplot2::aes(
              x = conf.high,
              label = annotation
            ),
            hjust = -0.08,
            size = 3,
            fontface = "bold",
            lineheight = 0.9
          )
      }
      
      
      # --------------------------------------------------------
      # Complete plot
      # --------------------------------------------------------
      
      forest_plot <-
        forest_plot +
        
        ggplot2::scale_x_log10(
          expand = ggplot2::expansion(
            mult = c(
              0.08,
              ifelse(
                forest_annotate,
                0.55,
                0.10
              )
            )
          )
        ) +
        
        ggplot2::facet_wrap(
          ggplot2::vars(outcome_label),
          ncol = plot_cols,
          nrow = plot_rows,
          scales = "free_x"
        ) +

        ggplot2::scale_y_discrete(labels = forest_term_labels, drop = FALSE) +
        
        ggplot2::labs(
          x = "Effect estimate",
          y = NULL,
          title = paste0(
            tools::toTitleCase(
              forest_source
            ),
            " regression results"
          ),
          subtitle = forest_subtitle
        ) +
        
        ggplot2::coord_cartesian(
          clip = "off"
        ) +
        
        ggplot2::theme_bw() +
        
        ggplot2::theme(
          text = ggplot2::element_text(face = "bold"),
          
          panel.grid.minor =
            ggplot2::element_blank(),
          
          panel.grid.major.y =
            ggplot2::element_blank(),
          
          strip.background =
            ggplot2::element_blank(),
          
          strip.text =
            ggplot2::element_text(
              face = "bold",
              size = 11
            ),
          
          axis.text.y =
            ggplot2::element_text(
              size = 9
            ),
          
          axis.text.x =
            ggplot2::element_text(
              size = 9
            ),
          
          axis.title.x =
            ggplot2::element_text(
              margin =
                ggplot2::margin(
                  t = 8
                )
            ),
          
          panel.spacing.x =
            grid::unit(
              1.5,
              "lines"
            ),
          
          plot.margin =
            ggplot2::margin(
              t = 10,
              r = 35,
              b = 10,
              l = 10
            )
        )

      # Reuse the completed plot style with each outcome's data.
      for (outcome_name in unique(as.character(forest_data$outcome))) {
        individual <- forest_plot
        individual$data <- forest_data[
          as.character(forest_data$outcome) == outcome_name, , drop = FALSE
        ]
        individual$data$term <- droplevels(individual$data$term)
        individual <- individual +
          ggplot2::facet_wrap(ggplot2::vars(outcome_label), ncol = 1) +
          ggplot2::labs(title = paste0(outcome_name, " — ", forest_source, " regression"))
        forest_plots[[outcome_name]] <- individual
      }

      if (incomplete_row) {
        row_indices <- split(seq_along(label_levels),
                             ceiling(seq_along(label_levels) / plot_cols))
        row_plots <- lapply(row_indices, function(indices) {
          row_plot <- forest_plot
          row_plot$data <- forest_data[
            forest_data$outcome_label %in% label_levels[indices], , drop = FALSE
          ]
          row_plot +
            ggplot2::facet_wrap(
              ggplot2::vars(outcome_label), ncol = length(indices), scales = "free_x"
            ) +
            ggplot2::labs(title = NULL, subtitle = NULL)
        })
        forest_plot <- patchwork::wrap_plots(row_plots, ncol = 1, heights = 1) +
          patchwork::plot_annotation(
            title = paste0(tools::toTitleCase(forest_source), " regression results"),
            subtitle = forest_subtitle,
            theme = ggplot2::theme(
              text = ggplot2::element_text(face = "bold")
            )
          )
      }
      log_progress("Created combined forest plot and ", length(forest_plots),
                   " individual outcome plots")
    }
  }
  
  
  # ============================================================
  # 7. RESULTS OBJECT
  # ============================================================
  
  results <- list(

    settings = list(
      conf_level = conf_level,
      cluster = cluster,
      cohort_filters = cohort_filters,
      cohort_counts = cohort_counts,
      forest_ncol = forest_ncol,
      forest_nrow = forest_nrow,
      forest_labels = forest_labels,
      verbose = verbose
    ),

    table1 = table1,
    
    univariate = list(
      tables = univariate_results,
      combined = univariate_combined,
      models = univariate_models
    ),
    
    multivariable = list(
      tables = multivariable_results,
      combined = multivariable_combined,
      models = multivariable_models
    ),
    
    forest = list(
      data = forest_data,
      plot = forest_plot,
      plots = forest_plots
    ),
    
    saved_files = NULL
  )
  
  
  # ============================================================
  # 8. SAVE
  # ============================================================
  
  if (isTRUE(save)) {
    log_progress("Saving analysis outputs")
    
    if (is.null(output_dir)) {
      stop(
        "'output_dir' must be specified when save = TRUE."
      )
    }
    
    
    # ----------------------------------------------------------
    # Main output folder
    # ----------------------------------------------------------
    
    if (!dir.exists(output_dir)) {
      
      dir.create(
        output_dir,
        recursive = TRUE,
        showWarnings = FALSE
      )
    }
    
    
    if (!dir.exists(output_dir)) {
      
      stop(
        paste0(
          "Unable to create output directory: ",
          output_dir
        )
      )
    }
    
    
    # ----------------------------------------------------------
    # Subfolders
    # ----------------------------------------------------------
    
    dirs <- list(
      
      table1 =
        file.path(
          output_dir,
          "table1"
        ),
      
      univariate =
        file.path(
          output_dir,
          "univariate"
        ),
      
      multivariable =
        file.path(
          output_dir,
          "multivariable"
        ),
      
      forest =
        file.path(
          output_dir,
          "forest"
        ),
      
      models =
        file.path(
          output_dir,
          "models"
        )
    )
    
    
    for (folder in dirs) {
      
      if (!dir.exists(folder)) {
        
        dir.create(
          folder,
          recursive = TRUE,
          showWarnings = FALSE
        )
      }
      if (!dir.exists(folder)) {
        stop("Unable to create output directory: ", folder)
      }
    }
    
    
    # ==========================================================
    # Save Table 1
    # ==========================================================
    
    table1_to_save <- list(overall = table1$overall)
    for (s in names(table1$by_strata)) {
      table1_to_save[[paste0("by_strata_", s)]] <- table1$by_strata[[s]]
    }

    for (nm in names(table1_to_save)) {
      
      # RDS
      saveRDS(
        table1_to_save[[nm]],
        file = file.path(
          dirs$table1,
          paste0(
            prefix,
            "_table1_",
            nm,
            ".rds"
          )
        )
      )
      
      
      # CSV
      table1_csv <-
        tryCatch(
          
          gtsummary::as_tibble(
            table1_to_save[[nm]],
            col_labels = FALSE
          ),
          
          error = function(e) {
            warning("Unable to format Table 1 '", nm, "' for CSV: ",
                    conditionMessage(e), ". Exporting the table body instead.")
            table1_to_save[[nm]]$table_body
          }
        )
      
      
      utils::write.csv(
        table1_csv,
        file = file.path(
          dirs$table1,
          paste0(
            prefix,
            "_table1_",
            nm,
            ".csv"
          )
        ),
        row.names = FALSE,
        na = ""
      )
    }
    
    
    # ==========================================================
    # Save univariate results
    # ==========================================================
    
    if (length(univariate_results) > 0) {
      
      for (
        nm in names(
          univariate_results
        )
      ) {
        
        utils::write.csv(
          univariate_results[[nm]],
          file = file.path(
            dirs$univariate,
            paste0(
              prefix,
              "_univariate_",
              nm,
              ".csv"
            )
          ),
          row.names = FALSE,
          na = ""
        )
      }
      
      
      if (!is.null(univariate_combined)) {
        
        utils::write.csv(
          univariate_combined,
          file = file.path(
            dirs$univariate,
            paste0(
              prefix,
              "_univariate_all_outcomes.csv"
            )
          ),
          row.names = FALSE,
          na = ""
        )
      }
    }
    
    
    # ==========================================================
    # Save multivariable results
    # ==========================================================
    
    if (length(multivariable_results) > 0) {
      
      for (
        nm in names(
          multivariable_results
        )
      ) {
        
        utils::write.csv(
          multivariable_results[[nm]],
          file = file.path(
            dirs$multivariable,
            paste0(
              prefix,
              "_multivariable_",
              nm,
              ".csv"
            )
          ),
          row.names = FALSE,
          na = ""
        )
      }
      
      
      if (!is.null(multivariable_combined)) {
        
        utils::write.csv(
          multivariable_combined,
          file = file.path(
            dirs$multivariable,
            paste0(
              prefix,
              "_multivariable_all_outcomes.csv"
            )
          ),
          row.names = FALSE,
          na = ""
        )
      }
    }
    
    
    # ==========================================================
    # Save model objects
    # ==========================================================
    
    saveRDS(
      univariate_models,
      file = file.path(
        dirs$models,
        paste0(
          prefix,
          "_univariate_models.rds"
        )
      )
    )
    
    
    saveRDS(
      multivariable_models,
      file = file.path(
        dirs$models,
        paste0(
          prefix,
          "_multivariable_models.rds"
        )
      )
    )
    
    
    # ==========================================================
    # Save forest data
    # ==========================================================
    
    if (
      !is.null(forest_data) &&
      nrow(forest_data) > 0
    ) {
      
      utils::write.csv(
        forest_data,
        file = file.path(
          dirs$forest,
          paste0(
            prefix,
            "_forest_data.csv"
          )
        ),
        row.names = FALSE,
        na = ""
      )
    }
    
    
    # ==========================================================
    # Save forest plot
    # ==========================================================
    
    if (!is.null(forest_plot)) {
      
      n_outcomes_plot <-
        dplyr::n_distinct(
          forest_data$outcome
        )
      
      
      n_terms_plot <-
        dplyr::n_distinct(
          forest_data$term
        )

      
      # More horizontal room when annotations are included
      plot_width <-
        max(
          8,
          ifelse(
            forest_annotate,
            4.3,
            3.2
          ) *
            plot_cols +
            2
        )
      
      
      plot_height <-
        max(
          5,
          (0.45 * n_terms_plot + 2) * plot_rows
        )
      
      
      # PDF
      ggplot2::ggsave(
        filename = file.path(
          dirs$forest,
          paste0(
            prefix,
            "_forest_plot.pdf"
          )
        ),
        plot = forest_plot,
        width = plot_width,
        height = plot_height,
        limitsize = FALSE
      )
      
      
      # PNG
      ggplot2::ggsave(
        filename = file.path(
          dirs$forest,
          paste0(
            prefix,
            "_forest_plot.png"
          )
        ),
        plot = forest_plot,
        width = plot_width,
        height = plot_height,
        dpi = 300,
        limitsize = FALSE
      )

      # Save a separate PDF and PNG for every plotted outcome.
      for (outcome_name in names(forest_plots)) {
        safe_name <- gsub("[^A-Za-z0-9_-]+", "_", outcome_name)
        individual_data <- forest_data[
          as.character(forest_data$outcome) == outcome_name, , drop = FALSE
        ]
        individual_width <- if (isTRUE(forest_annotate)) 8.5 else 7
        individual_height <- max(5, 0.45 *
          dplyr::n_distinct(individual_data$term) + 2)
        for (ext in c("pdf", "png")) {
          ggplot2::ggsave(
            filename = file.path(dirs$forest, paste0(
              prefix, "_forest_", safe_name, ".", ext
            )),
            plot = forest_plots[[outcome_name]],
            width = individual_width,
            height = individual_height,
            dpi = 300,
            limitsize = FALSE
          )
        }
        log_progress("Saved individual forest plot: ", outcome_name)
      }
    }
    
    
    # ==========================================================
    # Saved-file information
    # ==========================================================
    
    complete_results_file <-
      file.path(
        output_dir,
        paste0(
          prefix,
          "_complete_results.rds"
        )
      )
    
    
    results$saved_files <- list(
      
      output_dir =
        normalizePath(
          output_dir,
          winslash = "/",
          mustWork = FALSE
        ),
      
      complete_results =
        complete_results_file,
      
      table1_dir =
        dirs$table1,
      
      univariate_dir =
        dirs$univariate,
      
      multivariable_dir =
        dirs$multivariable,
      
      forest_dir =
        dirs$forest,
      
      models_dir =
        dirs$models
    )
    
    
    # Save complete results object
    saveRDS(
      results,
      complete_results_file
    )
    
    
    message(
      "Analysis completed and saved to: ",
      normalizePath(
        output_dir,
        winslash = "/",
        mustWork = FALSE
      )
    )
  }

  log_progress("Completed analysis")
  
  
  # ============================================================
  # 9. RETURN
  # ============================================================
  
  return(
    results
  )
}
