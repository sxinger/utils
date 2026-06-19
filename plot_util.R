####################################
## utility functions for plotting ##
####################################

## generate transformation value and formula for overlayed plots with dual-axis
## compatiable with dplyr and ggplot
convert_scale<-function(from,to){
  # convert value
  range_rt<-(max(to) - min(to))/(max(from) - min(from))
  new_from<-((from - min(from))*range_rt) + min(to)
  
  # generate formula to recover axis (revert the new_from formula)
  revert_formula<-c(paste0("~ (. - ",min(to),") * ",1/range_rt),
                    rep(NA,length(from)-1))
  
  # results
  return(list(val=new_from,
              formula=revert_formula))
}

plot_dualy<-function(
  df,   # data.frame
  x,    # x-axis
  y,    # primary y-axis - bar plot
  y_sec # secondary y-axis - line plot
){
  # require(tidyverse)
  # covert secondary y to the same scale of primary y
  df_cp<-df %>%
      mutate(new_col = convert_scale(y_sec,y)[["val"]],
             axis_formula = convert_scale(y_sec,y)[["formula"]])

  # ggplot object
  plt<-ggplot(df_cp,aes(x=x)) +
    geom_col(aes(y=y), size = 1, color = "darkblue", fill = "white")+
    geom_line(aes(y=new_col), size = 1.5, color="red", group = 1)+
    scale_y_continuous(sec.axis = sec_axis(as.formula(df_cp$axis_formula[1]), name = "sec_axis"))

  # compatible with ggplot syntax
  return(plt)
}

## forestplot
forestplot.HR <- function (
  df,  # long table
  x_idx1="vari", # 1st layer index
  x_idx2="vari_cat", # 2nd layer index
  y_idx="endpt", # 1st layer y index
  est="est", # estimates
  lower="lower", # 95% CI lower bound
  upper="upper", # 95% CI upper bound
  pval="pval", # p value
  plt_par = list(), # other plotting parameters passed in forest function
  ny = 1, # number of y groups
  idx_display = "Variable",
  lbl_pval = TRUE, # whether to use p-value for the label column
  tm = forest_theme(
    arrow_type = "closed",
    arrow_label_just = "end"
  ) # theme parameters for forest function
){
  # require(tidyverse,grid,forestploter)
  # https://cran.r-project.org/web/packages/forestploter/vignettes/forestploter-intro.html

  plt_par$ci_column=Filter(Negate(is.null),list(plt_par$ci_column,2*seq_len(ny)))[[1]]
  plt_par$ref_line = Filter(Negate(is.null),list(plt_par$ref_line,rep(1, ny)))[[1]]
  plt_par$xlab = Filter(Negate(is.null),list(plt_par$xlab,rep("HR", ny)))[[1]]
  plt_par$nudge_y = c(plt_par$nudge_y,0.2)[1]
  plt_par$xlim = Filter(Negate(is.null),list(plt_par$xlim,rep(list(c(0, 1.5)),ny)))[[1]]
  plt_par$ticks_at = Filter(Negate(is.null),list(plt_par$ticks_at,rep(list(c(0.1, 0.5, 1, 1.5)),ny)))[[1]]

  # optional
  plt_par$vert_line = plt_par$vert_line
  plt_par$arrow_lab = plt_par$arrow_lab
  plt_par$x_trans = plt_par$x_trans

  # change to internal names for easy reference
  nm_map<-data.frame(
    ext_nm=c(x_idx1,x_idx2,y_idx,est,lower,upper,pval)) %>%
    mutate(int_nm=c("x_idx1","x_idx2","y_idx","est","lower","upper","pval"))
  plt_df<-df %>%
    select(all_of(nm_map$ext_nm)) %>%
    rename_at(vars(nm_map$ext_nm), ~ nm_map$int_nm)

  # add an empty column for HR plots and a label column
  plt_df %<>%
    # CI with asterisks
    mutate(
      pvalstar = case_when(pval > 0.1 ~ "",
                          pval <= 0.1 & pval > 0.05 ~ "*",
                          pval <= 0.05 & pval > 0.01 ~ "**",
                          pval <= 0.01 & pval > 0.001 ~ "***",
                          TRUE ~ "****"),
      `..` = paste(rep(" ", 20), collapse = " "),
      `HR (95% CI)` = sprintf("%.2f (%.2f-%.2f) %s",est, lower, upper, pvalstar)
    ) %>%
    # P-values
    mutate(
      pvalrd = case_when(
        pval >= 0.01 ~ as.character(round(pval,2)),
        pval < 0.0005 ~ '0.000',
        TRUE ~ as.character(round(pval,3))
      ),
      `..` = paste(rep(" ", 20), collapse = " "),
      `P-value` = pvalrd
    )

  # pivot wide
  y_grp<-unique(plt_df$y_idx)
  n_grp<-length(y_grp)
  if(lbl_pval){
    plt_df %<>%
      pivot_wider(
        id_cols = c("x_idx1","x_idx2"),
        names_from = "y_idx",
        values_from = c(`..`,est,lower,upper,`P-value`),
        names_glue = "{y_idx}.{.value}"
      )
  }else{
    plt_df %<>%
      pivot_wider(
        id_cols = c("x_idx1","x_idx2"),
        names_from = "y_idx",
        values_from = c(`..`,est,lower,upper,`HR (95% CI)`),
        names_glue = "{y_idx}.{.value}"
      )
  }

  # add header rows
  if(lbl_pval){
    tidy_col<-c("x_idx2",tidyr::expand_grid(!!!list(b=y_grp,a=c("...",".P-value"))) |> pmap_chr(paste0))
  }else{
    tidy_col<-c("x_idx2",tidyr::expand_grid(!!!list(b=y_grp,a=c("...",".HR (95% CI)"))) |> pmap_chr(paste0))
  }
  
  plt_header<-plt_df[FALSE,tidy_col] %>%
    bind_rows(data.frame(
      x_idx1 = unique(plt_df$x_idx1),
      x_idx2 = unique(plt_df$x_idx1),
      idx=0
    )) %>%
    replace(is.na(.), " ")

  plt_df %<>%
    mutate_at(tidy_col, ~replace_na(.," ")) %>%
    mutate(x_idx2 = paste0("  ",x_idx2)) %>%
    bind_rows(plt_header) %>%
    arrange(x_idx1,idx) %>% select(-idx)

  # collect list for HR segment plots
  est_lst<-list()
  lower_lst<-list()
  upper_lst<-list()
  for(i in seq_along(y_grp)){
    y<-y_grp[i]
    est_lst[[i]]<-unlist(plt_df[,paste0(y,".est")])
    lower_lst[[i]]<-unlist(plt_df[,paste0(y,".lower")])
    upper_lst[[i]]<-unlist(plt_df[,paste0(y,".upper")])
  }

  # plot forest
  plt_df %<>% rename(!!idx_display:=x_idx2)
  p <- forest(plt_df[,c(idx_display,tidy_col[-1])],
              est = est_lst,
              lower = lower_lst,
              upper = upper_lst,
              ci_column = plt_par$ci_column,
              ref_line = plt_par$ref_line,
              vert_line = plt_par$vert_line,
              arrow_lab = plt_par$arrow_lab,
              xlim = plt_par$xlim,
              ticks_at = plt_par$ticks_at,
              xlab = plt_par$xlab,
              nudge_y = plt_par$nudge_y,
              theme = tm)
  return(p)
}

# require(tidyverse, scales, patchwork, rlang)
plot_catBysites <- function(
    df,
    title = "Demographic distribution across sites",
    site_col = "site",
    cat_col  = "cat",
    cnt_col  = "cnt",
    site_order = NULL,
    cat_order = NULL,
    min_label_prop = 0.00,
    palette = NULL,
    show_count_axis = TRUE,
    label_offset_count_frac = 0.006,
    label_offset_prop = 0.006,
    count_label_size = 3,
    prop_label_size = 3,
    bar_height = 0.75
) {
  
  site_col <- sym(site_col)
  cat_col  <- sym(cat_col)
  cnt_col  <- sym(cnt_col)
  
  if (is.null(palette)) {
    palette <- c(
      "#E69F00", "#56B4E9", "#009E73", "#F0E442",
      "#0072B2", "#D55E00", "#CC79A7", "#999999"
    )
  }
  
  dat <- df %>%
    filter(
      !is.na(!!site_col),
      !is.na(!!cat_col),
      !is.na(!!cnt_col)
    ) %>%
    mutate(
      .site = as.character(!!site_col),
      .cat  = as.character(!!cat_col),
      .cnt  = as.numeric(!!cnt_col)
    ) %>%
    group_by(.site) %>%
    mutate(
      site_total = sum(.cnt, na.rm = TRUE),
      prop = .cnt / site_total
    ) %>%
    ungroup()
  
  if (is.null(site_order)) {
    site_order <- dat %>%
      distinct(.site, site_total) %>%
      arrange(site_total) %>%
      pull(.site)
  }
  
  if (is.null(cat_order)) {
    cat_order <- dat %>%
      distinct(.cat) %>%
      pull(.cat)
  }
  
  dat <- dat %>%
    mutate(
      .site = factor(.site, levels = site_order),
      .cat  = factor(.cat, levels = cat_order),
      .site_y = as.numeric(.site)
    ) %>%
    arrange(.site, .cat) %>%
    group_by(.site) %>%
    mutate(
      cnt_xmin = lag(cumsum(.cnt), default = 0),
      cnt_xmax = cumsum(.cnt),
      
      prop_xmin = lag(cumsum(prop), default = 0),
      prop_xmax = cumsum(prop)
    ) %>%
    ungroup()
  
  max_site_total <- max(dat$site_total, na.rm = TRUE)
  count_offset <- max_site_total * label_offset_count_frac
  
  dat <- dat %>%
    mutate(
      cnt_label_x = cnt_xmin + count_offset,
      prop_label_x = prop_xmin + label_offset_prop,
      
      cnt_label = if_else(
        prop >= min_label_prop,
        comma(.cnt),
        ""
      ),
      
      prop_label = if_else(
        prop >= min_label_prop,
        percent(prop, accuracy = 0.1),
        ""
      ),
      
      ymin = .site_y - bar_height / 2,
      ymax = .site_y + bar_height / 2
    )
  
  n_cat <- length(cat_order)
  
  p_cnt <- ggplot(dat) +
    geom_rect(
      aes(
        xmin = cnt_xmin,
        xmax = cnt_xmax,
        ymin = ymin,
        ymax = ymax,
        fill = .cat
      ),
      color = "white",
      linewidth = 0.25
    ) +
    geom_text(
      data = dat %>% filter(cnt_label != ""),
      aes(
        x = cnt_label_x,
        y = .site_y,
        label = cnt_label
      ),
      hjust = 0,
      size = count_label_size,
      color = "black"
    ) +
    scale_y_continuous(
      breaks = seq_along(site_order),
      labels = site_order,
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_x_continuous(
      labels = comma,
      limits = c(0, max_site_total * 1.08),
      expand = expansion(mult = c(0, 0.02))
    ) +
    scale_fill_manual(
      values = rep(palette, length.out = n_cat),
      drop = FALSE
    ) +
    labs(
      title = "Count scale",
      subtitle = "Bar length reflects site denominator size",
      x = "Count",
      y = NULL,
      fill = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.title = element_text(face = "bold")
    )
  
  p_prop <- ggplot(dat) +
    geom_rect(
      aes(
        xmin = prop_xmin,
        xmax = prop_xmax,
        ymin = ymin,
        ymax = ymax,
        fill = .cat
      ),
      color = "white",
      linewidth = 0.25
    ) +
    geom_text(
      data = dat %>% filter(prop_label != ""),
      aes(
        x = prop_label_x,
        y = .site_y,
        label = prop_label
      ),
      hjust = 0,
      size = prop_label_size,
      color = "black"
    ) +
    scale_y_continuous(
      breaks = seq_along(site_order),
      labels = site_order,
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_x_continuous(
      labels = percent_format(accuracy = 1),
      limits = c(0, 1.04),
      breaks = seq(0, 1, by = 0.25),
      expand = expansion(mult = c(0, 0.02))
    ) +
    scale_fill_manual(
      values = rep(palette, length.out = n_cat),
      drop = FALSE
    ) +
    labs(
      title = "Proportion scale",
      subtitle = "Bar length standardized to 100%",
      x = "Within-site proportion",
      y = NULL,
      fill = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(face = "bold")
    )
  
  if (!show_count_axis) {
    p_cnt <- p_cnt +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
  }
  
  p_cnt + p_prop +
    plot_layout(widths = c(1.25, 1), guides = "collect") +
    plot_annotation(title = title) &
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.title = element_text(face = "bold")
    )
}

# require(tidyverse, ggrepel, rlang, tidyr)
plot_numBysites <- function(
    df,
    selected_site,
    site_col = "site",
    rate_col = "incident_rate",
    event_col = "event",
    title = "Incidence rates by event",
    y_label = "Incidence rate",
    x_label = "Event",
    violin_alpha = 0.35,
    violin_width = 0.95,
    point_alpha = 0.65,
    point_size = 2.2,
    selected_point_size = 3.8,
    label_size = 3.4,
    jitter_width = 0.12,
    box_width = 0.35,
    show_box_outliers = FALSE,
    label_digits = 3,
    annotate_quartiles = TRUE,
    quartile_label_size = 2.8,
    selected_label_nudge_x = -0.28,
    quartile_label_nudge_x = 0.32,
    selected_label_fill = "white",
    quartile_label_fill = "gray97",
    text_size = 15
) {
  
  
  
  site_col  <- sym(site_col)
  rate_col  <- sym(rate_col)
  event_col <- sym(event_col)
  
  dat <- df %>%
    filter(
      !is.na(!!site_col),
      !is.na(!!rate_col),
      !is.na(!!event_col)
    ) %>%
    mutate(
      .site = as.character(!!site_col),
      .event = as.factor(!!event_col),
      .rate = as.numeric(!!rate_col),
      .is_selected = .site %in% selected_site
    )
  
  selected_dat <- dat %>%
    filter(.is_selected) %>%
    mutate(
      .label = paste0(
        .site,
        "\n",
        formatC(.rate, format = "f", digits = label_digits)
      )
    )
  
  quartile_dat <- dat %>%
    group_by(.event) %>%
    summarise(
      q1 = quantile(.rate, probs = 0.25, na.rm = TRUE),
      q3 = quantile(.rate, probs = 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = c(q1, q3),
      names_to = "quartile",
      values_to = ".rate"
    ) %>%
    mutate(
      quartile = toupper(quartile),
      .label = paste0(
        quartile,
        ": ",
        formatC(.rate, format = "f", digits = label_digits)
      )
    )
  
  p <- ggplot(dat, aes(x = .event, y = .rate)) +
    
    geom_violin(
      aes(fill = .event),
      width = violin_width,
      alpha = violin_alpha,
      trim = TRUE,
      color = "gray35",
      linewidth = 0.45,
      show.legend = FALSE
    ) +

    geom_boxplot(
      width = box_width,
      outlier.shape = ifelse(show_box_outliers, 19, NA),
      alpha = 0.85,
      color = "gray15",
      fill = "white",
      linewidth = 0.65
    ) +
    
    geom_jitter(
      aes(color = .event),
      width = jitter_width,
      height = 0,
      alpha = point_alpha,
      size = point_size,
      show.legend = FALSE
    )
  
  # Add Q1/Q3 first so selected-site labels can appear on top
  if (annotate_quartiles) {
    
    p <- p +
      geom_point(
        data = quartile_dat,
        aes(x = .event, y = .rate),
        inherit.aes = FALSE,
        shape = 23,
        fill = "gray90",
        color = "black",
        size = 2.6,
        stroke = 0.7
      ) +
      
      geom_label_repel(
        data = quartile_dat,
        aes(
          x = .event,
          y = .rate,
          label = .label
        ),
        inherit.aes = FALSE,
        size = quartile_label_size,
        fontface = "bold",
        label.size = 0.18,
        label.padding = unit(0.12, "lines"),
        label.r = unit(0.08, "lines"),
        fill = quartile_label_fill,
        color = "gray20",
        alpha = 0.88,
        nudge_x = quartile_label_nudge_x,
        direction = "y",
        box.padding = 0.35,
        point.padding = 0.20,
        min.segment.length = 0,
        segment.size = 0.22,
        segment.alpha = 0.55,
        max.overlaps = Inf
      )
  }
  
  # Add selected-site point and label last so they stay visually dominant
  p +
    geom_point(
      data = selected_dat,
      aes(x = .event, y = .rate),
      inherit.aes = FALSE,
      shape = 21,
      fill = "white",
      color = "black",
      stroke = 1.3,
      size = selected_point_size
    ) +
    
    geom_label_repel(
      data = selected_dat,
      aes(x = .event, y = .rate, label = .label),
      inherit.aes = FALSE,
      size = label_size,
      fontface = "bold",
      label.size = 0.35,
      label.padding = unit(0.20, "lines"),
      label.r = unit(0.12, "lines"),
      fill = selected_label_fill,
      color = "black",
      nudge_x = selected_label_nudge_x,
      direction = "y",
      box.padding = 0.65,
      point.padding = 0.35,
      min.segment.length = 0,
      segment.size = 0.35,
      segment.alpha = 0.85,
      max.overlaps = Inf
    ) +
    
    labs(
      title = title,
      x = x_label,
      y = y_label
    ) +
    
    coord_cartesian(ylim = c(0, NA)) +
    
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = text_size, face = "bold", angle = 35, hjust = 1),
      plot.title = element_text(face = "bold"),
      plot.margin = margin(10, 35, 10, 35)
    )
}

# require(survival, tidycmprsk, ggsurvfit, ggplot2, ggrepel, scales, broom, rlang, viridisLite)
plot_str_cif <- function(
    dat,
    time_days = "time_days",
    status = "status",
    strata = "strata",
    outcome = NULL,
    censor_level = NULL,
    status_levels = NULL,
    days_per_month = 365.25 / 12,
    x_break_months = 3,
    label_times_months = NULL,
    label_every_months = 12,
    conf.level = 0.95,
    ci_type = c("ribbon", "lines"),
    ci_alpha = 0.18,
    palette = NULL,
    label_ci = FALSE,
    label_accuracy = 0.1,
    axis_percent_accuracy = 1,
    linewidth = 0.9,
    title = NULL,
    subtitle = NULL,
    caption = NULL,
    
    # Highlight one observed stratum
    highlight_strata = NULL,
    highlight_linewidth = 1.8,
    highlight_linetype = "solid",
    highlight_alpha = 1,
    highlight_color = NULL,
    
    # Overlay pooled population-level CIF
    overlay_population = FALSE,
    population_label = "Overall",
    population_color = "black",
    population_linewidth = 2.4,
    population_linetype = "solid",
    population_alpha = 1,
    population_ci = FALSE,
    population_ci_alpha = 0.10,
    population_ci_fill = NULL,
    population_annual_labels = TRUE,
    population_label_prefix = TRUE,
    population_in_legend = TRUE,
    
    label_size = 3.2,
    label_fill = "white",
    legend_title = "Strata",
    xlab = "Follow-up time, months",
    ylab = "Cumulative incidence",
    x_limits = NULL,
    y_limits = c(0, NA),
    seed = 123,
    repel_max_overlaps = Inf,
    return_type = c("plot", "list")
) {
  ci_type <- match.arg(ci_type)
  return_type <- match.arg(return_type)
  
  required_pkgs <- c(
    "survival", "tidycmprsk", "ggsurvfit", "ggplot2",
    "ggrepel", "scales", "broom", "rlang", "viridisLite"
  )
  
  missing_pkgs <- required_pkgs[
    !vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
  ]
  
  if (length(missing_pkgs) > 0) {
    stop(
      "Please install missing package(s): ",
      paste(missing_pkgs, collapse = ", "),
      call. = FALSE
    )
  }
  
  check_logical_scalar <- function(x, nm) {
    if (!is.logical(x) || length(x) != 1 || is.na(x)) {
      stop("`", nm, "` must be TRUE or FALSE.", call. = FALSE)
    }
  }
  
  check_logical_scalar(label_ci, "label_ci")
  check_logical_scalar(overlay_population, "overlay_population")
  check_logical_scalar(population_ci, "population_ci")
  check_logical_scalar(population_annual_labels, "population_annual_labels")
  check_logical_scalar(population_label_prefix, "population_label_prefix")
  check_logical_scalar(population_in_legend, "population_in_legend")
  
  if (!is.data.frame(dat)) {
    dat <- as.data.frame(dat)
  }
  
  if (
    !is.numeric(days_per_month) ||
    length(days_per_month) != 1 ||
    !is.finite(days_per_month) ||
    days_per_month <= 0
  ) {
    stop("`days_per_month` must be a single positive numeric value.", call. = FALSE)
  }
  
  if (
    !is.numeric(x_break_months) ||
    length(x_break_months) != 1 ||
    !is.finite(x_break_months) ||
    x_break_months <= 0
  ) {
    stop("`x_break_months` must be a single positive numeric value.", call. = FALSE)
  }
  
  if (
    !is.numeric(label_every_months) ||
    length(label_every_months) != 1 ||
    !is.finite(label_every_months) ||
    label_every_months <= 0
  ) {
    stop("`label_every_months` must be a single positive numeric value.", call. = FALSE)
  }
  
  if (
    !is.numeric(conf.level) ||
    length(conf.level) != 1 ||
    !is.finite(conf.level) ||
    conf.level <= 0 ||
    conf.level >= 1
  ) {
    stop("`conf.level` must be a single number between 0 and 1.", call. = FALSE)
  }
  
  if (!is.null(x_limits)) {
    if (!is.numeric(x_limits) || length(x_limits) != 2) {
      stop("`x_limits` must be NULL or a numeric vector of length 2.", call. = FALSE)
    }
  }
  
  if (!is.null(y_limits)) {
    if (!is.numeric(y_limits) || length(y_limits) != 2) {
      stop("`y_limits` must be NULL or a numeric vector of length 2.", call. = FALSE)
    }
  }
  
  if (!is.null(highlight_color) && length(highlight_color) != 1) {
    stop("`highlight_color` must be NULL or a single color value.", call. = FALSE)
  }
  
  if (isTRUE(overlay_population)) {
    population_label <- as.character(population_label)
    
    if (
      length(population_label) != 1 ||
      is.na(population_label) ||
      !nzchar(population_label)
    ) {
      stop("`population_label` must be a single non-empty value.", call. = FALSE)
    }
    
    if (length(population_color) != 1) {
      stop("`population_color` must be a single color value.", call. = FALSE)
    }
    
    if (!is.null(population_ci_fill) && length(population_ci_fill) != 1) {
      stop("`population_ci_fill` must be NULL or a single color value.", call. = FALSE)
    }
  }
  
  time_nm <- rlang::as_name(rlang::ensym(time_days))
  status_nm <- rlang::as_name(rlang::ensym(status))
  strata_nm <- rlang::as_name(rlang::ensym(strata))
  
  missing_cols <- setdiff(c(time_nm, status_nm, strata_nm), names(dat))
  
  if (length(missing_cols) > 0) {
    stop(
      "dat is missing column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  
  time_raw <- dat[[time_nm]]
  status_raw <- dat[[status_nm]]
  strata_raw <- dat[[strata_nm]]
  
  time_num <- suppressWarnings(as.numeric(time_raw))
  
  if (any(is.na(time_num) & !is.na(time_raw))) {
    stop("`time_days` could not be converted cleanly to numeric.", call. = FALSE)
  }
  
  if (any(time_num < 0, na.rm = TRUE)) {
    stop("`time_days` contains negative follow-up times.", call. = FALSE)
  }
  
  status_chr <- as.character(status_raw)
  observed_status <- unique(status_chr[!is.na(status_chr)])
  
  # ---- Establish status levels with censoring first ----
  if (!is.null(status_levels)) {
    status_levels_chr <- as.character(status_levels)
    
    if (length(status_levels_chr) < 2) {
      stop(
        "`status_levels` must include at least one censoring level and one event level.",
        call. = FALSE
      )
    }
    
    if (anyDuplicated(status_levels_chr)) {
      stop("`status_levels` contains duplicated values.", call. = FALSE)
    }
  } else if (is.factor(status_raw)) {
    status_levels_chr <- levels(status_raw)
    
    if (length(status_levels_chr) < 2) {
      stop(
        "`status` must have at least one censoring level and one event level.",
        call. = FALSE
      )
    }
  } else {
    if (is.null(censor_level)) {
      stop(
        "For non-factor `status`, supply `censor_level`, e.g. ",
        "`censor_level = 0` or `censor_level = 'censor'`. ",
        "Alternatively, make `status` a factor with censoring as the first level.",
        call. = FALSE
      )
    }
    
    censor_chr_tmp <- as.character(censor_level)
    
    if (!censor_chr_tmp %in% observed_status) {
      stop(
        "`censor_level` was not found among observed status values. ",
        "If it is a valid but unobserved status level, supply `status_levels` explicitly.",
        call. = FALSE
      )
    }
    
    status_levels_chr <- c(
      censor_chr_tmp,
      setdiff(observed_status, censor_chr_tmp)
    )
  }
  
  if (is.null(censor_level)) {
    censor_chr <- status_levels_chr[1]
  } else {
    censor_chr <- as.character(censor_level)
    
    if (!censor_chr %in% status_levels_chr) {
      stop(
        "`censor_level` must be included in `status_levels` or in the factor levels of `status`.",
        call. = FALSE
      )
    }
    
    status_levels_chr <- c(
      censor_chr,
      setdiff(status_levels_chr, censor_chr)
    )
  }
  
  event_levels <- setdiff(status_levels_chr, censor_chr)
  
  if (length(event_levels) < 1) {
    stop("No event levels remain after defining the censoring level.", call. = FALSE)
  }
  
  if (is.null(outcome)) {
    outcome_chr <- event_levels[1]
    message("`outcome` not supplied; using first non-censoring level: ", outcome_chr)
  } else {
    outcome_chr <- as.character(outcome)
  }
  
  if (length(outcome_chr) != 1) {
    stop("This function expects a single `outcome`.", call. = FALSE)
  }
  
  if (!outcome_chr %in% event_levels) {
    stop(
      "`outcome` must be one of the non-censoring status levels: ",
      paste(event_levels, collapse = ", "),
      call. = FALSE
    )
  }
  
  status_fac <- factor(status_chr, levels = status_levels_chr)
  
  unknown_status <- unique(status_chr[is.na(status_fac) & !is.na(status_chr)])
  
  if (length(unknown_status) > 0) {
    stop(
      "These observed status values are not in `status_levels`: ",
      paste(unknown_status, collapse = ", "),
      call. = FALSE
    )
  }
  
  # ---- Strata ----
  strata_chr <- as.character(strata_raw)
  
  strata_levels <- if (is.factor(strata_raw)) {
    levels(strata_raw)
  } else {
    unique(strata_chr[!is.na(strata_chr)])
  }
  
  strata_fac <- factor(strata_chr, levels = strata_levels)
  
  base_df <- data.frame(
    .time_months = time_num / days_per_month,
    .status = status_fac,
    .strata = strata_fac
  )
  
  stratified_complete_idx <- stats::complete.cases(
    base_df[, c(".time_months", ".status", ".strata"), drop = FALSE]
  )
  
  population_complete_idx <- stats::complete.cases(
    base_df[, c(".time_months", ".status"), drop = FALSE]
  )
  
  if (!all(stratified_complete_idx)) {
    warning(
      "Dropping ",
      sum(!stratified_complete_idx),
      " row(s) from stratified curves because of missing time, status, or strata.",
      call. = FALSE
    )
  }
  
  plot_df <- base_df[stratified_complete_idx, , drop = FALSE]
  plot_df$.strata <- droplevels(plot_df$.strata)
  
  population_df <- base_df[
    population_complete_idx,
    c(".time_months", ".status"),
    drop = FALSE
  ]
  
  if (isTRUE(overlay_population)) {
    n_population_extra <- sum(population_complete_idx & !stratified_complete_idx)
    
    if (n_population_extra > 0) {
      message(
        "Population CIF includes ",
        n_population_extra,
        " row(s) with complete time/status but missing strata; ",
        "these row(s) are not included in the stratified curves."
      )
    }
  }
  
  strata_levels_obs <- levels(plot_df$.strata)
  
  if (nrow(plot_df) == 0) {
    stop("No complete rows remain for the stratified curves.", call. = FALSE)
  }
  
  if (isTRUE(overlay_population) && nrow(population_df) == 0) {
    stop("No complete rows remain for the pooled population curve.", call. = FALSE)
  }
  
  if (length(strata_levels_obs) < 1) {
    stop("No observed strata remain after preprocessing.", call. = FALSE)
  }
  
  if (!any(plot_df$.status == outcome_chr, na.rm = TRUE)) {
    stop(
      "No complete stratified-analysis rows have the requested outcome: ",
      outcome_chr,
      call. = FALSE
    )
  }
  
  if (
    isTRUE(overlay_population) &&
    !any(population_df$.status == outcome_chr, na.rm = TRUE)
  ) {
    stop(
      "No complete population-analysis rows have the requested outcome: ",
      outcome_chr,
      call. = FALSE
    )
  }
  
  max_month_stratified <- max(plot_df$.time_months, na.rm = TRUE)
  
  max_month_population <- if (isTRUE(overlay_population)) {
    max(population_df$.time_months, na.rm = TRUE)
  } else {
    NA_real_
  }
  
  max_month <- if (isTRUE(overlay_population)) {
    max(max_month_stratified, max_month_population, na.rm = TRUE)
  } else {
    max_month_stratified
  }
  
  # ---- Highlight validation ----
  highlight_chr <- NULL
  
  if (!is.null(highlight_strata)) {
    highlight_chr <- as.character(highlight_strata)
    
    if (length(highlight_chr) != 1) {
      stop("`highlight_strata` must be a single strata value.", call. = FALSE)
    }
    
    if (!highlight_chr %in% strata_levels_obs) {
      stop(
        "`highlight_strata` must be one of the observed strata values: ",
        paste(strata_levels_obs, collapse = ", "),
        call. = FALSE
      )
    }
  }
  
  # ---- Color-blind-friendly palette ----
  okabe_ito <- c(
    "#0072B2", # blue
    "#D55E00", # vermillion
    "#009E73", # bluish green
    "#CC79A7", # reddish purple
    "#E69F00", # orange
    "#56B4E9", # sky blue
    "#F0E442", # yellow
    "#000000"  # black
  )
  
  n_strata <- length(strata_levels_obs)
  
  if (is.null(palette)) {
    pal_vals <- if (n_strata <= length(okabe_ito)) {
      okabe_ito[seq_len(n_strata)]
    } else {
      viridisLite::viridis(n_strata, option = "D")
    }
    
    names(pal_vals) <- strata_levels_obs
  } else {
    pal_vals <- palette
    
    if (!is.null(names(pal_vals))) {
      missing_pal <- setdiff(strata_levels_obs, names(pal_vals))
      
      if (length(missing_pal) > 0) {
        stop(
          "Named `palette` is missing color(s) for strata: ",
          paste(missing_pal, collapse = ", "),
          call. = FALSE
        )
      }
      
      pal_vals <- pal_vals[strata_levels_obs]
    } else {
      if (length(pal_vals) < n_strata) {
        stop(
          "`palette` has fewer colors than the number of observed strata.",
          call. = FALSE
        )
      }
      
      pal_vals <- pal_vals[seq_len(n_strata)]
      names(pal_vals) <- strata_levels_obs
    }
  }
  
  population_key <- "..pooled_population_cif.."
  
  color_scale_vals <- pal_vals
  color_scale_breaks <- names(pal_vals)
  color_scale_labels <- names(pal_vals)
  
  legend_title_eff <- legend_title
  
  if (isTRUE(overlay_population)) {
    color_scale_vals <- c(
      color_scale_vals,
      stats::setNames(population_color, population_key)
    )
    
    if (isTRUE(population_in_legend)) {
      color_scale_breaks <- c(color_scale_breaks, population_key)
      color_scale_labels <- c(color_scale_labels, population_label)
      
      if (identical(legend_title, "Strata")) {
        legend_title_eff <- "Curve"
      }
    }
  }
  
  # ---- Small internal helpers ----
  make_time_label <- function(x) {
    is_exact_year <- abs(x / 12 - round(x / 12)) < 1e-8
    
    out <- paste0(
      scales::number(x, accuracy = 0.1),
      " mo"
    )
    
    out[is_exact_year] <- paste0(
      "Year ",
      scales::number(x[is_exact_year] / 12, accuracy = 1)
    )
    
    out
  }
  
  format_cif_labels <- function(df) {
    if (nrow(df) == 0) {
      return(df)
    }
    
    df$.time_label <- make_time_label(df$time)
    
    df$.estimate_label <- scales::percent(
      df$estimate,
      accuracy = label_accuracy
    )
    
    if (isTRUE(label_ci) && all(c("conf.low", "conf.high") %in% names(df))) {
      ci_level_label <- scales::percent(conf.level, accuracy = 1)
      
      ci_label <- ifelse(
        is.na(df$conf.low) | is.na(df$conf.high),
        "",
        paste0(
          "\n", ci_level_label, " CI ",
          scales::percent(df$conf.low, accuracy = label_accuracy),
          "–",
          scales::percent(df$conf.high, accuracy = label_accuracy)
        )
      )
      
      df$label <- paste0(
        df$.time_label,
        ": ",
        df$.estimate_label,
        ci_label
      )
    } else {
      df$label <- paste0(
        df$.time_label,
        ": ",
        df$.estimate_label
      )
    }
    
    df
  }
  
  add_step_endpoints <- function(
    curve,
    max_t,
    strata_label = NULL,
    curve_key = NULL,
    include_ci = FALSE
  ) {
    tol <- 1e-10
    
    if (!is.finite(max_t) || max_t < 0) {
      max_t <- 0
    }
    
    empty_curve <- function() {
      out <- data.frame(
        time = c(0, max_t),
        estimate = c(0, 0),
        outcome = rep(outcome_chr, 2)
      )
      
      if (!is.null(strata_label)) {
        out$strata <- factor(
          rep(as.character(strata_label), 2),
          levels = strata_levels_obs
        )
      }
      
      if (!is.null(curve_key)) {
        out$.curve_key <- rep(as.character(curve_key), 2)
      }
      
      if (isTRUE(include_ci)) {
        out$conf.low <- c(0, 0)
        out$conf.high <- c(0, 0)
      }
      
      out
    }
    
    if (
      nrow(curve) == 0 ||
      !all(c("time", "estimate") %in% names(curve))
    ) {
      return(empty_curve())
    }
    
    curve <- curve[
      !is.na(curve$time) & !is.na(curve$estimate),
      ,
      drop = FALSE
    ]
    
    if (nrow(curve) == 0) {
      return(empty_curve())
    }
    
    if (isTRUE(include_ci)) {
      if (!"conf.low" %in% names(curve)) {
        curve$conf.low <- NA_real_
      }
      
      if (!"conf.high" %in% names(curve)) {
        curve$conf.high <- NA_real_
      }
    }
    
    if (!"outcome" %in% names(curve)) {
      curve$outcome <- outcome_chr
    }
    
    if (!is.null(strata_label)) {
      curve$strata <- factor(
        as.character(strata_label),
        levels = strata_levels_obs
      )
    }
    
    if (!is.null(curve_key)) {
      curve$.curve_key <- as.character(curve_key)
    }
    
    curve <- curve[order(curve$time), , drop = FALSE]
    curve <- curve[!duplicated(curve$time, fromLast = TRUE), , drop = FALSE]
    
    if (!any(abs(curve$time - 0) < tol)) {
      start_row <- curve[1, , drop = FALSE]
      start_row$time <- 0
      start_row$estimate <- 0
      
      if (isTRUE(include_ci)) {
        start_row$conf.low <- 0
        start_row$conf.high <- 0
      }
      
      curve <- rbind(start_row, curve)
    } else if (isTRUE(include_ci)) {
      idx0 <- which(abs(curve$time - 0) < tol)
      
      bad_low <- idx0[is.na(curve$conf.low[idx0])]
      bad_high <- idx0[is.na(curve$conf.high[idx0])]
      
      if (length(bad_low) > 0) {
        curve$conf.low[bad_low] <- 0
      }
      
      if (length(bad_high) > 0) {
        curve$conf.high[bad_high] <- 0
      }
    }
    
    if (!any(abs(curve$time - max_t) < tol)) {
      end_row <- curve[nrow(curve), , drop = FALSE]
      end_row$time <- max_t
      curve <- rbind(curve, end_row)
    }
    
    curve <- curve[order(curve$time), , drop = FALSE]
    rownames(curve) <- NULL
    
    curve
  }
  
  make_step_ribbon <- function(curve) {
    if (!all(c("time", "conf.low", "conf.high") %in% names(curve))) {
      return(data.frame())
    }
    
    curve <- curve[
      is.finite(curve$time) &
        is.finite(curve$conf.low) &
        is.finite(curve$conf.high),
      ,
      drop = FALSE
    ]
    
    if (nrow(curve) < 2) {
      return(data.frame())
    }
    
    curve <- curve[order(curve$time), , drop = FALSE]
    curve <- curve[!duplicated(curve$time, fromLast = TRUE), , drop = FALSE]
    
    if (nrow(curve) < 2) {
      return(data.frame())
    }
    
    ribbon_list <- lapply(seq_len(nrow(curve) - 1), function(i) {
      data.frame(
        time = c(curve$time[i], curve$time[i + 1]),
        conf.low = rep(curve$conf.low[i], 2),
        conf.high = rep(curve$conf.high[i], 2)
      )
    })
    
    out <- do.call(rbind, ribbon_list)
    rownames(out) <- NULL
    
    out
  }
  
  # ---- Fit stratified cumulative incidence curves ----
  ci_fit <- tidycmprsk::cuminc(
    survival::Surv(.time_months, .status) ~ .strata,
    data = plot_df,
    conf.level = conf.level
  )
  
  # ---- Fit pooled population-level CIF, if requested ----
  population_fit <- NULL
  population_curve <- data.frame()
  population_ribbon <- data.frame()
  
  if (isTRUE(overlay_population)) {
    population_fit <- tidycmprsk::cuminc(
      survival::Surv(.time_months, .status) ~ 1,
      data = population_df,
      conf.level = conf.level
    )
    
    population_curve <- broom::tidy(
      population_fit,
      conf.int = isTRUE(population_ci),
      conf.level = conf.level
    )
    
    population_curve <- population_curve[
      as.character(population_curve$outcome) == outcome_chr &
        !is.na(population_curve$estimate),
      ,
      drop = FALSE
    ]
    
    population_curve <- add_step_endpoints(
      curve = population_curve,
      max_t = max_month_population,
      curve_key = population_key,
      include_ci = isTRUE(population_ci)
    )
    
    if (isTRUE(population_ci)) {
      population_ribbon <- make_step_ribbon(population_curve)
    }
  }
  
  # ---- Prepare highlighted stratum curve ----
  highlight_curve <- data.frame()
  
  if (!is.null(highlight_chr)) {
    curve_tidy <- broom::tidy(
      ci_fit,
      conf.int = FALSE
    )
    
    highlight_curve <- curve_tidy[
      as.character(curve_tidy$outcome) == outcome_chr &
        as.character(curve_tidy$strata) == highlight_chr,
      ,
      drop = FALSE
    ]
    
    max_highlight_month <- max(
      plot_df$.time_months[as.character(plot_df$.strata) == highlight_chr],
      na.rm = TRUE
    )
    
    if (!is.finite(max_highlight_month)) {
      max_highlight_month <- max_month_stratified
    }
    
    highlight_curve <- add_step_endpoints(
      curve = highlight_curve,
      max_t = max_highlight_month,
      strata_label = highlight_chr,
      include_ci = FALSE
    )
  }
  
  # ---- Time breaks and label times ----
  x_max <- ceiling(max_month / x_break_months) * x_break_months
  
  if (!is.finite(x_max) || x_max <= 0) {
    x_max <- x_break_months
  }
  
  if (is.null(x_limits)) {
    x_limits_eff <- c(0, x_max)
  } else {
    x_limits_eff <- x_limits
  }
  
  x_break_upper <- x_limits_eff[2]
  
  if (is.na(x_break_upper) || !is.finite(x_break_upper)) {
    x_break_upper <- x_max
  }
  
  x_breaks <- seq(
    from = 0,
    to = x_break_upper,
    by = x_break_months
  )
  
  if (length(x_breaks) == 0) {
    x_breaks <- 0
  }
  
  if (is.null(label_times_months)) {
    label_upper <- floor(max_month / label_every_months) * label_every_months
    
    label_times_months <- if (label_upper >= label_every_months) {
      seq(
        from = label_every_months,
        to = label_upper,
        by = label_every_months
      )
    } else {
      numeric(0)
    }
  }
  
  label_times_months <- sort(unique(as.numeric(label_times_months)))
  label_times_months <- label_times_months[is.finite(label_times_months)]
  label_times_months <- label_times_months[
    label_times_months >= 0 & label_times_months <= max_month
  ]
  
  # ---- Annual cumulative-incidence labels for stratified curves ----
  annual_lab <- data.frame()
  
  if (!is.null(highlight_strata)&length(label_times_months) > 0) {
    annual_lab <- broom::tidy(
      ci_fit,
      times = label_times_months,
      conf.int = TRUE,
      conf.level = conf.level
    )
    
    annual_lab <- annual_lab %>%
      filter(
        outcome == outcome_chr &
        !is.na(annual_lab$estimate) &
        strata == highlight_strata  
      )
    
    if (nrow(annual_lab) > 0) {
      annual_lab$strata <- factor(
        as.character(annual_lab$strata),
        levels = strata_levels_obs
      )
      
      annual_lab <- format_cif_labels(annual_lab)
    }
  }
  
  
  # ---- Annual cumulative-incidence labels for population curve ----
  population_lab <- data.frame()
  
  if (
    isTRUE(overlay_population) &&
    isTRUE(population_annual_labels) &&
    length(label_times_months) > 0
  ) {
    population_lab <- broom::tidy(
      population_fit,
      times = label_times_months,
      conf.int = TRUE,
      conf.level = conf.level
    )
    
    population_lab <- population_lab[
      as.character(population_lab$outcome) == outcome_chr &
        !is.na(population_lab$estimate),
      ,
      drop = FALSE
    ]
    
    if (nrow(population_lab) > 0) {
      population_lab <- format_cif_labels(population_lab)
      
      if (isTRUE(population_label_prefix)) {
        population_lab$label <- paste0(
          population_label,
          "\n",
          population_lab$label
        )
      }
      
      population_lab$.curve_key <- population_key
    }
  }
  
  # ---- Combine label data so ggrepel repels all labels together ----
  label_dfs <- list()
  
  if (nrow(annual_lab) > 0) {
    tmp <- annual_lab[, c("time", "estimate", "label"), drop = FALSE]
    tmp$.curve_key <- as.character(annual_lab$strata)
    label_dfs[[length(label_dfs) + 1]] <- tmp
  }
  
  if (nrow(population_lab) > 0) {
    tmp <- population_lab[, c("time", "estimate", "label", ".curve_key"), drop = FALSE]
    label_dfs[[length(label_dfs) + 1]] <- tmp
  }
  
  label_plot_df <- if (length(label_dfs) > 0) {
    do.call(rbind, label_dfs)
  } else {
    data.frame()
  }
  
  # ---- Plot ----
  p <- ggsurvfit::ggcuminc(
    ci_fit,
    outcome = outcome_chr,
    linewidth = linewidth
  ) +
    ggsurvfit::add_confidence_interval(
      type = ci_type,
      alpha = ci_alpha
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      caption = caption,
      x = xlab,
      y = ylab
    ) +
    ggplot2::scale_color_manual(
      values = color_scale_vals,
      breaks = color_scale_breaks,
      labels = color_scale_labels,
      name = legend_title_eff
    ) +
    ggplot2::scale_fill_manual(
      values = pal_vals,
      guide = "none"
    ) +
    ggplot2::scale_x_continuous(
      breaks = x_breaks,
      limits = x_limits_eff,
      name = xlab,
      guide = ggplot2::guide_axis(n.dodge = 2)
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(accuracy = axis_percent_accuracy),
      limits = y_limits,
      expand = ggplot2::expansion(mult = c(0, 0.12)),
      name = ylab
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      plot.margin = ggplot2::margin(5.5, 35, 5.5, 5.5)
    )
  
  # Optional pooled population CI ribbon.
  if (isTRUE(overlay_population) && isTRUE(population_ci) && nrow(population_ribbon) > 0) {
    population_ci_fill_eff <- if (is.null(population_ci_fill)) {
      population_color
    } else {
      population_ci_fill
    }
    
    p <- p +
      ggplot2::geom_ribbon(
        data = population_ribbon,
        mapping = ggplot2::aes(
          x = .data$time,
          ymin = .data$conf.low,
          ymax = .data$conf.high,
          group = 1
        ),
        inherit.aes = FALSE,
        fill = population_ci_fill_eff,
        alpha = population_ci_alpha,
        show.legend = FALSE
      )
  }
  
  # Draw the highlighted stratum curve on top of the base curves and CI.
  if (nrow(highlight_curve) > 0) {
    if (is.null(highlight_color)) {
      p <- p +
        ggplot2::geom_step(
          data = highlight_curve,
          mapping = ggplot2::aes(
            x = .data$time,
            y = .data$estimate,
            color = .data$strata
          ),
          inherit.aes = FALSE,
          linewidth = highlight_linewidth,
          linetype = highlight_linetype,
          alpha = highlight_alpha,
          show.legend = FALSE
        )
    } else {
      p <- p +
        ggplot2::geom_step(
          data = highlight_curve,
          mapping = ggplot2::aes(
            x = .data$time,
            y = .data$estimate
          ),
          inherit.aes = FALSE,
          color = highlight_color,
          linewidth = highlight_linewidth,
          linetype = highlight_linetype,
          alpha = highlight_alpha,
          show.legend = FALSE
        )
    }
  }
  
  # Draw the pooled population curve last so it is visually emphasized.
  if (isTRUE(overlay_population) && nrow(population_curve) > 0) {
    p <- p +
      ggplot2::geom_step(
        data = population_curve,
        mapping = ggplot2::aes(
          x = .data$time,
          y = .data$estimate,
          color = .data$.curve_key
        ),
        inherit.aes = FALSE,
        linewidth = population_linewidth,
        linetype = population_linetype,
        alpha = population_alpha,
        show.legend = isTRUE(population_in_legend)
      )
  }
  
  # Add repelled annual labels after all curve layers.
  if (nrow(label_plot_df) > 0) {
    p <- p +
      ggrepel::geom_label_repel(
        data = label_plot_df,
        mapping = ggplot2::aes(
          x = .data$time,
          y = .data$estimate,
          label = .data$label,
          color = .data$.curve_key
        ),
        inherit.aes = FALSE,
        fill = label_fill,
        label.size = 0.2,
        size = label_size,
        min.segment.length = 0,
        seed = seed,
        max.overlaps = repel_max_overlaps,
        show.legend = FALSE
      )
  }
  
  if (return_type == "list") {
    return(
      list(
        plot = p,
        fit = ci_fit,
        population_fit = population_fit,
        annual_labels = annual_lab,
        population_labels = population_lab,
        label_plot_data = label_plot_df,
        highlight_curve = highlight_curve,
        population_curve = population_curve,
        population_ribbon = population_ribbon,
        analysis_data = plot_df,
        population_analysis_data = population_df,
        outcome = outcome_chr,
        censor_level = censor_chr,
        status_levels = status_levels_chr,
        palette = pal_vals,
        color_scale_values = color_scale_vals,
        highlight_strata = highlight_chr,
        overlay_population = overlay_population,
        population_label = if (isTRUE(overlay_population)) population_label else NULL
      )
    )
  } else {
    return(p)
  }
}




