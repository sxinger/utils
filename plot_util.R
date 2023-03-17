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
  y,    # primary y-axis
  y_sec # secondary y-axis
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
  tm = forest_theme(
    arrow_type = "closed",
    arrow_label_just = "end"
  ) # theme parameters for forest function
){
  # require(tidyverse,grid,forestploter)
  # https://cran.r-project.org/web/packages/forestploter/vignettes/forestploter-intro.html

  plt_par$ci_column=Filter(Negate(is.null),list(plt_par$ci_column,2*seq_len(ny)))[[1]]
  plt_par$ref_line = Filter(Negate(is.null),list(plt_par$ref_line,rep(1, ny)))[[1]]
  plt_par$vert_line = Filter(Negate(is.null),list(plt_par$vert_line,rep(list(0.3, 1.4),ny)))[[1]]
  plt_par$arrow_lab = Filter(Negate(is.null),list(plt_par$arrow_lab,rep(list("Low", "High"),ny)))[[1]]
  plt_par$xlim = Filter(Negate(is.null),list(plt_par$xlim,rep(list(0, 3),ny)))[[1]]
  plt_par$x_trans = Filter(Negate(is.null),list(plt_par$x_trans,rep("none",ny)))[[1]]
  plt_par$ticks_at = Filter(Negate(is.null),list(plt_par$ticks_at,rep(list(0.1, 0.5, 1, 2.5),ny)))[[1]]
  plt_par$xlab = Filter(Negate(is.null),list(plt_par$xlab,rep("HR", ny)))[[1]]
  plt_par$nudge_y = c(plt_par$nudge_y,0.2)[1]

  # change to internal names for easy reference
  nm_map<-data.frame(
    ext_nm=c(x_idx1,x_idx2,y_idx,est,lower,upper,pval)
  ) %>%
    mutate(int_nm=deparse(substitute(ext_nm)))
  plt_df<-df %>%
    select(all_of(nm_map$ext_nm)) %>%
    rename_at(vars(nm_map$ext_nm), ~ nm_map$int_nm)
  
  # add an empty column for HR plots and a label column
  plt_df %<>%
    mutate(
      pvalstar = case_when(pval > 0.1 ~ "",
                           pval <= 0.1 & pval > 0.05 ~ "*",
                           pval <= 0.05 & pval > 0.01 ~ "**",
                           pval <= 0.01 & pval > 0.001 ~ "***",
                           TRUE ~ "****"),
      Outcome = paste(rep(" ", 20), collapse = " "),
      `HR (95% CI)` = sprintf("%.2f (%.2f to %.2f) %s",est, lower, upper, pvalstar)
    )
  
  # pivot wide
  y_grp<-unique(plt_df$y_idx)
  n_grp<-length(y_grp)
  plt_df %<>%
    pivot_wider(
      id_cols = c("x_idx1","x_idx2"),
      names_from = "y_idx",
      names_sep = ".",
      values_from = c(Outcome,est,lower,upper,`HR (95% CI)`),
    ) %>%
    group_by(x_idx1,x_idx2) %>%
    mutate(idx=order(x_idx2)) %>%
    ungroup
  
  # add header rows
  plt_df %<>%
    bind_rows(data.frame(idx=0,x_idx1 = unique(plt_df$x_idx1))) %>%
    mutate(x_idx2 = paste0("  ",x_idx1)) %>%
    arrange(x_idx1,idx) 

  # collect list for HR segment plots
  est_lst<-list()
  lower_lst<-list()
  upper_lst<-list()
  for(y in seq_along(y_grp)){
    est_lst<-list(est_lst,plt_df[,paste0(y,".est")])
    lower_lst<-list(lower_lst,plt_df[,paste0(y,".lower")])
    upper_lst<-list(upper_lst,plt_df[,paste0(y,".upper")])
  }

  # plot forest
  tidy_col<-c("x_idx2",list(a=c("Outcome","HR (95% CI)"),b=y_grp) |> cross() |> map_chr(\(...) paste0(..., collapse = ".")))
  p <- forest(plt_df[,tidy_col],
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