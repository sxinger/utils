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

