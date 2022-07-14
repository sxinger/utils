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
## example
# ggplot(df %<>%
#          mutate(new_col = convert_scale(enc_pat_ratio,enc)[["val"]],
#                 axis_formula = convert_scale(enc_pat_ratio,enc)[["formula"]]),
#        aes(x=site)) +
#   geom_col(aes(y=enc), size = 1, color = "darkblue", fill = "white")+
#   geom_line(aes(y=new_col), size = 1.5, color="red", group = 1)+
#   scale_y_continuous(sec.axis = sec_axis(as.formula(df$axis_formula[1]), name = "sec_axis"))+
#   theme(axis.text.x = element_text(angle = 45),text=element_text(face="bold")) +
#   facet_wrap(~ enc_type,ncol=2,scales = "free")
