###############################
## general utility functions ##
###############################
##some initial setup
# startup_config<-function(){
#   #create .Renviron file 
# }

##install and load multiple libraries (equivalent to pacman::p_load)
require_libraries<-function(package_list,verb=T){
  for (lib in package_list) {
    chk_install<-!(lib %in% installed.packages()[,"Package"])
    if(chk_install){
      install.packages(lib)
    }
    library(lib, character.only=TRUE,lib.loc=.libPaths())
    if(verb){
      cat("\n", lib, " loaded.", sep="") 
    }
  }
}

## dual-axis utility function
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

