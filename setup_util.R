###############################
## general utility functions ##
###############################
##some initial setup
# startup_config<-function(){
#   #create .Renviron file 
# }

## install and load multiple libraries (equivalent to pacman::p_load)
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

## remove all 
remove_libraries<-function(pkg_lst=c()){
  # create a list of all installed packages
  ip <- as.data.frame(installed.packages()) %>%
    # if you use MRO, make sure that no packages in this library will be removed
    filter(!grepl("MRO", LibPath)) %>%
    # we don't want to remove base or recommended packages
    filter(Priority %in% c("base", "recommended"))
  if (length(pkg_lst) > 0 ){
    ip %<>% filter(Package %in% pkg_lst)
  }
  # remove the packages
  sapply(ip$Package, remove.packages, lib = ip$LibPath)
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

