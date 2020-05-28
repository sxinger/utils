###############################
## general utility functions ##
###############################
##some initial setup
startup_config<-function(){
  #create .Renviron file 
  
}


##install and load multiple libraries
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

unzip_redcap<-function(req_nm){
  #unzip and save
  zipF<-paste0("./",req_nm,"-raw.zip")
  unzip(zipF)
  
  #load date
  out<-list(pat_tbl=read.csv(paste0("./",req_nm,"/",req_nm,"-patient.csv"),stringsAsFactors = F,na.strings = c(""," ")),
            data_tbl=read.csv(paste0("./",req_nm,"/",req_nm,"-data.csv"),stringsAsFactors = F,na.strings = c(""," ")),
            code_info=read.csv(paste0("./",req_nm,"/",req_nm,"-code.info.csv"),stringsAsFactors = F,na.strings = c(""," ")))
  
  return(out)  
}

