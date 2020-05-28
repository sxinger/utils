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

##################################################
## utility functions for database interrogation ##
##################################################
#----collect facts from i2b2 observation_fact table----
collect_i2b2_obs<-function(conn,
                           code_vec=c(),
                           regexp_str="",
                           col_out=c("patient_num",
                                     "encounter_num",
                                     "concept_cd",
                                     "units_cd",
                                     "nval_num",
                                     "tval_char",
                                     "modifier_cd",
                                     "start_date",
                                     "end_date"),
                           key_col=c("patient_num"),
                           schema=c("blueherondata"),
                           pat_num){
  
  col_out<-col_out[!col_out %in% key_col]
  
  match_key<-c()
  for(i in seq_along(key_col)){
    match_key<-c(match_key,paste(paste0(c("p.","f."),key_col[i]),collapse = "="))
  }
  if(length(key_col)>1){
    match_key<-paste(match_key,collapse = " and ")
  }
  
  i2b2_obs_lst<-list()
  
  for(i in seq_along(schema)){
    sql<-paste0("select distinct ",
                paste0(paste(paste0("p.",key_col,collapse = ","),",",
                             paste(paste0("f.",col_out,collapse = ",")),
                             paste0(" from ", pat_num, " p"),
                             paste0(" join ", schema, ".observation_fact f"),
                             " on ",match_key)))
    
    if(length(code_vec)>0&nchar(regexp_str)>0){
      sql<-paste0(sql," and ",
                  "(f.concept_cd in ",paste0("('",paste(code_vec,collapse="','"),"')")," or ",
                  "regexp_like(f.concept_cd,'",regexp_str,"','i'))")
    }else if(length(code_vec)==0&nchar(regexp_str)>0){
      sql<-paste0(sql," and ",
                  "regexp_like(f.concept_cd,'",regexp_str,"','i')")
    }else if(length(code_vec)>0&nchar(regexp_str)==0){
      sql<-paste0(sql," and ",
                  "f.concept_cd in ",paste0("('",paste(code_vec,collapse="','"),"')"))
    }else{
      stop("either code_vec or regexp_str should be specified for filtering concept_cd!")
    }
    
    i2b2_obs_lst[[schema[i]]]<-DBI::dbGetQuery(conn,sql)
  }
  
  return(i2b2_obs_lst)
}

#----collect concepts from i2b2 concept_dimension table----
collect_i2b2_cd<-function(conn,
                          exact_match=T,
                          cd_prefix=NULL,
                          col_out=c("concept_cd",
                                    "name_char",
                                    "concept_path"),
                          key_col=c("CONCEPT_CD"),
                          schema=c("blueherondata"),
                          concept){
  
  col_out<-col_out[!col_out %in% key_col]
  
  match_key<-c()
  for(i in seq_along(key_col)){
    if(exact_match){
      match_key<-c(match_key,paste(paste0(c("cd.","f."),key_col[i]),collapse = "="))
    }else{
      match_key<-c(match_key,paste0("regexp_like(f.",key_col[i],",('(' || cd.",key_col[i]," || ')+'),'i')"))
    }
    
    if(!is.null(cd_prefix)){
      match_key<-paste0(match_key," and regexp_like(f.",key_col[i],",'^(",cd_prefix,")+')")
    }else{
      
    }
  }
  
  if(length(key_col)>1){
    match_key<-paste(match_key,collapse = " and ")
  }
  
  i2b2_obs_lst<-list()
  
  for(i in seq_along(schema)){
    sql<-paste0("select distinct ",
                paste0(paste(paste0("cd.",key_col,collapse = ",")," ICD_FUZZY,",
                             paste(paste0("f.",col_out,collapse = ",")),
                             paste0(" from ", concept, " cd"),
                             paste0(" join ", schema, ".concept_dimension f"),
                             " on ",match_key)))
    
    i2b2_obs_lst[[schema[i]]]<-DBI::dbGetQuery(conn,sql)
  }
  
  return(i2b2_obs_lst)
}


#----collect data from one of the CDM tables----
collect_cdm<-function(conn,
                      code_vec=c(),
                      str_vec=c(),
                      col_out=NULL,
                      key_col_schema=c("PATID"),
                      key_col_pat=c("PATIENT_NUM"),
                      schema=c("PCORNET_CDM_C7R2"),
                      tbl="DEMOGRAPHIC",
                      pat_num){
  if(is.null(col_out)){
    col_out<-colnames(DBI::dbGetQuery(conn,
                                      paste0("select * from ",schema[1],".",tbl," where 1=0")))
  }
  col_out<-col_out[!col_out %in% key_col_schema]
  
  
  match_key<-c()
  for(i in seq_along(key_col_pat)){
    match_key<-c(match_key,paste(c(paste0("p.",key_col_pat[i]),
                                   paste0("f.",key_col_schema[i])),
                                 collapse = "="))
  }
  if(length(key_col_pat)>1){
    match_key<-paste(match_key,collapse = " and ")
  }
  
  cdm_obs_lst<-list()
  
  for(i in seq_along(schema)){
    sql<-paste0("select distinct ",
                paste0(paste(paste0("p.",key_col_pat,collapse = ","),",",
                             paste(paste0("f.",col_out,collapse = ",")),
                             paste0(" from ", pat_num, " p"),
                             paste0(" join ", schema[i], ".",tbl," f"),
                             " on ",match_key)))
    if(tbl=="PROCEDURES"){
      sql<-paste0(sql," and",
                  " f.PX in ",paste0("('",paste(code_vec,collapse="','"),"')"))
    }else if(tbl=="PRESCRIBING"){
      if(length(str_vec)>0){
        sql<-paste0(sql," and ",
                    "(f.RXNORM_CUI in ",paste0("('",paste(code_vec,collapse="','"),"')")," or ",
                    " regexp_like(f.RAW_RX_MED_NAME,",paste0("'(",paste(str_vec,collapse = ")|("),")+'"),",'i'))")
      }else{
        sql<-paste0(sql," and ",
                    "f.RXNORM_CUI in ",paste0("('",paste(code_vec,collapse="','"),"')"))
      }
    }else if(tbl=="DISPENSING"){
      sql<-paste0(sql," and ",
                  "(f.NDC in ",paste0("('",paste(code_vec,collapse="','"),"')"))
    }else{
      sql<-sql
    }
    
    cdm_obs_lst[[schema[i]]]<-DBI::dbGetQuery(conn,sql)
  }
  
  return(cdm_obs_lst)
  
}

