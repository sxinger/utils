
#######################################################
## project-specific functions for HERON data request ##
#######################################################
#Note: all following functions are specific to projects with dependencies

#purpose: collect data from target CDM table based on provided valueset
#dependency: 1) patient set in database; 2) valueset specified in json file;
#            3) packages: DBI,jsonlite

collect_cdm<-function(conn,
                      patset = "PAT_ELIG",
                      key_col_pat = c("PATID"),
                      cdm_schema = "CMS_PCORNET_CDM",
                      key_col_schema=c("PATID"),
                      cdm_tbl = "DIAGNOSIS",
                      idx_dt_col = "INDEX_DATE",
                      ep_dt_col = "EVENT_DATE",
                      vs_url = "",
                      vs_template = "curated-dx",
                      vs_name = "",
                      col_out = NULL,
                      write_back = FALSE,
                      home_schema = "PUBLIC",
                      write_tbl = "TEMP",
                      dry_run = TRUE){
  
  # get column name of the target cdm table
  if(is.null(col_out)){
    col_out<-colnames(DBI::dbGetQuery(
      conn,paste0("select * from ",cdm_schema,".",cdm_tbl," where 1=0")
    ))
  }
  col_out<-col_out[!col_out %in% key_col_schema]

  # construct matching phrase
  match_key<-c()
  for(i in seq_along(key_col_pat)){
    match_key<-c(match_key,paste(c(paste0("p.",key_col_pat[i]),
                                   paste0("f.",key_col_schema[i])),
                                 collapse = "="))
  }
  if(length(key_col_pat)>1){
    match_key<-paste(match_key,collapse = " and ")
  }
  
  # construct sql baseline 
  sql<-paste0("select distinct ",
              paste0(paste(paste0("p.",key_col_pat,collapse = ","),",",
                           paste(paste0("f.",col_out,collapse = ",")),
                           paste0(",datediff('day',p.",idx_dt_col,",f.",ep_dt_col,") as DAYS_SINCE_INDEX"),
                           paste0(" from ", home_schema, ".", patset, " p"),
                           paste0(" join ", cdm_schema, ".", cdm_tbl," f"),
                           " on ",match_key)))
  
  vs_file_type<-gsub(".*\\.","",vs_url)
  # specify "where" clause based on valuset type
  if (vs_file_type == "json"){
    vs<-jsonlite::fromJSON(vs_url)[[vs_name]]
    # curated valueset of icd diagnosis codes
    if (vs_template == "curated-dx"){
      sql_where<-""      
      for (icd_vrsn in c("icd9","icd10")){
        for(lev_i in 0:2){
          icd_lev<-paste0("lev",lev_i)
          if (lev_i==0&&length(vs[[icd_vrsn]][[icd_lev]])>0){
            sql_where<-paste0(sql_where," or split_part(f.DX,'.',1) in ", 
                         paste0("('",paste(vs[[icd_vrsn]][[icd_lev]],collapse="','"),"')"))
          }else if(length(vs[[icd_vrsn]][[icd_lev]])>0){
            sql_where<-paste0(sql_where," or substr(f.DX,1,",4+lev_i,") in ", 
                         paste0("('",paste(vs[[icd_vrsn]][[icd_lev]],collapse="','"),"')"))
          }
        }
      }
      sql<-paste0(sql," and (",gsub("^( or)+","",sql_where),")")
      
    # curated valueset of cpt or icd procedure codes
    } else if (vs_template == "curated-px"){
      # icd pcs codes
      sql_where<-"" 
      for(icd_vrsn in c("icd9-pcs","icd10-pcs")){
        for(lev_i in 0:2){
          icd_lev<-paste0("lev",lev_i)
          if(lev_i==0&&length(vs[[icd_vrsn]][[icd_lev]])>0){
            sql_where<-paste0(sql_where," or split_part(f.PX,'.',1) in ", 
                         paste0("('",paste(vs[[icd_vrsn]][[icd_lev]],collapse="','"),"')"))
          }else if(length(vs[[icd_vrsn]][[icd_lev]])>0){
            sql_where<-paste0(sql_where," or substr(f.PX,1,",3+lev_i,") in ", 
                         paste0("('",paste(vs[[icd_vrsn]][[icd_lev]],collapse="','"),"')"))
          }
        }
      # cpt codes
      for(type in c("range","list")){
          if(type=="range"&&length(vs[[icd_vrsn]][[icd_lev]])>0){
            rg2vec<-lapply(vs[["hcpcs"]][["range"]],
                           function(x)seq(as.numeric(strsplit(x,"-")[[1]][1]),as.numeric(strsplit(x,"-")[[1]][2])))
            sql_where<-paste0(sql_where," or f.PX in ", 
                         paste0("('",paste(unlist(rg2vec),collapse="','"),"')"))
          }else if(length(vs[[icd_vrsn]][[icd_lev]])>0){
            sql_where<-paste0(sql_where," or f.PX in ", 
                         paste0("('",paste(vs[[icd_vrsn]][[icd_lev]],collapse="','"),"')"))
          }
        }
      }
      sql<-paste0(sql," and (",gsub("^( or)+","",sql_where),")")
    
    # eCQM published valuesets
    } else if (vs_template == "ecqm") {
      for(cdtype in names(vs[["codelist"]])){
        codelist<-paste0("('",paste(vs[["codelist"]][[cdtype]][,1],collapse="','"),"')")
        # diagnosis
        if(grepl("^(ICD)+",cdtype)){
          if(cdm_table=="DIAGNOSIS"){
            sql<-paste0(sql," and f.DX in ",codelist)
          }else if(cdm_table=="CONDITION"){
            sql<-paste0(sql," and f.CONDITION in ",codelist)
          }
        # procedures
        }else if((grepl("^((CPT)|(HCPCS))+",cdtype))){
          if(cdm_table=="PROCEDURES"){
            sql<-paste0(sql," and f.PX in ",codelist)
          }
        # lab
        }else if(grepl("^(LOINC)+",cdtype)){
          if(cdm_table=="LAB_RESULT_CM"){
            sql<-paste0(sql," and f.LAB_LOINC in ",codelist)
          }else if(cdm_table=="OBS_CLIN"){
            sql<-paste0(sql," and f.OBSLINC_CODE in ",codelist,
                        "where OBSCLIN_TYPE ='LC'")
          }
        # meds
        }else if(grepl("^(RXNORM)+",cdtype)){
          if(cdm_table=="PRESCRIBING"){
            ssql<-paste0(sql," and f.RXNORM_CUI in ",codelist)
          }else if(cdm_table=="MED_ADMIN"){
            sql<-paste0(sql," and f.MEDADMIN_CODE in ",codelist,
                      "where MEDADMIN_TYPE = 'RX'")
          }          
        # obs_clin
        }else if(grepl("^(SNOMED)+",cdtype)){
          if(cdm_table=="DIAGNOSIS"){
            sql<-paste0(sql," and f.DX in ",codelist,
                      "where DX_TYPE = 'SM'")
          }else if(cdm_table=="CONDITION"){
            sql<-paste0(sql," and f.CONDITION in ",codelist,
                      "where CONDITION_TYPE = 'SM'")
          }else if(cdm_table=="OBS_CLIN"){
            sql<-paste0(sql," and f.OBSLINC_CODE in ",codelist,
                        "where OBSCLIN_TYPE ='SM'")
          }  
        }
      }
      
    # CCDA published valuesets
    } else if (vs_template == "ccda"){
      # TODO
    
    # bioportal search results
    } else if (vs_template == "ncbo"){
      # TODO
      
    # UMLS API search results
    } else if (vs_template == "umls"){
      # TODO
    } else {
      stop("valueset file type not supported!")
    }
  } else if (vs_file_type == "csv"){
    #TODO
  }
  
  # run query
  if(dry_run==TRUE){
    return(sql)
  }else{
    if(write_back==TRUE){
      DBI::dbSendQuery(conn,
                       paste0("create or replace table ",
                              home_schema,".",write_tbl," as ",sql))
    }else{
      DBI::dbGetQuery(conn,sql)
    }
  }
}
