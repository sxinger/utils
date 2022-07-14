
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
                      cdm_tbl = "DIAGNOSIS",
                      vs_url = "",
                      vs_type = "curated-dx",
                      vs_file_type = "json",
                      vs_name = "",
                      key_col_schema=c("PATID"),
                      col_out = NULL,
                      dry_run = TRUE){
  # get column name of the target cdm table
  cdm_tbl_nm<-paste0("identifier($",cdm_tbl,")")
  if(is.null(col_out)){
    col_out<-colnames(DBI::dbGetQuery(
      conn,paste0("select * from ",cdm_schema,".",cdm_tbl_nm," where 1=0")
    ))
  }
  col_out<-col_out[!col_out %in% key_col_schema]
  # contruct matching phrase
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
                            paste0(" from ", patset, " p"),
                            paste0(" join ", cdm_schema, ".",cdm_tbl_nm," f"),
                            " on ",match_key)))
  # specify where clause based on valuset type
  if (vs_file_type == "json"){
    vs<-jsonlite::fromJSON(file=vs_url)[[vs_name]]
    # curated valueset of icd diagnosis codes
    if (vs_type == "curated-dx"){
      sql_where<-"("      
      for (icd_vrsn in c("icd9","icd10")){
        for(lev_i in 0:2){
          icd_lev<-paste0("lev",lev_i)
          if (lev_i==0){
            sql_where<-paste0(sql_where,"split_part(f.DX,'.',1) in ", 
                         paste0("('",paste(vs[[icd_vrsn]][[icd_lev]],collapse="','"),"')")," or ")
          }else{
            sql_where<-paste0(sql_where,"substr(f.DX,1,",4+lev_i,") in ", 
                         paste0("('",paste(vs[[icd_vrsn]][[icd_lev]],collapse="','"),"')")," or ")
          }
        }
      }
      sql<-paste0(sql," and ",sql_where,")")
    # curated valueset of cpt or icd procedure codes
    } else if (vs_type == "curated-px"){
      # icd pcs codes
      sql_where<-"(" 
      for(icd_vrsn in c("icd9","icd10")){
        for(lev_i in 0:2){
          icd_lev<-paste0("lev",lev_i)
          if(lev_i==0&&length(vs[[icd_vrsn]][[icd_lev]])>0){
            sql_where<-paste0(sql_where,"split_part(f.DX,'.',1) in ", 
                         paste0("('",paste(vs[[icd_vrsn]][[icd_lev]],collapse="','"),"')")," or ")
          }else if(length(vs[[icd_vrsn]][[icd_lev]])>0){
            sql_where<-paste0(sql_where,"substr(f.DX,1,",4+lev_i,") in ", 
                         paste0("('",paste(vs[[icd_vrsn]][[icd_lev]],collapse="','"),"')")," or ")
          }else{
            stop("valueset is empty!")
          }
        }
      # cpt codes
      for(type in c("range","list")){
          if(type=="range"&&length(vs[[icd_vrsn]][[icd_lev]])>0){
            rg2vec<-lapply(vs[["hcpcs"]][["range"]],
                           function(x)seq(as.numeric(strsplit(x,"-")[[1]][1]),as.numeric(strsplit(x,"-")[[1]][2])))
            sql_where<-paste0(sql_where,"f.PX in ", 
                         paste0("('",paste(unlist(rg2vec),collapse="','"),"')")," or ")
          }else if(length(vs[[icd_vrsn]][[icd_lev]])>0){
            sql_where<-paste0(sql_where,"f.PX in ", 
                         paste0("('",paste(vs[[icd_vrsn]][[icd_lev]],collapse="','"),"')")," or ")
          }else{
            stop("valueset is empty!")
          }
        }
      }
      sql<-paste0(sql," and ",sql_where,")")
    } else if (vs_type == "ecqm") {
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
        }else{
          stop("ontology not used in CDM!")
        }
      }
    } else if (vs_type == "ccda"){
      # TODO
    } else if (vs_type == "ncbo"){
      # TODO
    } else if (vs_type == "umls"){
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
    DBI::dbGetQuery(conn,sql)
  }
}












