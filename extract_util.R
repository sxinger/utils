#######################################################
# search CDM table with imported valuesets
#######################################################

cdm_code_type_map<-function(cdtype){
  if(grepl("icd9",tolower(cdtype))){
    return("09")
  }else if(grepl("icd10",tolower(cdtype))){
    return("10")
  }else if(grepl("loinc",tolower(cdtype))){
    return("LC")
  }else if(grepl("rxn",tolower(cdtype))){
    return("RX")
  }else if(grepl("ndc",tolower(cdtype))){
    return("ND")
  }else if(grepl("snom",tolower(cdtype))){
    return("SM")
  }else if(grepl("hcpcs",tolower(cdtype))){
    return("CH")
  }else if(grepl("cpt",tolower(cdtype))){
    return("CH")
  }else{
    stop("Standardized value doesn't exist for code type:",cdtype)
  }
}

load_valueset.ncbo<-function(vs_url = "",vs_name_str = ""){
  # load valueset in json
  vs_file_type<-gsub(".*\\.","=",vs_url)
  vs_file<-jsonlite::fromJSON(vs_url)
  
  # initialize lookup table
  lookup_tbl<-data.frame(CODE_TYPE=as.character(),
                         CODE_TYPE_CDM=as.character(),
                         CODE_SUBTYPE=as.character(),
                         CODE=as.character(),
                         CODE_GRP=as.character(),
                         stringsAsFactors=F)
  
  # main code body for parsing json file
  vs_name_list<-names(vs_file)
  vs_name_dist<-stringdist(tolower(vs_name_str),vs_name_list, method="jw")
  vs_name_match<-vs_name_list[which.min(vs_name_dist)]
  vs<-vs_file[[vs_name_match]]
  for(cd_type_idx in seq_along(vs[["code_type"]])){
    lookup_tbl %<>%
      bind_rows(data.frame(CODE_TYPE=vs[["code_type"]][[cd_type_idx]],
                           CODE_TYPE_CDM=cdm_code_type_map(vs[["code_type"]][cd_type_idx]),
                           CODE_SUBTYPE="exact",
                           CODE=as.character(vs[["code"]][[cd_type_idx]]),
                           CODE_GRP=vs_name_match,
                           stringsAsFactors = F))
  }
  # return data.frame
  return(lookup_tbl)
}

load_valueset.curated<-function(vs_url = "",vs_name_str = ""){
  # load valueset in json
  vs_file_type<-gsub(".*\\.","=",vs_url)
  vs_file<-jsonlite::fromJSON(vs_url)
  
  # initialize lookup table
  lookup_tbl<-data.frame(CODE_TYPE=as.character(),
                         CODE_TYPE_CDM=as.character(),
                         CODE_SUBTYPE=as.character(),
                         CODE=as.character(),
                         CODE_GRP=as.character(),
                         stringsAsFactors=F)
  
  # main code body for parsing json file
  vs_name_list<-names(vs_file)
  vs_name_dist<-stringdist(tolower(vs_name_str),vs_name_list, method="jw")
  vs_name_match<-vs_name_list[which.min(vs_name_dist)]
  vs<-vs_file[[vs_name_match]]
  for(cd_type in names(vs)){
    vs_sub<-vs[[cd_type]]
    for(cd_subtype in names(vs_sub)){
      # skip if empty
      if(length(vs_sub[[cd_subtype]])==0) next
      # subtype of icd9-cm, icd10-cm codes are: lev0, lev1, lev2
      if(cd_subtype %in% c('lev0','lev1','lev2')){
        cd_lst<-vs_sub[[cd_subtype]]
        # subtypes of hcpcs, icd10-pcs codes are: range, list
      }else if(cd_subtype=="range"){
        rg2vec<-lapply(vs_sub[[cd_subtype]],
                       function(x)seq(as.numeric(strsplit(x,"-")[[1]][1]),as.numeric(strsplit(x,"-")[[1]][2])))
        cd_lst<-unlist(rg2vec)
        cd_subtype<-"exact"
      }else if(cd_subtype %in% c("icd10-pcs","exact")){
        cd_lst<-vs_sub[[cd_subtype]]
      }
      lookup_tbl %<>%
        bind_rows(data.frame(CODE_TYPE=cd_type,
                             CODE_TYPE_CDM=cdm_code_type_map(cd_type),
                             CODE_SUBTYPE=cd_subtype,
                             CODE=as.character(cd_lst),
                             CODE_GRP=vs_name_match,
                             stringsAsFactors = F))
    }
  }
  # return data.frame
  return(lookup_tbl)
}

load_valueset.ecqm<-function(vs_url = "",vs_name_str = ""){
  # load valueset in json
  vs_file_type<-gsub(".*\\.","=",vs_url)
  vs_file<-jsonlite::fromJSON(vs_url)
  
  # initialize lookup table
  lookup_tbl<-data.frame(CODE_TYPE=as.character(),
                         CODE_TYPE_CDM=as.character(),
                         CODE_SUBTYPE=as.character(),
                         CODE=as.character(),
                         CODE_GRP=as.character(),
                         stringsAsFactors=F)
  
  # main code body for parsing json file
  vs_name_list<-names(vs_file)
  vs_name_dist<-stringdist(tolower(vs_name_str),vs_name_list, method="jw")
  vs_name_match<-vs_name_list[which.min(vs_name_dist)]
  vs<-vs_file[[vs_name_match]]
  for(cd_type in names(vs[["codelist"]])){
    cd_lst<-vs[["codelist"]][[cd_type]][,1]
    lookup_tbl %<>%
      bind_rows(data.frame(CODE_TYPE=cd_type,
                           CODE_TYPE_CDM=cdm_code_type_map(cd_type),
                           CODE_SUBTYPE="exact",
                           CODE=as.character(cd_lst),
                           CODE_GRP=vs_name_match,
                           stringsAsFactors = F))
    }
  # return data.frame
  return(lookup_tbl)
}

load_valueset<-function(vs_template = "curated",
                        vs_url = "",
                        vs_name_str = "",
                        dry_run = TRUE,
                        conn=NULL,
                        write_to_schema = "PUBLIC",
                        write_to_tbl = "TEMP",
                        overwrite=TRUE){
  vs_load_func<-get(paste0("load_valueset.",vs_template))
  lookup_tbl<-vs_load_func(vs_url=vs_url,vs_name_str=vs_name_str)
  
  # run query
  if(dry_run==TRUE){
    return(lookup_tbl)
  }else{
    if(is.null(conn)){
      stop("connection needs to be specified!")
    }else{
      # write valueset table to target db
      DBI::dbWriteTable(conn,
                        SQL(paste0(write_to_schema,".",write_to_tbl)),
                        lookup_tbl,
                        overwrite=overwrite,
                        append=!overwrite)
    }
  }
}

cdm_search_key<-function(cdm_tbl){
  if(cdm_tbl=="DIAGNOSIS"){
    return(c('DX','DX_TYPE'))
  }else if(cdm_tbl=="CONDITION"){
    return(c('CONDITION','CONDITION_TYPE'))
  }else if(cdm_tbl=="PROCEDURES"){
    return(c('PX','PX_TYPE'))
  }else if(cdm_tbl=="LAB_RESULT_CM"){
    return('LAB_LOINC')
  }else if(cdm_tbl=="OBS_CLIN"){
    return(c('OBSCLN_CODE','OBSCLN_TYPE'))
  }else if(cdm_tbl=="PRESCRIBING"){
    return('RXNORM_CUI')
  }else if(cdm_tbl=="DISPENSING"){
    return('NDC')
  }else if(cdm_tbl=="MED_ADMIN"){
    return(c('MEDADMIN_CODE','MEDADMIN_TYPE'))
  }else{
    stop("target CDM table is not searchable!")
  }
}

collect_cdm<-function(conn,
                      patset = "PAT_ELIG",
                      key_col_pat = c("PATID"),
                      idx_dt_col = "INDEX_DATE",
                      cdm_schema = "CMS_PCORNET_CDM",
                      key_col_schema=c("PATID"),
                      cdm_tbl = "DIAGNOSIS",
                      ep_dt_col = "EVENT_DATE",
                      vs_template = "curated",
                      vs_url = "",
                      vs_name_str = "",
                      col_out = NULL,
                      write_back = FALSE,
                      write_to_schema = "PUBLIC",
                      write_to_tbl = "TEMP",
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
                           paste0(" from ", write_to_schema, ".", patset, " p"),
                           paste0(" join ", cdm_schema, ".", cdm_tbl," f"),
                           " on ",match_key)))
  
  # load valueset as a data.frame
  lookup_tbl<-load_valueset(vs_url = vs_url,
                            vs_template = vs_template,
                            vs_name_str = vs_name_str)
  
  # specify "where" clause based on valuset type
  search_key<-cdm_search_key(cdm_tbl)
  dx_ind<-as.numeric(grepl("DX",search_key[1]))
  cdtype_scan<-lookup_tbl %>% select(CODE_TYPE) %>% unique %>% unlist
  cdsubtype_scan<-lookup_tbl %>% select(CODE_SUBTYPE) %>% unique %>% unlist
  sql_where<-list()
  for(type in cdtype_scan){
    sql_where_sub<-""
    for(subtype in cdsubtype_scan){
      search_val<-lookup_tbl %>%
        filter(CODE_TYPE==type & CODE_SUBTYPE==subtype) %>%
        select(CODE, CODE_TYPE_CDM) %>% unique
      if(nrow(search_val)==0){
        next
      }else{
        search_key_val=search_val$CODE_TYPE_CDM[1]
      }
      if(subtype=="lev0"){
        sql_where_sub<-paste0(sql_where_sub," or split_part(f.",search_key[1],",'.',1) in ", 
                              paste0("('",paste(search_val$CODE,collapse="','"),"')"))
      }else if(subtype %in% c("lev1","lev2")){
        lev_i<-as.numeric(gsub("lev","",subtype))
        sql_where_sub<-paste0(sql_where_sub," or substr(f.",search_key[1],",1,",3+dx_ind+lev_i,") in ", 
                              paste0("('",paste(search_val$CODE,collapse="','"),"')"))
      }else{
        sql_where_sub<-paste0(sql_where_sub,"f.",search_key[1]," in ",
                              paste0("('",paste(search_val$CODE,collapse="','"),"')"))
      }
    }
    sql_where[[type]]<-gsub("^( or)+","",sql_where_sub)
    if(length(search_key)>1){
      sql_where[[type]]<-paste0(sql_where[[type]]," and f.",search_key[2],"='",search_key_val,"'")
    }
  }
  sql<-paste0(sql," where ",paste0("(",paste(unlist(sql_where),collapse=") or ("),")"))
  
  # run query
  if(dry_run==TRUE){
    return(sql)
  }else{
    if(write_back==TRUE){
      DBI::dbSendQuery(conn,
                       paste0("create or replace table ",
                              write_to_schema,".",write_tbl," as ",sql))
    }else{
      DBI::dbGetQuery(conn,sql)
    }
  }
}
