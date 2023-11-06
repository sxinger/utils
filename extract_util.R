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

load_valueset.ncbo<-function(
  vs_url = "",
  vs_name_str = ""
){
  # load valueset in json
  vs_file_type<-gsub(".*\\.","=",vs_url)
  vs_file<-jsonlite::fromJSON(vs_url)
  
  # initialize lookup table
  lookup_tbl<-data.frame(
    CODE_TYPE=as.character(),
    CODE_TYPE_CDM=as.character(),
    CODE_SUBTYPE=as.character(),
    CODE=as.character(),
    CODE_GRP=as.character(),
    stringsAsFactors=F
    )
  
  # main code body for parsing json file
  vs_name_list<-names(vs_file)
  vs_name_dist<-stringdist(tolower(vs_name_str),vs_name_list, method="jw")
  vs_name_match<-vs_name_list[which.min(vs_name_dist)]
  vs<-vs_file[[vs_name_match]]
  for(cd_type_idx in seq_along(vs[["code_type"]])){
    # skip if empty
    if(length(vs[["code_list"]])==0) next
    lookup_tbl %<>%
      bind_rows(data.frame(CODE_TYPE=vs[["code_type"]][[cd_type_idx]],
                           CODE_TYPE_CDM=cdm_code_type_map(vs[["code_type"]][cd_type_idx]),
                           CODE_SUBTYPE="exact",
                           CODE=vs[["code_list"]][[cd_type_idx]][["code"]],
                           CODE_GRP=vs_name_match,
                           stringsAsFactors = F))
  }
  # return data.frame   
  return(lookup_tbl)
}

load_valueset.rxnav<-function(
  vs_url = "",
  vs_name_str = ""
){
  # load valueset in json
  vs_file_type<-gsub(".*\\.","=",vs_url)
  vs_file<-jsonlite::fromJSON(vs_url)
  
  # initialize lookup table
  lookup_tbl<-data.frame(RXCUI=as.character(),
                         LABEL=as.character(),
                         NDC=as.character(),
                         stringsAsFactors=F)
  
  # main code body for parsing json file
  vs_name_list<-names(vs_file)
  vs_name_dist<-stringdist(tolower(vs_name_str),vs_name_list, method="jw")
  vs_name_match<-vs_name_list[which.min(vs_name_dist)]
  vs<-vs_file[[vs_name_match]] %>% 
    filter(ndc!="character(0)") %>%
    unnest_wider(ndc,names_sep="_") %>%
    gather(ndc_idx,ndc,-rxcui,-label) %>%
    filter(!is.na(ndc)) %>% select(-ndc_idx)

  # return data.frame
  colnames(vs)<-toupper(colnames(vs))
  lookup_tbl %<>% bind_rows(data.frame(vs))
  return(lookup_tbl)
}

load_valueset.curated<-function(
  vs_url = "",
  vs_name_str = "",
  add_meta = FALSE
){
  # load valueset in json
  vs_file_type<-gsub(".*\\.","=",vs_url)
  vs_file<-jsonlite::fromJSON(vs_url)
  
  # initialize lookup table
  lookup_tbl<-data.frame(
    CODE_TYPE=as.character(),
    CODE_TYPE_CDM=as.character(),
    CODE=as.character(),
    CODE_GRP=as.character(),
    stringsAsFactors=F
  )
  if(add_meta){
    meta_tbl<-c()
  }

  # search closest concept key
  vs_name_list<-names(vs_file)
  if(vs_name_str != ""){
    vs_name_dist<-stringdist(tolower(vs_name_str),vs_name_list, method="jw")
    vs_name_match<-vs_name_list[which.min(vs_name_dist)]
  }else{
    vs_name_match<-vs_name_list
  }

  for(key in vs_name_match){
    vs<-vs_file[[key]]
    if(add_meta){
      # metadata
      meta_tbl %<>% 
        bind_rows(
          cbind(
            CODE_GRP=key,
            as.data.frame(do.call(cbind, vs[["meta"]]))
          )
        )
    }

    # valuesets
    for(cd_type in names(vs)){
      if(cd_type == "meta") next
      # exact list without decimal points
      cd_lst<-vs[[cd_type]] 
      # stack up
      lookup_tbl %<>%
        bind_rows(data.frame(
          CODE_TYPE=cd_type,
          CODE_TYPE_CDM=cdm_code_type_map(cd_type),
          CODE=as.character(cd_lst),
          CODE_GRP=key,
          stringsAsFactors = F
        ))
    }
  }
  # return data.frame
  out<-lookup_tbl %>% 
    inner_join(meta_tbl,by="CODE_GRP")
  return(out)
}

load_valueset.ecqm<-function(
  vs_url = "",
  vs_name_str = ""
){
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

load_valueset.vsac<-function(
  vs_url = "",
  vs_name_str = ""
){
  # load valueset in json
  vs_file_type<-gsub(".*\\.","=",vs_url)
  lookup_tbl<-jsonlite::fromJSON(vs_url) %>%
    unnest(vs) %>%
    rename(
      "CODEGRP" = "oid" ,
      "CODEGRP_LABEL" = "displayName",
      "CODE" = "@code",
      "CODE_LABEL" = "@displayName",
      "CODE_TYPE" = "@codeSystemName"
    ) %>%
    mutate(CODE_TYPE_CDM = unlist(lapply(CODE_TYPE, function(x) cdm_code_type_map(x)))) %>%
    select(CODE_TYPE,CODE_TYPE_CDM,CODE,CODE_LABEL,CODEGRP,CODEGRP_LABEL) 

  # return data.frame
  return(lookup_tbl)
}

load_valueset<-function(
  vs_template = c("curated",
                  "ecqm",
                  "ncbo",
                  "rxnav",
                  "vsac"),
  vs_url = "",
  vs_name_str = "",
  add_meta = TRUE,
  dry_run = TRUE,
  conn=NULL,
  write_to_schema = "PUBLIC",
  write_to_tbl = "TEMP",
  overwrite=TRUE,
  file_encoding ="latin-1"
){
  vs_load_func<-get(paste0("load_valueset.",vs_template))
  if(vs_template=="curated"){
    lookup_tbl<-vs_load_func(vs_url=vs_url,vs_name_str=vs_name_str)
  }else{
    lookup_tbl<-vs_load_func(vs_url=vs_url,vs_name_str=vs_name_str,add_meta = add_meta)
  }
  
  # run query
  if(dry_run==TRUE){
    return(lookup_tbl)
  }else{
    if(is.null(conn)){
      stop("connection needs to be specified!")
    }else{
      # specify field.types to accommodate long strings
      max_str<-rep("varchar(500)",ncol(lookup_tbl))
      names(max_str) <- names(lookup_tbl)
      if(!overwrite){
        max_str = NULL
      }
      # write valueset table to target db
      DBI::dbWriteTable(
        conn,
        SQL(paste0(write_to_schema,".",write_to_tbl)),
        lookup_tbl,
        overwrite = overwrite,
        append = !overwrite,
        file_encoding = file_encoding,
        field.types = max_str
      )
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

load_mapping.resdac<-function(mp_name=c("SSA_STATE",
                                        "PRVDR_SPCLTY",
                                        "CMS_TYPE_SRVC")){
  if(grepl("SSA_STATE",mp_name)){
    tbl<-read.table("https://resdac.org/sites/datadocumentation.resdac.org/files/State%20Table.txt",
                    sep="=", skip = 2) %>%
      rename(CODE = V1, LABEL = V2) %>%
      mutate(CODE = str_pad(trimws(CODE),2,"left","0"),
             LABEL = trimws(gsub(" \\(.*","",LABEL)))

  }else if(grepl("PRVDR_SPCLTY",mp_name)){
    tbl<-readLines("https://resdac.org/sites/datadocumentation.resdac.org/files/CMS_PRVDR_SPCLTY_TB_rev01242018_0.txt") %>% as.data.frame
    colnames(tbl)<-"line"
    tbl %<>%
      filter(trimws(line) != "") %>% slice(-1) %>%
      separate(line,c("CODE","LABEL"),sep="=",extra="merge",fill="left") %>%
      fill(CODE,.direction = "down") %>%
      group_by(CODE) %>%
      summarise(LABEL = paste(trimws(LABEL),collapse = " "),.groups="drop")
    
  }else if(grepl("CMS_TYPE_SRVC",mp_name)){
    tbl<-read.table("https://resdac.org/sites/datadocumentation.resdac.org/files/CMS%20Type%20of%20Service%20Table.txt",
                         sep="=", skip = 2) %>%
      rename(CODE = V1, LABEL = V2) %>%
      mutate(CODE = trimws(CODE),
             LABEL = trimws(gsub(" \\(.*","",LABEL)))
  }
  return(tbl)
}

load_mapping<-function(mp_template = "resdac",mp_name){
  mp_load_func<-get(paste0("load_mapping.",mp_template))
  mp_tbl<-mp_load_func(mp_name)
  return(mp_tbl)
}