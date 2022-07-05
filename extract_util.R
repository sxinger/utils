
#######################################################
## project-specific functions for HERON data request ##
#######################################################
#Note: all following functions are specific to projects with dependencies

#purpose: AKI cohort extraction
#dependency: sql scripts of difference types
extract_cohort<-function(conn,
                         cdm_db_name,
                         cdm_db_schema,
                         start_date="2010-01-01",
                         end_date="2018-12-31",
                         verb=T){
  
  #check if DBMS type is currently supported
  if(!(attr(conn,"DBMS_type") %in% c("Oracle","tSQL","PostgreSQL"))){
    stop("DBMS_type=",attr(conn,"DBMS_type"),"is not currently supported \n(should be one of 'Oracle','tSQL','PostgreSQL', case-sensitive)")
  }
  
  #execute(write) the following sql snippets according to the specified order
  statements<-paste0(
    paste0("./src/",attr(conn,"DBMS_type")),
    c("/cohort_initial.sql",
      "/cohort_all_SCr.sql",
      "/cohort_enc_SCr.sql",
      "/cohort_baseline_SCr.sql",
      "/cohort_exclude.sql",
      "/cohort_eligib.sql",
      "/cohort_AKI_staging.sql",
      "/cohort_final.sql")
  )
  
  execute_batch_sql(conn,statements,verb,
                    cdm_db_name=cdm_db_name,
                    cdm_db_schema=cdm_db_schema,
                    start_date=start_date,
                    end_date=end_date)
  
  #collect attrition info
  sql<-parse_sql(paste0("./src/",attr(conn,"DBMS_type"),"/consort_diagram.sql"))
  attrition<-execute_single_sql(conn,
                                statement=sql$statement,
                                write=(sql$action=="write"),
                                table_name=toupper(sql$tbl_out))
  
  #read Table1
  tbl1<-parse_sql(statements[length(statements)])$tbl_out
  aki_enc<-dbGetQuery(conn,paste("select * from",tbl1))
  
  #clean out intermediate tables
  for(i in 1:(length(statements)-1)){
    parse_out<-parse_sql(statements[i])
    if(parse_out$action=="write"){
      drop_tbl(conn,toupper(parse_out$tbl_out))
    }else{
      warning("no temporary table was created by this statment!")
    }
    if(verb){
      cat("temp table",toupper(parse_out$tbl_out),"dropped. \n")
    }
  }
  
  #output
  out<-list(aki_enc=aki_enc,
            attrition=attrition)
  
  return(out)
}

#purpose: AKI Cohort consort diagram
#dependency: AKI attritrion table
consort_diag<-function(consort_tbl){
  require_libraries("diagram")
  tbl<-data.frame(CNT_TYPE=c("Initial",
                             "Has_at_least_2_SCr",
                             "Initial_GFR_below_15",
                             "RRT_within_48hr",
                             "Burn_patients",
                             "Pre_ESRD",
                             "Pre_RRT",
                             "Total",
                             "nonAKI",
                             "AKI1",
                             "nonAKI_to_AKI2",
                             "AKI1_to_AKI2",
                             "nonAKI_to_AKI3",
                             "nonAKI_to_AKI2_to_AKI3",
                             "AKI1_to_AKI2_to_AKI3"),
                  label_txt=c("Inpatient visit with LOS >= 2\nand of age >= 18",
                              "Has at least 2 SCr record",
                              "Excluded: Initial eGFR below 15",
                              "Excluded: RRT with 48 hours since \nadmission",
                              "Excluded: Burn Patients",
                              "Excluded: Pre-existance of \nESRD",
                              "Excluded: Pre-existance of \ndialysis and renal transplantation",
                              "Total eligible encounters",
                              "Non-AKI",
                              "AKI1",
                              "AKI2",
                              "AKI1 to AKI2",
                              "AKI3",
                              "AKI2 to AKI3",
                              "AKI1 to AKI2 to AKI3"),
                  stringsAsFactors=F) %>%
    left_join(consort_tbl, by="CNT_TYPE") %>%
    replace_na(list(ENC_CNT=0)) %>%
    mutate(cnt_ref=ifelse(CNT_TYPE %in% c("Initial","Has_at_least_1_SCr","Total"),ENC_CNT,NA)) %>%
    fill(cnt_ref,.direction="down") %>%
    mutate(cnt_ref=ifelse(CNT_TYPE=="Has_at_least_1_SCr",lag(cnt_ref,n=1L),cnt_ref)) %>%
    mutate(ENC_PROP=round(ENC_CNT/cnt_ref,4)) %>%
    mutate(label_val=paste0("(",ENC_CNT,",",ENC_PROP*100,"%)")) %>%
    mutate(label=paste(label_txt,"\n",label_val)) %>%
    mutate(node_id=c(2,5,7,9,10,12,13,17,18,22,23,25,24,26,28))
  
  #prepare canvas
  par(mfrow=c(1,1))
  par(mar=c(0,0,0,0))
  openplotmat()
  
  ##number of elements per row
  elpos<-coordinates(rep(3,10))
  fromto<-matrix(ncol=2,byrow=T,
                 c(2,5,
                   5,8,
                   8,7,
                   8,9,
                   8,11,
                   11,10,
                   11,12,
                   11,14,
                   14,13,
                   14,17,
                   17,18,
                   17,20,
                   20,19,
                   20,21,
                   19,22,
                   20,23,
                   21,24,
                   22,25,
                   23,26,
                   25,28
                 ))
  ##draw arrows
  arrpos <- matrix(ncol = 2, nrow = nrow(fromto))
  for (i in 1:nrow(fromto)){
    arrpos[i, ] <- straightarrow (to = elpos[fromto[i, 2], ],
                                  from = elpos[fromto[i, 1], ],
                                  lwd = 1, arr.pos = 0.6, arr.length = 0.3)
  }
  
  ##draw nodes
  for(i in 1:nrow(tbl)){
    textrect(elpos[tbl$node_id[i],],
             radx=0.15,
             rady=0.05,
             lab=tbl$label[i],
             font=4,
             cex=0.7)
  }
}



## render report
render_report<-function(which_report="./report/AKI_CDM_EXT_VALID_p1_QA.Rmd",
                        DBMS_type,driver_type,remote_CDM=F,
                        start_date,end_date=as.character(Sys.Date())){
  
  # to avoid <Error in unlockBinding("params", <environment>) : no binding for "params">
  # a hack to trick r thinking it's in interactive environment --not work!
  # unlockBinding('interactive',as.environment('package:base'))
  # assign('interactive',function() TRUE,envir=as.environment('package:base'))
  
  rmarkdown::render(input=which_report,
                    params=list(DBMS_type=DBMS_type,
                                driver_type=driver_type,
                                remote_CDM=remote_CDM,
                                start_date=start_date,
                                end_date=end_date),
                    output_dir="./output/",
                    knit_root_dir="../")
}


## compress dataframe into a condensed format
compress_df<-function(dat,tbl=c("DEMO","VITAL","LAB","DX","PX","MED","DRG"),save=F){
  if(tbl=="DEMO"){
    tbl_zip<-dat %>% 
      filter(key %in% c("AGE","HISPANIC","RACE","SEX")) 
    
    idx_map<-tbl_zip %>% dplyr::select(key) %>%
      mutate(idx=paste0("demo",dense_rank(key))) %>% 
      unique %>% arrange(idx)
    
    tbl_zip %<>%
      spread(key,value,fill=0) %>% #impute 0 for alignment
      unite("fstr",c("AGE","HISPANIC","RACE","SEX"),sep="_")
  }else if(tbl=="VITAL"){
    tbl_zip<-dat %>%
      filter(key %in% c("HT","WT","BMI",
                        "BP_SYSTOLIC","BP_DIASTOLIC",
                        "SMOKING","TOBACCO","TOBACCO_TYPE")) %>%
      mutate(key=recode(key,
                        HT="1HT",
                        WT="2WT",
                        BMI="3BMI",
                        SMOKING="4SMOKING",
                        TOBACCO="5TOBACCO",
                        TOBACCO_TYPE="6TOBACCO_TYPE",
                        BP_SYSTOLIC="7BP_SYSTOLIC",
                        BP_DIASTOLIC="8BP_DIASTOLIC")) %>%
      mutate(add_time=difftime(timestamp,format(timestamp,"%Y-%m-%d"),units="mins")) %>%
      mutate(dsa=dsa+round(as.numeric(add_time)/(24*60),2)) %>%
      arrange(key,dsa) 
    
    idx_map<-tbl_zip %>% dplyr::select(key) %>%
      mutate(idx=paste0("vital",dense_rank(key))) %>% 
      unique %>% arrange(idx)
    
    tbl_zip %<>% unique %>%
      unite("val_date",c("value","dsa"),sep=",") %>%
      group_by(ENCOUNTERID,key) %>%
      dplyr::summarize(fstr=paste(val_date,collapse=";")) %>%
      ungroup %>%
      spread(key,fstr,fill=0) %>% #impute 0 for alignment
      unite("fstr",c("1HT","2WT","3BMI",
                     "4SMOKING","5TOBACCO","6TOBACCO_TYPE",
                     "7BP_SYSTOLIC","8BP_DIASTOLIC"),sep="_")
  }else if(tbl=="LAB"){
    tbl_zip<-dat %>%
      mutate(idx=paste0("lab",dense_rank(key)))
    
    idx_map<-tbl_zip %>% dplyr::select(key,idx) %>%
      unique %>% arrange(idx)
    
    tbl_zip %<>%
      arrange(ENCOUNTERID,idx,dsa) %>%
      unite("val_unit_date",c("value","unit","dsa"),sep=",") %>%
      group_by(ENCOUNTERID,idx) %>%
      dplyr::summarize(fstr=paste(val_unit_date,collapse=";")) %>%
      ungroup %>%
      unite("fstr2",c("idx","fstr"),sep=":") %>%
      group_by(ENCOUNTERID) %>%
      dplyr::summarize(fstr=paste(fstr2,collapse="_")) %>%
      ungroup
  }else if(tbl=="DRG"){
    tbl_zip<-dat %>%
      mutate(idx=paste0("dx",dense_rank(key2))) 
    
    idx_map<-tbl_zip %>% dplyr::select(key2,idx) %>%
      unique %>% arrange(idx) %>% dplyr::rename(key=key2)
    
    tbl_zip %<>%
      group_by(ENCOUNTERID,key1,idx) %>%
      dplyr::summarize(dsa=paste(dsa,collapse=",")) %>%
      ungroup %>%
      unite("fstr",c("idx","dsa"),sep=":") %>%
      group_by(ENCOUNTERID,key1) %>%
      dplyr::summarize(fstr=paste(fstr,collapse="_")) %>%
      ungroup %>%
      spread(key1,fstr,fill=0) %>%
      unite("fstr",c("ADMIT_DRG","COMMORB_DRG"),sep="|") %>%
      unique
  }else if(tbl=="DX"){
    tbl_zip<-dat %>%
      group_by(ENCOUNTERID,key) %>%
      arrange(dsa) %>%
      dplyr::summarize(dsa=paste(dsa,collapse=",")) %>%
      ungroup %>%
      mutate(idx=paste0("ccs",key))
    
    idx_map<-tbl_zip %>% dplyr::select(key,idx) %>%
      unique %>% arrange(key)
    
    tbl_zip %<>%
      unite("fstr",c("idx","dsa"),sep=":") %>%
      group_by(ENCOUNTERID) %>%
      dplyr::summarize(fstr=paste(fstr,collapse="_")) %>%
      ungroup %>% unique
  }else if(tbl=="PX"){
    tbl_zip<-dat %>%
      mutate(idx=paste0("px",dense_rank(key))) 
    
    idx_map<-tbl_zip %>% dplyr::select(key,idx) %>%
      unique %>% arrange(idx)
    
    tbl_zip %<>%
      group_by(ENCOUNTERID,idx) %>%
      arrange(dsa) %>%
      dplyr::summarize(dsa=paste(dsa,collapse=",")) %>%
      ungroup %>%
      unite("fstr",c("idx","dsa"),sep=":") %>%
      group_by(ENCOUNTERID) %>%
      dplyr::summarize(fstr=paste(fstr,collapse="_")) %>%
      ungroup %>% unique
  }else if(tbl=="MED"){
    tbl_zip<-dat %>%
      mutate(idx=paste0("med",dense_rank(key))) 
    
    idx_map<-tbl_zip %>% dplyr::select(key,idx) %>%
      unique %>% arrange(idx)
    
    tbl_zip %<>%
      transform(value=strsplit(value,","),
                dsa=strsplit(dsa,",")) %>%
      unnest %>%
      unite("val_date",c("value","dsa"),sep=",") %>%
      group_by(ENCOUNTERID,idx) %>%
      dplyr::summarize(fstr=paste(val_date,collapse=";")) %>%
      ungroup %>%
      unite("fstr2",c("idx","fstr"),sep=":") %>%
      group_by(ENCOUNTERID) %>%
      dplyr::summarize(fstr=paste(fstr2,collapse="_")) %>%
      ungroup
  }else{
    warning("data elements not considered!")
  }
  if(save){
    save(tbl_zip,file=paste0("./data/",tbl,"_zip.Rdata"))
  }
  
  zip_out<-list(tbl_zip=tbl_zip,idx_map=idx_map)
  return(zip_out)
}







