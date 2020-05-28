#######################################
## utility functions for database I/O##
#######################################
##connect to HERON database using ODB connector
#fill up the config file with your oracle username and password 
connect_to_oci<-function(config_file){
  conn<-DBI::dbConnect(ROracle::Oracle(),
                       config_file$username,
                       config_file$password,
                       paste0("//",file.path(config_file$access,config_file$sid)))
  return(conn)
}

##connect to HERON database using JDBC connector
#fill up the config file with your oracle username and password 
connect_to_jdbc<-function(config_file){
  drv<-JDBC(driverClass="oracle.jdbc.OracleDriver",
            classPath="./inst/ojdbc6.jar")
  url <- paste0("jdbc:oracle:thin:@localhost:1521:", 
                config_file$sid)
  conn <- RJDBC::dbConnect(drv=drv,
                           url=url, 
                           user=config_file$username, 
                           password=config_file$password)
  return(conn)
}


chunk_load<-function(conn,dataset="",chunk_size=1000,verb=T){
  dat<-c()
  i<-0
  error<-FALSE
  row_remain<-Inf
  while(!error&row_remain>0){
    #try load
    dat_add<-try(dbGetQuery(conn,
                            paste("select * from (",
                                  "select m.*, rownum r from",dataset," m)",
                                  "where r >= ",i+1,"and r < ",i+chunk_size)),
                 silent=T)
    
    #check unexpected errors (e.g. connection error)
    error<-grepl("error",tolower(dat_add))[1]
    
    #check remaining rows
    row_remain<-nrow(dat_add)
    
    #attach rows
    dat %<>% bind_rows(dat_add)
    
    #report progress
    if(verb){
      cat("row",i+1,"to","row",i+chunk_size,"loaded.\n") 
    }
    
    #loop updates
    i<-i+chunk_size
  }
  
  return(dat)
}


## parse Oracle sql lines
parse_sql<-function(file_path,...){
  param_val<-list(...)
  
  #read file
  con<-file(file_path,"r")
  
  #initialize string
  sql_string <- ""
  
  #intialize result holder
  params_ind<-FALSE
  tbl_out<-NULL
  action<-NULL
  
  while (TRUE){
    #parse the first line
    line <- readLines(con, n = 1)
    #check for endings
    if (length(line)==0) break
    #collect overhead info
    if(grepl("^(/\\*out)",line)){
      #output table name
      tbl_out<-trimws(gsub("(/\\*out\\:\\s)","",line),"both")
    }else if(grepl("^(/\\*action)",line)){
      #"write" or "query"(fetch) the output table
      action<-trimws(gsub("(/\\*action\\:\\s)","",line),"both")
    }else if(grepl("^(/\\*params)",line)){
      params_ind<-TRUE
      #breakdown global parameters
      params<-gsub(",","",strsplit(trimws(gsub("(/\\*params\\:\\s)","",line),"both")," ")[[1]])
      params_symbol<-params
      #normalize the parameter names
      params<-gsub("&&","",params) 
    }
    #remove the first line
    line<-gsub("\\t", " ", line)
    #translate comment symbol '--'
    if(grepl("--",line) == TRUE){
      line <- paste(sub("--","/*",line),"*/")
    }
    #attach new line
    if(!grepl("^(/\\*)",line)){
      sql_string <- paste(sql_string, line)
    }
  }
  close(con)
  
  #update parameters as needed
  if(params_ind){
    #align param_val with params
    params_miss<-params[!(params %in% names(param_val))]
    for(j in seq_along(params_miss)){
      param_val[params_miss[j]]<-list(NULL)
    }
    param_val<-param_val[which(names(param_val) %in% params)]
    param_val<-param_val[order(names(param_val))]
    params_symbol<-params_symbol[order(params)]
    params<-params[order(params)]
    
    #substitube params_symbol by param_val
    for(i in seq_along(params)){
      sql_string<-gsub(params_symbol[i],
                       ifelse(is.null(param_val[[i]])," ",
                              ifelse(params[i]=="cdm_db_link",
                                     paste0("@",param_val[[i]]),
                                     ifelse(params[i] %in% c("start_date","end_date"),
                                            paste0("'",param_val[[i]],"'"),
                                            param_val[[i]]))),
                       sql_string)
    }
  }
  #clean up excessive "[ ]." or "[@" in tSQL when substitute value is NULL
  sql_string<-gsub("\\[\\ ]\\.","",sql_string)
  sql_string<-gsub("\\[@","[",sql_string)
  
  out<-list(tbl_out=tbl_out,
            action=action,
            statement=sql_string)
  
  return(out)
}


## execute single sql snippet (JDBC)
execute_single_sql<-function(conn,statement,write,table_name,verb){
  if(write){
    try_tbl<-try(dbGetQuery(conn,
                            paste("select * from",table_name,"where 1=0")),
                 silent=T)
    if(is.null(attr(try_tbl,"condition"))){
      dbSendUpdate(conn,paste("drop table",table_name)) #in case there exists same table name
    }
    
    dbSendUpdate(conn,statement)
    
  }else{
    dat<-dbGetQuery(conn,statement)
    return(dat)
  }
  
  if(verb){
    cat(ifelse(write,"write","query"),"table",table_name,"on grouse USER tablespace.\n")
  }
}


## execute multiple sql snippets
#---statements have to be in correct logical order
execute_batch_sql<-function(conn,statements,verb,benchmark=F,...){
  if(benchmark){
    bm<-c()
    start<-Sys.time()
  }
  
  for(i in seq_along(statements)){
    sql<-parse_sql(file_path=statements[i],...)
    execute_single_sql(conn,
                       statement=sql$statement,
                       write=(sql$action=="write"),
                       table_name=toupper(sql$tbl_out),
                       verb=verb)
    if(verb){
      cat(statements[i],"has been executed and table",
          toupper(sql$tbl_out),"was created.\n")
    }
    
    if(benchmark){
      lapse<-Sys.time()-start
      bm<-c(bm,paste0(round(lapse,1),units(lapse)))
    }
  }
  
  if(benchmark){
    return(bm)
  }
}


## clean up intermediate tables
drop_tbl<-function(conn,table_name){
  dbSendUpdate(conn,paste("drop table",table_name,"purge"))
}


