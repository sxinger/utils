
###########################################
## utility functions for sampling frame  ##
###########################################

#exact matching with coarsening, with/without replacement
matched_sample<-function(case_ref,
                         ctrl_pool,
                         id_col,
                         exact=c(),
                         coarse=c(),
                         coarse_range=c(),
                         other_cov=c(),
                         weighting_scheme=c("weighted","hiearchical"),
                         wt=c(1),
                         match_rt=match_rt,
                         replace=FALSE,
                         verbose=TRUE){
  #only keep relavant columns
  case_ref %<>% select(c(id_col,exact,coarse,other_cov))
  ctrl_pool %<>% select(c(id_col,exact,coarse,other_cov))
  
  #--exact matching with coarsening (relative range)
  filter_cond_vec<-c()
  for(cond_i in seq_along(coarse)){
    filter_cond_vec<-c(filter_cond_vec,
                       paste0(coarse[cond_i],".y",">=",coarse[cond_i],".x-",coarse_range[cond_i],"&",
                              coarse[cond_i],".y","<=",coarse[cond_i],".x+",coarse_range[cond_i]))
  }
  filter_cond<-paste0("(",paste(filter_cond_vec,collapse = ")&("),")")
  
  #--apply exact matching conditions
  ctrl_match<-case_ref %>% unique %>%
    dplyr::rename("id"=id_col) %>% #for convenience
    inner_join(ctrl_pool,by=exact)
  
  #--apply relative/coarsening conditions
  if(length(coarse) > 0){
    ctrl_match %<>% filter(eval(parse(text=filter_cond)))
  }
  
  #checkpoint--if no additional matched samples can be found, break out
  if(nrow(ctrl_match)==0){
    stop("no matched sample can be found!")
  }
  
  #--rank multiple matched controls based on similarity metrics
  similarity_cond_vec<-c()
  if(weighting_scheme=="hiearchical"){
    wt<-rank(wt)
    for(cond_i in seq_along(coarse)){
      similarity_cond_vec<-c(similarity_cond_vec,
                             paste0("abs((",coarse[cond_i],".x","-",coarse[cond_i],".y",")*",10^(wt[cond_i]-1),")"))
    }
    similarity_formula<-paste0("rank(",paste(similarity_cond_vec,collapse = "+"),",ties.method=\"random\")")
    
  }else if(weighting_scheme=="weighted"){
    for(cond_i in seq_along(coarse)){
      similarity_cond_vec<-c(similarity_cond_vec,
                             paste0("abs((",coarse[cond_i],".x","-",coarse[cond_i],".y",")*",wt[cond_i],")"))
    }
    similarity_formula<-paste(similarity_cond_vec,collapse = "+")
    
  }else{
    stop("similarity measuring method between case and control is not currently supported!")
  }
  
  ctrl_match %<>%
    group_by(id) %>%
    mutate(similarity_rank=eval(parse(text=similarity_formula))) %>%
    ungroup 
  
  #sample without replacement
  max_rk<-max(ctrl_match$similarity_rank)
  case_n<-length(unique(ctrl_match$id))
  ctrl_match_undup<-c()
  
  match_rd<-1
  while(match_rd<=match_rt){
    #--identify a new batch of matched samples
    ctrl_match_undup_rd<-ctrl_match %>%
      filter(similarity_rank<=1)
    
    if(!replace){
      #--matched samples could be the same over different cases
      ctrl_match_undup_rd %<>%
        distinct(mrn,.keep_all=TRUE)
      
      ctrl_match %<>% 
        anti_join(ctrl_match_undup_rd,by=id_col) %>%
        group_by(id) %>%
        dplyr::mutate(similarity_rank=rank(similarity_rank,ties.method = "random")) %>% #move 2nd matched sample up
        ungroup
      
      #--go over the pool of matched samples until there is no duplicates
      matched_n<-length(unique(ctrl_match_undup_rd$id))
      while(matched_n<case_n&match_rd<=match_rt&match_rd<=max_rk) {
        case_unmatched<-ctrl_match %>%
          anti_join(ctrl_match_undup_rd,by="id") %>%
          filter(similarity_rank > 1) %>% #first matched sample already been picked
          group_by(id) %>%
          dplyr::mutate(similarity_rank=rank(similarity_rank)) %>% #move 2nd matched sample up
          ungroup
        
        #--when no more cases can be matched
        if(nrow(case_unmatched)==0) break
        
        ctrl_match_undup_rd %<>%
          bind_rows(case_unmatched %>%
                      filter(similarity_rank<=1))
        
        #--end-of-inner-loop updates
        matched_n<-length(unique(ctrl_match_undup_rd$id))
        ctrl_match %<>% 
          anti_join(ctrl_match_undup_rd,by=id_col) %>%
          group_by(id) %>%
          dplyr::mutate(similarity_rank=rank(similarity_rank,ties.method = "random")) %>% #move 2nd matched sample up
          ungroup
      }
    }
    
    #--end-of-outer-loop updates
    ctrl_match_undup %<>% 
      bind_rows(ctrl_match_undup_rd %>% mutate(similarity_rank=match_rd))
    
    if(verbose){
      cat("match round:",match_rd,"; ","matched samples:",matched_n,"\n")
    }
    
    match_rd<-match_rd+1
  }
  
  id_ctrl<-paste0(id_col,"_ctrl")
  ctrl_match_undup %<>%
    dplyr::rename(!!sym(id_ctrl) := id_col) %>%
    dplyr::rename(!!sym(id_col) := "id")
  
  
  return(ctrl_match_undup)
}

#one-hot coding is required!
#require (RANN, data.table)
strata_sample<-function(ref_dat, #reference dataset
                        match_dat, #matching dataset
                        keep_col="patient_num",
                        compare_metric=c("age","sex"), #matching criteria
                        update_match_metric=NULL, #copy values of matching sample from ref_dat
                        boots=5,
                        nnk=boots+1, #number of candidate neighbors, recommend:nnk>=boots
                        searchtype=c("standard", "priority", "radius"),
                        replace=F,
                        verb=T){
  #attach row_id
  ref_dat$row_id<-1:nrow(ref_dat)
  match_dat$row_id<-(nrow(ref_dat)+1):(nrow(ref_dat)+nrow(match_dat))
  
  boots_samp<-c()
  for(k in 1:boots){
    start_k<-Sys.time()
    if(verb){
      cat("bootstrapped sample:",k,"\n")
    }

    # identify k-nearest-neighbour
    sample_pos<-ref_dat$row_id
    sample_neg<-nn2(match_dat[,compare_metric],
                    ref_dat[,compare_metric],
                    k=nnk,
                    searchtype=searchtype)$nn.idx[,sample(seq_len(nnk),1)] #inject randomness
    sample_neg2<-match_dat[sample_neg,]$row_id
    idx_map<-data.frame(pos_idx=sample_pos,neg_idx=sample_neg2)
    
    # update certain field by copying matched value from ref_dat
    ref_match<-match_dat[sample_neg,c("row_id",update_match_metric)]
    if(!is.null(update_match_metric)){
      ref_match<-ref_match %>%
        inner_join(idx_map,by=c("row_id"="neg_idx")) %>%
        unique
      ref_match<-ref_match %>%
        select(-all_of(update_match_metric)) %>%
        inner_join(ref_dat %>% select(all_of(c("row_id",update_match_metric))),
                  by=c("pos_idx"="row_id")) %>%
        select(-pos_idx)
    }

    # reconstruct stratified samples
    col_sel<-c("row_id",compare_metric,keep_col)
    if(!is.null(update_match_metric)&&update_match_metric %in% col_sel){
      col_sel<-col_sel[!col_sel %in% update_match_metric]
    }
    sample_reconst<-ref_match %>%
      inner_join(match_dat %>% select(all_of(col_sel)),
                by="row_id") %>%
      select(-row_id) %>%
      mutate(boots_rnd=k)

    # stack bootstrapped sample
    boots_samp<-rbind(boots_samp,sample_reconst)
    if(verb){
      cat("Finish bootstrap sample ",k," in ",
          Sys.time()-start_k,units(Sys.time()-start_k),"\n")
    }
    
    # if without replacement, remove the current sample from sample pool
    if(!replace){
      match_dat<-match_dat[!(match_dat$row_id %in% unique(sample_neg2)),]
    }
  }
  return(boots_samp)
}
