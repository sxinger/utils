############################################
## utility functions for complex modeling ##
############################################
# require(tidyverse,magrittr)

# require(survival)
coxph_stratified<-function(dt,time_col,status_col,
                           expos_col="", # column of exposure/intervention
                           cov_col=list(), # covariate columns for full model
                           cols_strata=list(),
                           cols_excld=list() # 1-to-1 mapping with cols_strata
                           ){
  # create strata metadata file
  strata<-c()
  for(i in seq_len(cols_strata)){
    col_nm<-cols_strata[i]
    col_i<-dt[,col_nm] %>% unlist
    strata<-rbind(strata,
                  data.frame(
                    val=unique(col_i),
                    var=rep(col_nm,length(col_i)),
                    excld=rep(cols_excld[i],length(col_i))
                    )
                  )
  }
  # stack regression results from stratified models
  result<-c()
  for(i in seq_len(nrow(strata))){
    # curate variable list
    var_filter<-cov_col[grepl(cov_col %in% colnames(dt))]
    var_filter<-var_filter[!grepl(strata$excld[i],var_filter)]
    # form regression formula
    fit_frm<-formula(paste0("Surv(",time_col,",",status_col,") ~ ",
                          paste(c(var_filter,expos_col),collapse = "+")))
    fit_mort_cov<-coxph(fit_frm, 
                        data = dt %>% 
                          filter(.data[[strata$var[i]]]==strata$val[i]))
    # get all coefficients
    fit_summ<-summary(fit_mort_cov)$coefficients
    fit_var<-rownames(fit_summ)
    rownames(fit_summ)<-NULL
    result<-rbind(result,
                  cbind(stratum_var=strata$var[i],
                      stratum_val=strata$val[i],
                      fit_var=fit_var,
                      fit_summ))
  }                            
  return(result)
}

# require(glmnet,islasso)
iptw.lasso<-function(
  data_df, #
  id_col = 'PATID', # primary key
  yc = '', # column name of exposure at center
  yo_vec = c(""), # vector of other relevant exposures
  xo_vec = c(""), # vector of other covariates
  family = 'binomial', # ref to legal values for "glmnet"
  type.measure = "auc" # ref to legal values for "glmnet"
){
  # conversion to matrix
  x<-data.matrix(data_df[,xo_vec])
  
  # loop over yo_vec
  out<-list()
  tw_composit<-c()
  for(yo in seq_along c(yc,yo_vec)){
    out_yo<-list()
    # calculate weight with smoothing
    y<-unlist(data_df[,yo])
    fit_tw<-cv.glmnet(x=x,y=y,family=family,type.measure = type.measure)
    tw<-predict(fit_tw, newx = x, s = "lambda.min", type="response")
    tw_smth<-data.frame(id=id,tw=tw[,1]) %>% 
      mutate(idx=row_number()) %>%
      arrange(tw) %>%  
      mutate(tw_adj=rank(tw)/n()) %>%
      mutate(tw_adj=case_when(y==1 ~ tw_adj,
                              y==0&tw_adj==1 ~ 1-(tw_adj-0.0001),
                              TRUE ~ 1-tw_adj)) %>%
      arrange(idx) %>% 
      select(id,tw,tw_adj)
    colnames(tw_smth)<-c("id",yo,paste0(yo,"_adj"))

    # decompose propensity score
    fitx<-islasso(y~x,data=data_df,lambda=fit_tw$lambda.min)
    out1<-data.frame(summary(fitx,pval=0.1)$coefficients,stringsAsFactors = F) %>%
      rownames_to_column(var="varx")
    
    out_yo<-list(
      ps_tw = tw_smth,
      ps_decomp = out1
    )
    if(yo != yc){
      # mediator model
      xx<-data.matrix(data_df[,c(var_ps,yc)])
      fit_xx<-cv.glmnet(x=xx,y=y,family=family,type.measure = type.measure)
      out_xx<-islasso(y~xx,data=data_df,lambda=fit_xx$lambda.min)
      out2<-data.frame(summary(out_xx)$coefficients,stringsAsFactors = F) %>%
            rownames_to_column(var="varx")
      out_yo[['ps_medi']] = out2
    }

    # stack results
    out[[yo]]<-out_yo
    tw_composit %<>% inner_join(out_yo$ps_tw,by="id")
  }
  
  # calculate composit weight
  tw_composit %<>%
    mutate(tw_comp = eval(parse(paste(yo,collapse="*"))),
           idx=row_number()) %>%
    arrange(tw_comp) %>% 
    mutate(tw_comp_adj=rank(tw_comp)/n()) %>%
    mutate(tw_comp_adj=case_when(y==1 ~ tw_comp_adj,
                                 y==0&tw_comp_adj==1 ~ 1-(tw_comp_adj-0.0001),
                                 TRUE ~ 1-tw_comp_adj)) %>%
    arrange(idx) %>% select(id,tw_comp,tw_comp_adj)

    # stack result
    out[['composit']]<-list(
      ps_tw = tw_composit
    )

    return(out)
}

knockoff_grpreg<-function(){

}


boruta_grpreg<-function(){

}

# #TODO
# explain_model<-function(){
#   
# }
# 
# #TODO
# adjMMD<-function(){
#   
# }
