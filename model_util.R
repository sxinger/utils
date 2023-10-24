############################################
## utility functions for complex modeling ##
############################################

coxph_stratified<-function(
  dt,time_col,status_col,
  expos_col="", # column of exposure/intervention
  cov_col=list(), # covariate columns for full model
  cols_strata=list(),
  cols_excld=list() # 1-to-1 mapping with cols_strata
){
  # require(tidyverse,magrittr,survival)

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

ipw.lasso<-function(
  data_df, # data.frame including id_col, yc, yo_vec, xo_vec
  id_col = 'PATID', # primary key
  yc = 'TRT', # column name of exposure at center,
  yo_vec = c(""), # vector of other mediators likely on the pathway
  xo_vec = c(""), # vector of other covariates
  ycs = NULL, # column name of censoring indicator, if informative censoring needs to be controled
  family = 'binomial', # ref to legal values for "glmnet"
  type.measure = "class", # ref to legal values for "glmnet"
  verb = TRUE #verbose
){
  # require(tidyverse,magrittr,glmnet,islasso,scales)
  # loop over yo_vec
  out<-list()
  for(yo_i in seq_along(c(yc,ycs,yo_vec))){
    yo<-c(yc,ycs,yo_vec)[yo_i]
    ################################################################
    if(verb) print(sprintf("start propensity analysis for: %s",yo))
    ################################################################
    out_yo<-list()
    # calculate weight with smoothing
    y<-unlist(data_df[,yo])
    # form covariate matrix
    if(yo==yc){
      x<-data.matrix(data_df[,xo_vec])
    }else{
      # include yc on the pathway assuming yo_vec is mediated by yc
      x<-data.matrix(data_df[,c(yc,xo_vec)])
    }
    fit_tw<-cv.glmnet(x=x,y=y,family=family,type.measure = type.measure,alpha=1)
    tw<-predict(fit_tw, newx = x, s = "lambda.min", type="response")
    id<-unlist(data_df[,id_col])
    tw_smth<-data.frame(id=id,tw=tw[,1]) %>%  
      mutate(ipw = case_when(y==1&tw==0 ~ NA_real_,
                             y==1 ~ 1/tw,
                             y==0&tw==1 ~ NA_real_,
                             TRUE ~ 1/(1-tw)),
             ow = case_when(y==1&tw==1 ~ NA_real_,
                            y==1 ~ 1 - tw,
                            y==0&tw==0 ~ NA_real_,
                            TRUE ~ tw)) %>%
      mutate(across(c("ipw"), ~replace_na(.x, max(.x, na.rm = TRUE)))) %>% 
      mutate(across(c("ow"), ~replace_na(.x, min(.x, na.rm = TRUE)))) %>% 
      select(id,ipw,ow)
    colnames(tw_smth)<-c("id",yo,paste0(yo,"_ow"))
    ################################################################
    if(verb) print("...finish generating weights.")
    ################################################################   
    # decompose propensity score
    fitx<-islasso(y~x,data=data_df,lambda=fit_tw$lambda.min)
    out1<-data.frame(summary(fitx,pval=0.1)$coefficients,stringsAsFactors = F) %>%
      rownames_to_column(var="varx")
    ################################################################
    if(verb) print(paste("...finish decomposing propensity."))
    ################################################################
    # stack results
    out[[yo]]<-list(
      ps_tw = tw_smth,
      ps_decomp = out1
    )
    if(yo_i==1){
      tw_composit<-tw_smth
    }else{
      tw_composit %<>% inner_join(tw_smth,by="id")
    }
    ################################################################
    # mediator model 
    if(yo != yc){
      xx<-data.matrix(data_df[,c(xo_vec,yc)])
      fit_xx<-cv.glmnet(x=xx,y=y,family=family,type.measure = type.measure,alpha=1)
      out_xx<-islasso(y~xx,data=data_df,lambda=fit_xx$lambda.min)
      out2<-data.frame(summary(out_xx)$coefficients,stringsAsFactors = F) %>%
        rownames_to_column(var="varx")
      out_yo[['ps_medi']] = out2
    ################################################################
      if(verb) print(paste("...finish mediation analysis."))
    ################################################################
    }
  }
  # calculate composit weight
  expr<-paste0("tw_comp = ",paste(yo_vec,collapse="+"))
  tw_composit %<>%
    mutate(tw_comp := !!rlang::parse_expr(expr)) %>%
    mutate(tw_comp_adj = rescale(tw_comp, to=c(0.01,0.99)))
  ################################################################
  if(verb) print("finish generating composit weights.")
  ################################################################
  # stack result
  out[['composit']]<-list(
    ps_tw = tw_composit
  )
  return(out)
}

ipw.rf<-function(){

}

ipw.gbm<-function(){

}

ipw.ensemble<-function(){

}

fast_rfe.coxph<-function(
  data_df, # data.frame including yc, x_tw, xo_vec
  time_col='time', # time column reuiqred for Surv() object
  status_col='status', # stauts column reuiqred for Surv() object
  yc = 'TRT', # column name of exposure at center,
  x_wt = 'wt', # column with weights
  xo_vec = c(""), # vector of other covariates to be selected
  pval_threshold = 0.01,
  verb = TRUE # verbose
){
  # require(tidyverse,magrittr,survival,survminer)

  insig_n<-length(xo_vec) # init: all other covariates
  var_sel<-xo_vec # init: start with all covariates
  pval<-0 
  while(insig_n > 0 & length(var_sel) > 0 & pval <= pval_threshold){
    # build model
    fit_frm<-formula(paste0("Surv(",time_col,",",status_col,") ~ ",
                            paste(c(var_sel, yc), collapse = "+")))
    wt<-unlist(data_df[,x_wt])
    fit_mort_msm<-coxph(fit_frm, data = data_df, weights = wt)
    fit_mort_summ<-summary(fit_mort_msm)$coefficients
    
    # update significant feature list
    var_sel<-row.names(fit_mort_summ)[fit_mort_summ[,5]<=pval_threshold&!is.na(fit_mort_summ[,5])]
    insig_n<-nrow(fit_mort_summ) - length(var_sel) 
    pval<-fit_mort_summ[yc,5]
    
    # report progress when needed
    if(verb){
      print(paste0("significant variables:",length(var_sel),";",
                   "insignificant variables:",insig_n))
    }
  }
  # build final model with selected feature set
  fit_sel<-coxph(formula(paste0("Surv(",time_col,",",status_col,") ~ ",
                                 paste(unique(c(var_sel, yc)), collapse = "+"))),
                  data = data_df, weights = wt)
  # result list
  out<-list(
    var_sel = var_sel,
    fit_sel = fit_sel
  )

  return(out)
}

bayeopt_xgb<-function(
  dtrain,
  dtest,
  folds=5,
  params_bd=list(
    max_depth = c(4L, 10L),
    min_child_weight = c(2L,10L),
    subsample = c(0.5,0.8),
    colsample_bytree=c(0.3,0.8),
    eta=c(0.05,0.1)
  ),
  N_CL=1,
  verb=T
){
  # require(xgboost,ParBayesianOptimization)
  
  #--parallelization
  if(N_CL > 1){
    cl <- makeCluster(N_CL)
    registerDoParallel(cl)
    clusterExport(cl,'df') # copying data to clusters (note:xgb.DMatrix is not compatible with parallelization)
    clusterEvalQ(cl,expr= {                          
      library(xgboost) # copying model to clusters
    })
  }
  
  #--tune hyperparameter (less rounds, early stopping)
  xgb_cv_bayes <- function(
    max_depth=10L, min_child_weight=1L, subsample=0.7,
    eta=0.05,colsample_bytree=0.8,lambda=1,alpha=0,gamma=1
  ) {  
    cv <- xgb.cv(
      params = list(
        booster = "gbtree",
        max_depth = max_depth,
        min_child_weight = min_child_weight,
        subsample = subsample, 
        eta = eta,
        colsample_bytree = colsample_bytree,
        lambda = lambda,
        alpha = alpha,
        gamma = gamma,
        objective = "binary:logistic",
        eval_metric = "auc"
      ),
        data = dtrain,
        nround = 100,
        folds = folds,
        prediction = FALSE,
        # showsd = TRUE,
        early_stopping_rounds = 5,
        maximize = TRUE,
        verbose = 0
    )
    return(list(Score = cv$evaluation_log$test_auc_mean[cv$best_iteration]))
  }
  
  OPT_Res <- bayesOpt(
    FUN = xgb_cv_bayes,
    bounds = bounds,
    initPoints = 5,
    iters.n = 50,
    iters.k = N_CL,
    parallel = TRUE,
    acq = "ucb",
    kappa = 2.576,
    eps = 0.0,
    otherHalting = list(timeLimit = 18000) #--limit maximal running time for better efficiency-- <5hr
  )
  
  Best_Par<-getBestPars(OPT_Res)
  
  #--stop cluster
  if(N_CL > 1){
    stopCluster(cl)
    registerDoSEQ()
  }

  # returen best parameters
  return(Best_Par)
}

prune_xgb<-function(
  dtrain,
  dtest,
  params=list(
    booster = "gbtree",
    max_depth = 10,
    min_child_weight = 2,
    colsample_bytree = 0.8,
    subsample = 0.7,
    eta = 0.05,
    lambda = 1,
    alpha = 0,
    gamma = 1,
    objective = "binary:logistic",
    eval_metric = "auc"
  ),
  nround = 500,
  folds = folds,
  # nfold = 5,
  prediction = FALSE,
  showsd = FALSE,
  early_stopping_rounds = 50,
  maximize = TRUE,
  verbose = 1,
  print_every_n=50

){
  #--determine number of trees, or steps (more rounds, early stopping)
  bst <- xgb.cv(
    params = list(
      booster = params$booster,
      max_depth = params$max_depth,
      min_child_weight = params$min_child_weight,
      colsample_bytree = params$colsample_bytree,
      subsample = params$subsample,
      eta = params$eta,
      lambda = params$lambda,
      alpha = params$alpha,
      gamma = params$gamma,
      objective = params$objective,
      eval_metric = params$eval_metric
),
    data = dtrain,
    nround = nround,
    folds = folds,
    # nfold = nfold,
    prediction = prediction,
    showsd = showsd,
    early_stopping_rounds = early_stopping_rounds,
    maximize = maximize,
    verbose = verbose,
    print_every_n = print_every_n
  ) 
  # get optimal steps
  steps<-which(bst$evaluation_log$test_auc_mean==max(bst$evaluation_log$test_auc_mean))
  
  #--prune the optimal model
  xgb_tune<-xgb.train(
    data=dtrain,
    max_depth = params$max_depth,
    min_child_weight = params$min_child_weight,
    colsample_bytree = params$colsample_bytree,
    subsample = params$subsample,
    eta = params$eta,
    maximize = maximize,
    nrounds = steps,
    eval_metric = params$eval_metric,
    objective = params$objective,
    verbose = verbose
  )
  
  #--collect training results
  valid_tr<-data.frame(
    id = row.names(trainX),
    actual = getinfo(dtrain,"label"),
    pred = predict(xgb_tune,dtrain),
    stringsAsFactors = F
  )

  #--collect testing results
  valid_ts<-data.frame(
    id = row.names(testX),
    actual = getinfo(dtest,"label"),
    pred = predict(xgb_tune,dtest),
    stringsAsFactors = F
  )
  
  #--feature importance
  feat_imp<-xgb.importance(model=xgb_tune)
  
  #--save model and other results
  result<-list(
    model = xgb_tune,
    pred_tr = valid_tr,
    pred_ts = valid_ts,
    feat_imp = feat_imp
  )
  return(result)
}

explain_model<-function(
  X,y,
  xgb_rslt, # model, pred_df, feat_imp
  top_k = 20,
  var_lst = c(),
  boots = 10,
  nns = 30,
  shap_cond = NULL, # time index
  verb = TRUE
){
  
  #identify top k features
  var_imp<-xgb_rslt$feat_imp %>% 
    slice(seq_len(top_k)) %>%
    select(Feature, Gain)
  var_nm<-var_imp$Feature

  #add conditional index feature
  if(!is.null(shap_cond)){
    var_nm<-unique(c(shap_cond,var_nm))
  }
  
  #add custom feature set
  if(!is.null(var_lst)){
    var_nm<-unique(c(var_lst,var_nm))
  }

  #bootstrap CI for SHAP values 
  pred_brkdn_b<-c()
  x_val_b<-c()
  
  for(b in 1:boots){
    # stratified sampling
    n_idx<-which(y==0)
    p_idx<-which(y==1)
    nn<-length(n_idx)
    idxset<-c(p_idx,n_idx[sample(1:nn,nns,replace=F)])
    
    # get shap contribution
    contr <- predict(
      xgb_rslt$model,
      newdata=X[idxset,],
      predcontrib = TRUE
    )
    shap_sel<-contr[,which(colnames(contr) %in% var_nm)]
    
    # stack bootstrapping results
    pred_brkdn_b %<>%
      bind_rows(cbind(as.data.frame(shap_sel),
                      boot=b,idx=idxset))
    
    x_val_b %<>%
      bind_rows(cbind(as.data.frame(as.matrix(X[idxset,which(colnames(X) %in% var_nm)])),
                      boot=b,idx=idxset))
    
    # report progress
    if(verb){
      cat("bootstrapped shap values generated:",b,"./n")
    }
  }
  
  pred_brkdn<-c()
  var_lst<-colnames(shap_sel)
  for(v in seq_along(var_lst)){
    pred_brkdn_v<-pred_brkdn_b %>%
      select(all_of(c(var_lst[v],"boot","idx"))) 
  
    if(!is.null(shap_cond)){
      pred_brkdn_v %<>%
        mutate(
          val=round(x_val_b[,var_lst[v]],2),
          cond=x_val_b[,shap_cond]
        ) %>%
        group_by(boot,val,cond)
    }else{
      pred_brkdn_v %<>%
        mutate(
          val=round(x_val_b[,var_lst[v]],2)
        ) %>%
        group_by(boot,val)
    }
    
    pred_brkdn_v %<>%
      summarise(
        effect=mean(get(var_lst[v])),
        .groups = "drop") %>%
      mutate(var=var_lst[v])
      
    pred_brkdn %<>% bind_rows(pred_brkdn_v)

    # report progress
    if(verb){
      cat("shap value generated for:",var_lst[v],"./n")
    }
  }

  # result set
  return(pred_brkdn)
}

iptw_calc<-function(
  wt_long, # unit of obs-tgt-time per row
  id_col = 'id', # name of id column
  time_col = 'time', # name of time index column
  tgt_col = 'tgt', # name of propensity score target
  wt_col = 'wt', # name of the value column
  use_stablizer = TRUE,
  truncate = FALSE,
  truncate_lower = 0.01,
  truncate_upper = 0.99
){
  ext_nm<-c(id_col,time_col,tgt_col,wt_col)
  int_nm<-c('id','time','tgt','wt')

  # calculate cumulative product
  wt_df<-wt_long %>%
    rename_at(vars(ext_nm), ~ int_nm)
    group_by(id,tgt) %>%
    arrange(time) %>% 
    mutate(cum_wt = cumprod(wt)) %>%
    ungroup

  #  calculate numerator
  wt_num<-wt_df %>% 
      group_by(id,tgt,time) %>%
      summarise(num = mean(wt,.groups = 'drop')) 
  if(use_stablizer){
    wt_num<-wt_num %>%
      group_by(id,tgt) %>%
      mutate(num_lag = lag(num,n=1L,order_by=time)) %>% 
      mutate(num = num_lag) %>% select(-num_lag)
  }else{
    wt_num<-wt_num %>% 
      mutate(num = 1)
  }

  # calculate iptw ratio
  wt_df %<>%
    inner_join(wt_num,by=c("id","tgt","time")) %>%
    mutate(iptw = num/wt)

  # truncation
  if(truncate){
    wt_df %<>%
      arrange(iptw) %>%
      mutate(
        lb=quantile(iptw,probs=truncate_lower),
        ub=quantile(iptw,probs=truncate_upper)
      ) %>% 
      mutate(iptw = pmin(pmax(lb,iptw),ub))
  }

  # convert column name back
  wt_df %<>% 
    rename_at(vars(int_nm), ~ ext_nm)
  return(wt_df)
}

# bayesopt_rf<-function(

# ){
# # https://www.randomforestsrc.org/articles/survival.html

# }

# bayesopt_rf.surv<-function(){}

# generate_knockoff_copy<-function(dataX,y,burnin=1000,fdr_threshold=0.2){
#   N <- nrow(dataX)
#   factorIndi<- sapply(seq(1, ncol(dataX)), function(x){nlevels(dataX[,x])})
#   factorIndi[factorIndi==0]<-1
#   contNum<-sum(factorIndi==1)
#   cateNum<-sum(factorIndi>1)
  
  
#   cateVar<-colnames(dataX)[factorIndi>1]
#   contVar<-colnames(dataX)[factorIndi==1]

#   contX<-as.matrix(scale(dataX[,contVar]))
#   cateX<-dataX[,cateVar]
#   cateX<-matrix(unlist(lapply(1:cateNum, function(x){as.numeric(cateX[,x])})),ncol=cateNum)
  
#   cateDummyNum<-sum(factorIndi[factorIndi>1])
#   lengthZ<-contNum+cateDummyNum-cateNum
  
#   cumCateDummyNum<-c(0,cumsum(factorIndi[factorIndi>1]))
#   cumCateDummyNumInter<-c(0,cumsum(factorIndi[factorIndi>1]-1))
#   mu0<-rep(0,lengthZ)
#   tau0<-diag(rep(0.1/lengthZ,lengthZ))
#   Z<-mu0
#   Sigma<-rwish(lengthZ,tau0)
#   initsfun <-function() list(Z=Z, tau=1,Sigma=Sigma)
  
#   data <- list(N=N, contX=contX, cateX= cateX, contNum=contNum, cateNum=cateNum, cateDummyNum=cateDummyNum,
#                cumCateDummyNum=cumCateDummyNum, cumCateDummyNumInter=cumCateDummyNumInter, mu0=mu0, tau0=tau0,lengthZ=lengthZ)
#   fit_rstan <- stan(
#     file = "internalStan.stan",
#     data = data,
#     init=initsfun,
#     chains=1,
#     warmup=burnin,
#     iter=N+burnin
#   )
  
#   ZMCMC<-extract(fit_rstan,pars=paste0("Z[",1:lengthZ,"]"))
#   ZMCMC<-matrix(unlist(ZMCMC),ncol=lengthZ,byrow=FALSE) # ZMCMC is a N*d matrix. columns correspond to covariates.
#   tauMCMC<-unlist(extract(fit_rstan,pars="tau"))
#   samps<-cbind(ZMCMC,tauMCMC)
  
  
#   contMCMC<- ZMCMC[,1:contNum]
#   cateMCMC<- ZMCMC[,(contNum+1):lengthZ]
#   contSample<-matrix(unlist(lapply(1:N,function(x) mvrnorm(1,contMCMC[x,],diag(rep(tauMCMC[x],contNum))))),
#                      ncol=contNum,byrow=TRUE)
#   cateSample<-matrix(NA,nrow=N,ncol=cateDummyNum)
  
#   for(n in 1:N){
#     for(i in 1:cateNum){
#       temp<-c(1,cateMCMC[n,(1+cumCateDummyNumInter[i]):cumCateDummyNumInter[i+1]])
#       temp<-exp(temp)/sum(exp(temp))
#       cateSample[n,(1+cumCateDummyNum[i]):cumCateDummyNum[i+1]]<- rmultinom(1, 1, temp)
#     }
#   }
  
#   Group<-c(contVar,unlist(sapply(1:cateNum,function(x) rep(cateVar[x],factorIndi[factorIndi>1][x]))))
#   colnames(cateX)<-cateVar
#   Xoriginal<-cbind(contX,dummy_cols(cateX,select_columns = cateVar,remove_selected_columns = TRUE))
#   Xknockoff<-cbind(contSample,cateSample)
  
#   W.output <- stat.grpreg_penaltydiff(Xoriginal,Xknockoff, y, family='gaussian',Group)
#   t <- knockoff.group.threshold(W.output=W.output, fdr=fdr_threshold, offset=0)
#   SelectedVarNames<-c(contVar,cateVar)[t$subset.selected]
  
#   selectedResult<-list(SelectedVarNames=SelectedVarNames,samps=samps)
  
#   selectedResult
# }

# generate_boruta_copy<-function(){
#  
# }

# #TODO
# adjMMD<-function(){
#   
# }
