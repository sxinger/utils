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
  # conversion to matrix
  x<-data.matrix(data_df[,xo_vec])
  
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
    fit_tw<-cv.glmnet(x=x,y=y,family=family,type.measure = type.measure,alpha=1)
    tw<-predict(fit_tw, newx = x, s = "lambda.min", type="response")
    id<-unlist(data_df[,id_col])
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
  ################################################################
  # calculate composit weight
  expr<-paste0("tw_comp=",paste(yo_vec,collapse="*"))
  tw_composit %<>%
    mutate(tw_comp := !!rlang::parse_expr(expr)) %>%
    mutate(idx=row_number()) %>%
    arrange(tw_comp) %>% 
    mutate(tw_comp_adj=rank(tw_comp)/n()) %>%
    mutate(tw_comp_adj=case_when(y==1 ~ tw_comp_adj,
                                 y==0&tw_comp_adj==1 ~ 1-(tw_comp_adj-0.0001),
                                 TRUE ~ 1-tw_comp_adj)) %>%
    arrange(idx) %>% select(id,tw_comp,tw_comp_adj)
  ################################################################
  if(verb) print("finish generating composit weights.")
  ################################################################
  # stack result
  out[['composit']]<-list(
    ps_tw = tw_composit
  )
  return(out)
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
  insig_n<-length(xo_vec) # init: all other covariates
  var_sel<-var_ps # init: start with all covariates
  pval<-0 
  while(insig_n > 0 & length(var_sel) > 0 & pval <= pval_threshold){
    # build model
    fit_frm<-formula(paste0("Surv(",time_col,",",status_col,") ~ ",
                            paste(c(var_sel, yc), collapse = "+")))
    wt<-unlist(data_df[,x_wt])
    fit_mort_msm<-coxph(fit_frm, data = data_df, weights = wt)
    fit_mort_summ<-summary(fit_mort_msm)$coefficients
    
    # update significant feature list
    var_sel<-row.names(fit_mort_summ)[fit_mort_summ[,6]<=pval_threshold&!is.na(fit_mort_summ[,6])]
    insig_n<-nrow(fit_mort_summ) - length(var_sel) 
    pval<-fit_mort_summ[yc,6]
    
    # report progress when needed
    if(verb){
      print(paste0("significant variables:",length(var_sel),";",
                   "insignificant variables:",insig_n))
    }
  }
  return(var_sel)
}

# bayeopt_xgb<-function(
#   df_long,
#   params_bd=list(
#     max_depth = c(4L, 10L),
#     min_child_weight = c(2L,10L),
#     subsample = c(0.5,0.8),
#     colsample_bytree=c(0.3,0.8),
#     eta=c(0.05,0.1)
#   ),
#   N_CL=1,
#   verb=T

# ){
#   {
#   bm<-c()
#   bm_nm<-c()
#   start_tsk_i<-Sys.time()
  
#   #--parallelization
#   cl <- makeCluster(N_CL)
#   registerDoParallel(cl)
#   clusterExport(cl,'df') # copying data to clusters (note:xgb.DMatrix is not compatible with parallelization)
#   clusterEvalQ(cl,expr= {                          # copying model to clusters
#     library(xgboost)
#   })

#   y_df<-df[,which(colnames(df) %in% c("id","y"))]
#   X_df<-df[,which(!colnames(df) %in% c("id","y","fold"))]
  
#   #--covert to xgb data frame
#   dtrain<-xgb.DMatrix(data=X_df[,],
#                       label=y_df$y)
#   dtest<-xgb.DMatrix(data=X_df,
#                      label=y_df$y)
  
#   #-----------------------------------------------benchmark-----------------------------------------#
#   lapse_i<-Sys.time()-start_tsk_i
#   bm<-c(bm,paste0(round(lapse_i,1),units(lapse_i)))
#   bm_nm<-c(bm_nm,"transform data")
#   if(verb){
#     cat(paste0(c(pred_in_d,pred_task,fs_type),collapse = ","),
#         "...finish formatting training and testing sets.\n")
#   }
#   #-----------------------------------------------benchmark------------------------------------------#
  
  
#   #--tune hyperparameter (less rounds, early stopping)
#   xgb_cv_bayes <- function(max_depth=10L, min_child_weight=1L, subsample=0.7,
#                            eta=0.05,colsample_bytree=0.8,lambda=1,alpha=0,gamma=1) {
    
#     dtrain<-xgb.DMatrix(data=X_tr_sp,label=y_tr_sp$y)
#     cv <- xgb.cv(params = list(booster = "gbtree",
#                                max_depth = max_depth,
#                                min_child_weight = min_child_weight,
#                                subsample = subsample, 
#                                eta = eta,
#                                colsample_bytree = colsample_bytree,
#                                lambda = lambda,
#                                alpha = alpha,
#                                gamma = gamma,
#                                objective = "binary:logistic",
#                                eval_metric = "auc"),
#                  data = dtrain,
#                  nround = 100,
#                  folds = folds,
#                  prediction = FALSE,
#                  # showsd = TRUE,
#                  early_stopping_rounds = 5,
#                  maximize = TRUE,
#                  verbose = 0)
    
#     return(list(Score = cv$evaluation_log$test_auc_mean[cv$best_iteration]))
#   }
  
#   OPT_Res <- bayesOpt(
#     FUN = xgb_cv_bayes,
#     bounds = bounds,
#     initPoints = 5,
#     iters.n = 50,
#     iters.k = N_CL,
#     parallel = TRUE,
#     acq = "ucb",
#     kappa = 2.576,
#     eps = 0.0,
#     otherHalting = list(timeLimit = 18000) #--limit maximal running time for better efficiency-- <5hr
#   )
  
#   Best_Par<-getBestPars(OPT_Res)
  
#   #--stop cluster
#   stopCluster(cl)
#   registerDoSEQ()
  
#   #--determine number of trees, or steps (more rounds, early stopping)
#   dtrain<-xgb.DMatrix(data=X_tr_sp,label=y_tr_sp$y)
#   bst <- xgb.cv(params = list(booster = "gbtree",
#                               max_depth = Best_Par$max_depth,
#                               min_child_weight = Best_Par$min_child_weight,
#                               colsample_bytree = Best_Par$colsample_bytree,
#                               subsample=0.7,
#                               eta=0.05,
#                               lambda=1,
#                               alpha=0,
#                               gamma=1,
#                               objective = "binary:logistic",
#                               eval_metric = "auc"),
#                 data = dtrain,
#                 nround = 500,
#                 folds = folds,
#                 # nfold=5,
#                 prediction = FALSE,
#                 # showsd = TRUE,
#                 early_stopping_rounds = 50,
#                 maximize = TRUE,
#                 verbose = 1,
#                 print_every_n=50) 
  
#   steps<-which(bst$evaluation_log$test_auc_mean==max(bst$evaluation_log$test_auc_mean))
  
#   lapse_i<-Sys.time()-start_tsk_i
#   bm<-c(bm,paste0(round(lapse_i,1),units(lapse_i)))
#   bm_nm<-c(bm_nm,"tune model")
  
#   cat(paste0(c(pred_in_d,pred_task,fs_type),collapse = ","),
#       "...finish model tunning.\n")
  
#   #-----------validate model------------
#   start_tsk_i<-Sys.time() 
  
#   #--validation
#   xgb_tune<-xgb.train(data=dtrain,
#                       max_depth = Best_Par$max_depth,
#                       min_child_weight = Best_Par$min_child_weight,
#                       colsample_bytree = Best_Par$colsample_bytree,
#                       subsample=0.7,
#                       eta=0.05,
#                       maximize = TRUE,
#                       nrounds=steps,
#                       eval_metric="auc",
#                       objective="binary:logistic",
#                       verbose = 0)
  
#   valid<-data.frame(y_ts_sp,
#                     pred = predict(xgb_tune,dtest),
#                     stringsAsFactors = F)
  
#   #--feature importance
#   feat_imp<-xgb.importance(model=xgb_tune)
  
#   lapse_i<-Sys.time()-start_tsk_i
#   bm<-c(bm,paste0(round(lapse_i,1),units(lapse_i)))
#   bm_nm<-c(bm_nm,"validate model")
  
#   cat(paste0(c(pred_in_d,pred_task,fs_type),collapse = ","),
#       "...finish model validating.\n")
  
#   #-----------save model and other results--------
#   result<-list(hyper_param=c(Best_Par,steps),
#                model=xgb_tune,
#                pred_df=valid,
#                feat_imp=feat_imp)
  
#   #-------------------------------------------------------------------------------------------------------------
#   lapse_tsk<-Sys.time()-start_tsk
#   bm<-c(bm,paste0(round(lapse_tsk,1),units(lapse_tsk)))
#   bm_nm<-c(bm_nm,"complete task")
  
#   cat("\nFinish building reference models for task:",pred_task,"in",pred_in_d,"with",fs_type,",in",lapse_tsk,units(lapse_tsk),
#       ".\n--------------------------\n")
  
#   #benchmark
#   bm<-data.frame(bm_nm=bm_nm,bm_time=bm,
#                  stringsAsFactors = F)
# }

# }

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


# explain_model<-function(
# 
# ){
#   #load trained model
#   gbm_model<-gbm_ctnr$model
  
#   #identify top k features
#   var_imp<-xgb.importance(model=gbm_model) %>% dplyr::slice(k_seq) %>%
#     select(Feature, Gain)
  
#   #bootstrap CI for SHAP values
#   boots<-100
#   nns<-10000
  
#   pred_brkdn_b<-c()
#   x_val_b<-c()
  
#   for(b in 1:boots){
#     start_b<-Sys.time()
    
#     #stratified sampling
#     n_idx<-which(y_ts_sp$y==0)
#     p_idx<-which(y_ts_sp$y==1)
#     nn<-length(n_idx)
#     idxset<-c(p_idx,n_idx[sample(1:nn,nns,replace=F)])
    
#     contr <- predict(gbm_model,
#                       newdata=X_ts_sp[idxset,],
#                       predcontrib = TRUE)
    
#     shap_sel<-contr[,which(colnames(contr) %in% var_nm)]
    
#     #careful!! rows get re-ordered by xgb.plot.shap
#     # shap<-xgb.plot.shap(data=X_ts_sp[idxset,],
#     #                     shap_contrib=contr,
#     #                     model = gbm_model,
#     #                     top_n = 10,
#     #                     plot=F)
    
#     pred_brkdn_b %<>%
#       bind_rows(cbind(as.data.frame(shap_sel),
#                       boot=b,idx=idxset))
    
#     x_val_b %<>%
#       bind_rows(cbind(as.data.frame(as.matrix(X_ts_sp[idxset,which(colnames(X_ts_sp) %in% var_nm)])),
#                       boot=b,idx=idxset))
    
#     lapse<-Sys.time()-start_b
#     cat(paste0(c(params$site,pred_in_d,pred_task,fs_type),collapse = ","),
#         "...finish bootstrapped sample",b,"in",lapse,units(lapse),".\n")
#   }
  
#   pred_brkdn<-c()
#   var_lst<-colnames(shap_sel)
#   for(v in seq_along(var_lst)){
#     pred_brkdn %<>%
#       bind_rows(pred_brkdn_b %>%
#                   dplyr::select(var_lst[v],"boot","idx") %>%
#                   dplyr::mutate(val=round(x_val_b[,var_lst[v]],2)) %>%
#                   group_by(boot,val) %>%
#                   dplyr::summarise(effect=mean(get(var_lst[v]))) %>%
#                   ungroup %>%
#                   dplyr::mutate(var=var_lst[v]))
#   }
# }

# #TODO
# adjMMD<-function(){
#   
# }
