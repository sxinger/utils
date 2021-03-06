#########################################
## utility functions for data analysis ##
#########################################

# multiclass y is not supported yet!
# data_type should be a vector of "cat" or "num"
#require (purrr,broom,tibble)
univar_analysis_mixed<-function(id,grp,X,data_type,pretty=F){
  if(ncol(X)!=length(data_type)){
    stop("data types of X need to be specified")
  }
  
  #TODO: when there is only 1 category
  
  # anova
  df_num<-data.frame(cbind(id,grp,X[,(data_type=="num"),drop=F]),stringsAsFactors=F) %>%
    gather(var,val,-grp,-id) %>%
    mutate(grp=as.factor(grp)) %>%
    mutate(val=as.numeric(val))
  
  out_num<-df_num %>%
    group_by(var,grp) %>%
    dplyr::summarise(n=length(unique(id)),
                     val_miss=sum(is.na(val)),
                     val_mean=mean(val,na.rm=T),
                     val_sd=sd(val,na.rm=T),
                     val_med=median(val,na.rm=T),
                     val_q1=quantile(val,0.25,na.rm=T),
                     val_q3=quantile(val,0.75,na.rm=T),
                     val_min=min(val,na.rm=T),
                     val_max=max(val,na.rm=T)) %>% 
    ungroup %>%
    left_join(df_num %>%
                nest(-var) %>%
                mutate(fit=map(data, ~ aov(val~grp,data=.x)),
                       tidied=map(fit,tidy)) %>%
                unnest(tidied) %>% 
                filter(!is.na(p.value)) %>%
                select(var,p.value),
              by="var") %>%
    mutate(label=paste0(n,"; ",
                        # round(val_miss/n,2),"; ", #missing rate
                        round(val_mean,1),"(",round(val_sd,2),"); ",
                        val_med,"(",val_q1,",",val_q3,")"))
  
  
  # chi-sq
  df_cat<-data.frame(cbind(id,grp,X[,(data_type=="cat")]),stringsAsFactors=F) %>%
    gather(var,val,-grp,-id) %>%
    mutate(grp=as.factor(grp),val=as.factor(val))
  
  out_cat<-df_cat %>%
    group_by(grp) %>%
    dplyr::mutate(tot=length(unique(id))) %>%
    ungroup %>%
    group_by(var) %>%
    dplyr::mutate(val_miss=sum(is.na(val))) %>%
    ungroup %>% filter(!is.na(val)) %>%
    group_by(var,grp,tot,val_miss,val) %>%
    dplyr::summarise(n=length(unique(id))) %>%
    ungroup %>%
    mutate(prop=round(n/tot,4)) %>%
    left_join(df_cat %>%
                group_by(var) %>%
                dplyr::summarise(p.value=chisq.test(val,grp,simulate.p.value=T)$p.value) %>%
                ungroup,
              by="var") %>%
    mutate(label=paste0(n,"; ",
                        # round(val_miss/n,2),"; ", #missing rate
                        "(",prop*100,"%)"))
  
  #output
  if(pretty){
    out<-out_num %>% 
      select(n,grp) %>% unique %>%
      gather(var,val,-grp) %>% 
      mutate(val=as.character(val)) %>% 
      spread(grp,val) %>%
      bind_rows(out_num %>%
                  mutate(label2=paste0(round(val_mean,1)," (",round(val_sd,1),")"," [",round(val_miss/n,2),"]")) %>%
                  dplyr::select(var,grp,p.value,label2) %>% spread(grp,label2)) %>%
      bind_rows(out_cat %>%
                  unite("var",c("var","val"),sep="=") %>%
                  mutate(label2=paste0(n," (",round(prop*100,1),"%)"," [",round(val_miss/n,2),"]")) %>%
                  dplyr::select(var,grp,p.value,label2) %>% spread(grp,label2)) %>%
      mutate(p.value=round(p.value,4)) %>%
      separate("var",c("var","cat"),sep="=",extra="merge",fill="right") %>%
      mutate(cat=case_when(var=="n" ~ "",
                           is.na(cat) ~ "mean(sd) [miss]",
                           TRUE ~ paste0(cat,",n(%) [miss]")))
    
  }else{
    out<-list(out_num=out_num,
              out_cat=out_cat)
  }
  
  return(out)
}

get_perf_summ<-function(pred,real,keep_all_cutoffs=F){
  # various performace table
  pred_obj<-ROCR::prediction(pred,real)
  
  prc<-performance(pred_obj,"prec","rec")
  roc<-performance(pred_obj,"sens","spec")
  nppv<-performance(pred_obj,"ppv","npv")
  pcfall<-performance(pred_obj,"pcfall")
  acc<-performance(pred_obj,"acc")
  fscore<-performance(pred_obj,"f")
  mcc<-performance(pred_obj,"phi")
  
  perf_at<-data.frame(cutoff=prc@alpha.values[[1]],
                      prec=prc@y.values[[1]],
                      rec_sens=prc@x.values[[1]],
                      stringsAsFactors = F) %>% 
    arrange(cutoff) %>%
    left_join(data.frame(cutoff=nppv@alpha.values[[1]],
                         ppv=nppv@y.values[[1]],
                         npv=nppv@x.values[[1]],
                         stringsAsFactors = F),
              by="cutoff") %>%
    dplyr::mutate(prec_rec_dist=abs(prec-rec_sens)) %>%
    left_join(data.frame(cutoff=fscore@x.values[[1]],
                         fscore=fscore@y.values[[1]],
                         stringsAsFactors = F),
              by="cutoff") %>%
    left_join(data.frame(cutoff=roc@alpha.values[[1]],
                         spec=roc@x.values[[1]],
                         stringsAsFactors = F),
              by="cutoff") %>%
    dplyr::mutate(Euclid_meas=sqrt((1-rec_sens)^2+(0-(1-spec))^2),
                  Youden_meas=rec_sens+spec-1) %>%
    left_join(data.frame(cutoff=pcfall@x.values[[1]],
                         pcfall=pcfall@y.values[[1]],
                         stringsAsFactors = F),
              by="cutoff") %>%
    left_join(data.frame(cutoff=acc@x.values[[1]],
                         acc=acc@y.values[[1]],
                         stringsAsFactors = F),
              by="cutoff") %>%
    left_join(data.frame(cutoff=mcc@x.values[[1]],
                         mcc=mcc@y.values[[1]],
                         stringsAsFactors = F),
              by="cutoff") %>%
    filter(prec > 0 & rec_sens > 0 & spec > 0) %>%
    group_by(cutoff) %>%
    dplyr::mutate(size=n()) %>%
    ungroup
  
  # performance summary
  lab1<-pred[real==1]
  lab0<-pred[real==0]
  pr<-pr.curve(scores.class0 = lab1,
               scores.class1 = lab0,curve=F)
  roc_ci<-pROC::ci.auc(real,pred)
  
  perf_summ<-data.frame(overall_meas=c("roauc_low",
                                       "roauc",
                                       "roauc_up",
                                       "opt_thresh",
                                       "opt_sens",
                                       "opt_spec",
                                       "opt_ppv",
                                       "opt_npv",
                                       "prauc1",
                                       "prauc2",
                                       "opt_prec",
                                       "opt_rec",
                                       "opt_fscore"),
                        meas_val=c(roc_ci[[1]],
                                   roc_ci[[2]],
                                   roc_ci[[3]],
                                   perf_at$cutoff[which.min(perf_at$Euclid_meas)],
                                   perf_at$rec_sens[which.min(perf_at$Euclid_meas)],
                                   perf_at$spec[which.min(perf_at$Euclid_meas)],
                                   perf_at$ppv[which.min(perf_at$Euclid_meas)],
                                   perf_at$npv[which.min(perf_at$Euclid_meas)],
                                   pr$auc.integral,
                                   pr$auc.davis.goadrich,
                                   perf_at$prec[which.min(perf_at$prec_rec_dist)],
                                   perf_at$rec_sens[which.min(perf_at$prec_rec_dist)],
                                   perf_at$fscore[which.min(perf_at$prec_rec_dist)]),
                        stringsAsFactors = F) %>%
    bind_rows(perf_at %>% 
                dplyr::summarize(prec_m=mean(prec,na.rm=T),
                                 sens_m=mean(rec_sens,na.rm=T),
                                 spec_m=mean(spec,na.rm=T),
                                 ppv_m=mean(ppv,na.rm=T),
                                 npv_m=mean(npv,na.rm=T),
                                 acc_m=mean(acc,na.rm=T),
                                 fscore_m=mean(fscore,na.rm=T),
                                 mcc_m=mean(mcc,na.rm=T)) %>%
                gather(overall_meas,meas_val))
  
  out<-list(perf_summ=perf_summ)
  if(keep_all_cutoffs){
    out$perf_at<-perf_at
  }
  
  return(out)
}

get_calibr<-function(pred,real,n_bin=20){
  calib<-data.frame(pred=pred,
                    y=real) %>%
    arrange(pred) %>%
    dplyr::mutate(pred_bin = cut(pred,
                                 breaks=unique(quantile(pred,0:(n_bin)/(n_bin))),
                                 include.lowest=T,
                                 labels=F)) %>%
    ungroup %>% group_by(pred_bin) %>%
    dplyr::summarize(expos=n(),
                     bin_lower=min(pred),
                     bin_upper=max(pred),
                     bin_mid=median(pred),
                     y_agg = sum(y),
                     pred_p = mean(pred)) %>%
    dplyr::mutate(y_p=y_agg/expos) %>%
    dplyr::mutate(binCI_lower = pmax(0,pred_p-1.96*sqrt(y_p*(1-y_p)/expos)),
                  binCI_upper = pred_p+1.96*sqrt(y_p*(1-y_p)/expos))
  
  return(calib)
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

