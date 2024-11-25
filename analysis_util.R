#########################################
## utility functions for data analysis ##
#########################################

# find turning point 
get_turnpoints<-function(
  x,y
){
  # require(tidyverse)
  delta_y = diff(y)
  turns_ind = which(delta_y[-1] * delta_y[-length(delta_y)] < 0) + 1
  turn_df = data.frame(x=x[turns_ind],y=y[turns_ind])
  return(turn_df)
}


# p-value calculators based on sample summaries
pval_on_summ_2sampleprop<-function(
  v1, #=c(n1,p1)
  v2, #=c(n2,p2)
  pooled=TRUE,
  alternative="both"
){
  # require(tidyverse,magrittr)
  n1<-v1[1]
  p1<-v1[2]
  n2<-v2[1]
  p2<-v2[2]
  if(pooled){
    p<-(p1*n1+p2*n2)/(n1+n2)
    se<-sqrt(p*(1-p)*(1/n1+1/n2))
  }else{
    se<-sqrt(p1*(1-p1)/n1+p2*(1-p2)/n2)
  }
  z<-(p1-p2)/se
  
  if(alternative=="lower"){
    pval<-pnorm(z,lower.tail = TRUE)
  }else if(alternative=="upper"){
    pval<-pnorm(z,lower.tail = FALSE)
  }else{
    pval<-2*pnorm(abs(z),lower.tail = FALSE)
  }
  
  return(data.frame(z=z,pval=pval))
}

pval_on_summ_2samplemean<-function(
  v1, #=c(n1,m1,s1)
  v2, #=c(n2,m2,s2)
  pooled=FALSE,
  alternative="both"
){
  n1<-v1[1]
  m1<-v1[2]
  s1<-v1[3]
  n2<-v2[1]
  m2<-v2[2]
  s2<-v2[3]
  if(pooled){
    se<-sqrt(((n1-1)*s1^2+(n2-1)*s2^2)/(n1+n2-2))
  }else{
    se<-sqrt(s1^2/n1+s2^2/n2)
  }
  
  t<-(m1-m2)/se
  df<-n1+n2-2
  
  if(alternative=="lower"){
    pval<-pt(t,df,lower.tail = TRUE)
  }else if(alternative=="upper"){
    pval<-pt(t,df,lower.tail = FALSE)
  }else{
    pval<-2*pt(abs(t),df,lower.tail = FALSE)
  }
  
  return(data.frame(t=t,df=df,pval=pval))
}

# multiclass y is not supported yet!  
# data_type should be a vector of "cat" or "num"
# require (purrr,broom,tibble,kable,kabelExtra)
univar_analysis_mixed<-function(
  df,
  id_col="PATID",
  grp=1,
  var_lst,
  facvar_lst,
  pretty=F,
  var_lbl_df=data.frame(var=as.character(),
                        var_lbl=as.character()) # optional, pretty=T
){
  if(!all(facvar_lst %in% var_lst)){
    stop("facvar_lst must be subset of var_lst!")
  }
  X<-df[,var_lst]
  id<-unlist(df[,id_col])
  data_type<-rep("num",length(var_lst))
  data_type[which(var_lst %in% facvar_lst,arr.ind = T)]<-"cat"

  # line 1 summary
  out<-data.frame(cbind(id,grp,X[,1,drop=F]),stringsAsFactors=F) %>% 
    gather(var,val,-grp,-id) %>%
      mutate(grp=as.factor(grp)) %>%
      group_by(var,grp) %>%
      dplyr::summarise(n=length(unique(id)),.groups = "drop") %>% 
      dplyr::select(n,grp) %>% unique %>%
      gather(var,val,-grp) %>% 
      mutate(val=as.character(val)) %>% 
      spread(grp,val)

  # anova - numerical variables
  if(length(var_lst)-length(facvar_lst)>0){
    df_num<-data.frame(cbind(id,grp,X[,(data_type=="num"),drop=F]),stringsAsFactors=F) %>%
      gather(var,val,-grp,-id) %>%
      mutate(grp=as.factor(grp)) %>%
      mutate(val=as.numeric(val))
    
    out_num<-df_num %>%
      group_by(var,grp) %>%
      dplyr::summarise(
        n=length(unique(id)),
        val_miss=sum(is.na(val)),
        val_mean=mean(val,na.rm=T),
        val_sd=sd(val,na.rm=T),
        val_med=median(val,na.rm=T),
        val_q1=quantile(val,0.25,na.rm=T),
        val_q3=quantile(val,0.75,na.rm=T),
        val_min=min(val,na.rm=T),
        val_max=max(val,na.rm=T),
        .groups = "drop"
      ) 

      if(length(unique(grp))>1){
        out_num_i<-df_num %>%
          # filter complete missing
          anti_join(
            out_num %>% filter(val_miss==1),
            by = "var"
          ) %>%
          nest(data=!var) %>%
          mutate(
            fit=map(data, ~ aov(val~grp,data=.x)),
            tidied=map(fit,tidy)
          ) %>%
          unnest(tidied) %>% 
          filter(!is.na(p.value)) %>%
          dplyr::select(var,p.value)

        out_num %<>% 
          left_join(out_num_i,by="var") 
      } else {
        out_num %<>% mutate(p.value = 1)
      }
      out_num %<>% 
        mutate(
          label=paste0(
            n,"; ",
            # round(val_miss/n,2),"; ", #missing rate
            round(val_mean,1),"(",round(val_sd,2),"); ",
            round(val_med,1),"(",round(val_q1,1),",",round(val_q3,1),")"
          )
        )
      out %<>%
        bind_rows(
          out_num %>%
            mutate(label2=paste0(
              round(val_mean,1)," (",round(val_sd,1),");",
              val_med,"(",val_q1,",",val_q3,")",
              " [",round(val_miss/n,2),"]")
            ) %>%
            dplyr::select(var,grp,p.value,label2) %>% 
            spread(grp,label2)
        )
  }else{
    out_num<-c()
  }
                    
  # chi-sq - categorical variables
  if(length(facvar_lst)>0){
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
      dplyr::summarise(n=length(unique(id)),.groups="drop") %>%
      mutate(prop=round(n/tot,4)) 
      
      if(length(unique(grp))>1){
        out_cat_i<-df_cat %>%
          group_by(var) %>%
          dplyr::summarise(
            p.value=chisq.test(val,grp,simulate.p.value=T)$p.value,
            .groups = "drop"
        )
        out_cat %<>%
          left_join(out_cat_i,by="var") 
      }else{
        out_cat %<>% mutate(p.value = 1)
      }

      out_cat %<>%
        mutate(label=paste0(n,"; ","(",prop*100,"%)"))
      
      out %<>%
        bind_rows(
          out_cat %>%
            unite("var",c("var","val"),sep="=") %>%
            mutate(label2=paste0(n," (",round(prop*100,1),"%)"," [",round(val_miss/tot,2),"]")) %>%
            dplyr::select(var,grp,p.value,label2) %>% 
            spread(grp,label2) 
        ) 
      
  }else{
    out_cat<-c()
  }
  
  # organize into tabulated format
  out %<>% 
    mutate(p.value=round(p.value,4)) %>%
    separate("var",c("var","cat"),sep="=",extra="merge",fill="right") %>%
    mutate(cat=case_when(
      var=="n" ~ "",
      is.na(cat) ~ "mean(sd); med(iqr) [miss]",
      TRUE ~ paste0(cat,",n(%) [miss]"))
    ) 
    
  #output
  if(pretty){
      # collect variable label mapping if provided
      if(nrow(var_lbl_df)>0){
        out %<>%
          left_join(var_lbl_df,by="var") %>%
          mutate(var_lbl=coalesce(var_lbl,var))
      }else{
        out %<>%
          mutate(var_lbl=var)
      } 
      # convert to html table output using kable
      colnames(out)<-c("var","cat",paste0("exposure=",sort(unique(grp))),"p.value","var_lbl")
      out %<>% 
        mutate(var_fac=factor(var,ordered = TRUE, levels = c("n",gsub("-",".",var_lst)))) %>%
        arrange(var_fac,cat) %>%
        group_by(var,var_lbl) %>%
        mutate(keep1row = row_number()) %>%
        ungroup %>%
        mutate(
          var = case_when(keep1row==1 ~ var,
                          TRUE ~ ""),
          var_lbl = case_when(keep1row==1 ~ var_lbl,
                        TRUE ~ ""),
          p.value = case_when(keep1row==1 ~ as.character(p.value),
                              TRUE ~ "")
        ) %>%
        dplyr::select(-var_fac,-keep1row) %>%
        kbl() %>% kable_material(c("striped", "hover"))
  }
  return(out)
}

# require(ROCR,pROC,PRROC)
get_perf_summ<-function(
  pred,
  real,
  keep_all_cutoffs=F,
  boots = 10
){
  bt<-data.frame(
    pred = pred,
    real = real
  )

  perf_at<-c()
  for (i in seq_len(boots)){
    # bootstrapping
    bti<-bt[sample(seq_len(nrow(bt)), round(0.8*nrow(bt)), replace=TRUE), ]
    pred<-bti$pred
    real<-bti$real

    # ROC object
    pred_obj<-ROCR::prediction(pred,real)
    prc<-ROCR::performance(pred_obj,"prec","rec")
    roc<-ROCR::performance(pred_obj,"sens","spec")
    nppv<-ROCR::performance(pred_obj,"ppv","npv")
    pcfall<-ROCR::performance(pred_obj,"pcfall")
    acc<-ROCR::performance(pred_obj,"acc")
    fscore<-ROCR::performance(pred_obj,"f")
    mcc<-ROCR::performance(pred_obj,"phi")

    # full performance metrics
    perf_at_i<-data.frame(
      cutoff=prc@alpha.values[[1]],
      prec=prc@y.values[[1]],
      rec_sens=prc@x.values[[1]],
      stringsAsFactors = F
    ) %>% 
    arrange(cutoff) %>%
    left_join(data.frame(
      cutoff=nppv@alpha.values[[1]],
      ppv=nppv@y.values[[1]],
      npv=nppv@x.values[[1]],
      stringsAsFactors = F
    ),
    by="cutoff") %>%
    mutate(prec_rec_dist=abs(prec-rec_sens)) %>%
    left_join(data.frame(
      cutoff=fscore@x.values[[1]],
      fscore=fscore@y.values[[1]],
      stringsAsFactors = F
    ),
    by="cutoff") %>%
    left_join(data.frame(
      cutoff=roc@alpha.values[[1]],
      spec=roc@x.values[[1]],
      stringsAsFactors = F
    ),
    by="cutoff") %>%
    mutate(
      Euclid_meas=sqrt((1-rec_sens)^2+(0-(1-spec))^2),
      Youden_meas=rec_sens+spec-1
    ) %>%
    left_join(data.frame(
      cutoff=pcfall@x.values[[1]],
      pcfall=pcfall@y.values[[1]],
      stringsAsFactors = F
    ),
    by="cutoff") %>%
    left_join(data.frame(
      cutoff=acc@x.values[[1]],
      acc=acc@y.values[[1]],
      stringsAsFactors = F
    ),
    by="cutoff") %>%
    left_join(data.frame(
      cutoff=mcc@x.values[[1]],
      mcc=mcc@y.values[[1]],
      stringsAsFactors = F
    ),
    by="cutoff") %>%
    filter(prec > 0 & rec_sens > 0 & spec > 0) %>%
    group_by(cutoff) %>%
    mutate(size=n()) %>%
    ungroup
  
  # performance summary
  lab1<-pred[real==1]
  lab0<-pred[real==0]
  pr<-PRROC::pr.curve(
    scores.class0 = lab1,
    scores.class1 = lab0,
    curve=F
  )
  roc_ci<-pROC::ci.auc(
    real,pred,
    levels = c(0, 1), direction = "<",
    quite = TRUE
  )
  
  perf_summ_i<-data.frame(
    overall_meas=c(
      "roauc_low",
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
      "opt_fscore"
    ),
    meas_val=c(
      roc_ci[[1]],
      roc_ci[[2]],
      roc_ci[[3]],
      perf_at_i$cutoff[which.min(perf_at_i$Euclid_meas)],
      perf_at_i$rec_sens[which.min(perf_at_i$Euclid_meas)],
      perf_at_i$spec[which.min(perf_at_i$Euclid_meas)],
      perf_at_i$ppv[which.min(perf_at_i$Euclid_meas)],
      perf_at_i$npv[which.min(perf_at_i$Euclid_meas)],
      pr$auc.integral,
      pr$auc.davis.goadrich,
      perf_at_i$prec[which.min(perf_at_i$prec_rec_dist)],
      perf_at_i$rec_sens[which.min(perf_at_i$prec_rec_dist)],
      perf_at_i$fscore[which.min(perf_at_i$prec_rec_dist)]
    ),
    stringsAsFactors = F) %>%
    bind_rows(
      perf_at_i %>% 
      summarize(
        prec_m = mean(prec,na.rm=T),
        sens_m = mean(rec_sens,na.rm=T),
        spec_m = mean(spec,na.rm=T),
        ppv_m = mean(ppv,na.rm=T),
        npv_m = mean(npv,na.rm=T),
        acc_m = mean(acc,na.rm=T),
        fscore_m = mean(fscore,na.rm=T),
        mcc_m = mean(mcc,na.rm=T),
        .groups = "drop"
      ) %>%
      pivot_longer(
        cols = everything(),
        names_to = "overall_meas",
        values_to = "meas_val"
      )
    ) %>% 
    mutate(bn = i)

    # stack
    perf_at %<>% bind_rows(perf_at_i) 
    perf_summ %<>% bind_rows(perf_summ_i)
  }

  # summary over bootstraps
  perf_summ %<>% 
    group_by(overall_meas) %>%
    summarise(
      meas_val_m  = median(meas_val,na.rm=T),
      meas_val_lb = quantile(meas_val,0.025,na.rm=T),
      meas_val_ub = quantile(meas_val,0.975,na.rm=T),
      .groups = "drop"  
    )
  out<-list(perf_summ=perf_summ)
  if(keep_all_cutoffs){
    perf_at %<>% 
      pivot_longer(
        cols = -c(cutoff),
        names_to = "meas",
        values_to = "meas_val"
      ) %>%
      group_by(cutoff,meas) %>% 
      summarise(
        meas_val_m = median(meas_val,na.rm=T), 
        meas_val_lb = quantile(meas_val,0.025,na.rm=T),
        meas_val_ub = quantile(meas_val,0.975,na.rm=T),
        .groups = "drop"
      )
    out$perf_at<-perf_at
  }
  return(out)
}

get_calibr<-function(
  pred,
  real,
  n_bin=20,
  test=TRUE
){
  # require("ResourceSelection")
  calib<-data.frame(pred=pred,
                    y=real) %>%
    arrange(pred) %>%
    dplyr::mutate(
      pred_bin = cut(
        pred,
        breaks=unique(quantile(pred,0:(n_bin)/(n_bin))),
        include.lowest=T,
        labels=F
      )
    ) %>%
    ungroup %>% group_by(pred_bin) %>%
    dplyr::summarize(
      expos=n(),
      bin_lower=min(pred),
      bin_upper=max(pred),
      bin_mid=median(pred),
      y_agg = sum(y),
      pred_p = mean(pred),
      .groups = "drop"
    ) %>%
    dplyr::mutate(y_p=y_agg/expos) %>%
    dplyr::mutate(
      binCI_lower = pmax(0,pred_p-1.96*sqrt(y_p*(1-y_p)/expos)),
      binCI_upper = pred_p+1.96*sqrt(y_p*(1-y_p)/expos)
    )

  out<-list(calib = calib)
  if(test){
    # brier score
    brier<-mean((real-pred)^2)
    # hosmer-lemeshow
    hl<-hoslem.test(real, pred, g = n_bin)
    # recalibration test
    fit<-lm(y_p~pred_p,data=calib)
    sfit<-summary(fit)
    t_b1<-(1-sfit$coefficients[2,1])/sfit$coefficients[2,2]
    pval_b1<- 2 * pt(abs(t_b1), df = df.residual(fit), lower.tail = FALSE)
    t_b0<-sfit$coefficients[1,3]
    pval_b0<-sfit$coefficients[1,4]
    # put together
    out[["test"]]<-data.frame(
      test = c('HL','Re-intx','Re-slope','Br'),
      statistics = c(hl$statistic,t_b1,t_b0,brier),
      pval = c(hl$p.value,pval_b1,pval_b0,NA)
    )
    rownames(out[["test"]])<-NULL
  }
  return(out)
}

# require(survival,SurvMetrics)
get_perf_summ.surv<-function(
  model,   # coxph model object
  data_ts, # testing data set
  time_col,
  status_col,
  eval_times = seq(90,90*4*3,by=90) # time of interests
){
  # cindex1 
  data_ts[,"pred"]<-predict(model,
                            newdata = data_ts,
                            type="expected")
  calib<-coxph(formula(paste0("Surv(",time_col,',',status_col,") ~ pred")),
               data = data_ts)
  cindex1<-summary(calib)$concordance
 
  # cindex2
  pred_med<-predict(model,
                    newdata = data_ts,
                    type="expected")
  surv_obj = Surv(unlist(data_ts[,time_col]), unlist(data_ts[,status_col]))
  cindex2<-Cindex(surv_obj, predicted = pred_med)

  # Brier score
  survfit_obj<-survfit(formula(paste0("Surv(",time_col,',',status_col,") ~ 1")),
                     data = data_ts)
  t_star<-survfit_obj$time[which.min(abs(survfit_obj$surv-0.5))]
  brier<-Brier(surv_obj, pre_sp = pred_med, t_star)

  # Integrated Brier score 
  mat_pred<-c()
  for(t in eval_times){
    data_ts[,time_col]<-t
    pred<-predict(model,
                  newdata=ts,
                  type="expected")
    mat_pred<-cbind(mat_pred,pred)
  }
  int_brier<-IBS(surv_obj, sp_matrix = mat_pred, eval_times)

  # return results
  out<-list(
    cindex1 = cindex1,
    cindex2 = cindex2,
    brier =  brier,
    int_brier = int_brier
  )
  return(out)
}

# require(SurvMetrics)
get_calibr.surv<-function(
  model,   # coxph model object
  data_ts, # testing data set
  time_col,
  status_col,
  eval_times = seq(90,90*4*3,by=90) # time of interests
){
  # observed 
  fit_frm<-formula(paste0("Surv(",time_col,",",status_col,") ~ 1"))
  surfit_obs<-summary(survfit(fit_frm,data=data_ts),times=eval_times)
  # predicted
  mat_pred<-c()
  for(t in eval_times){
    data_ts[,time_col]<-t
    pred<-predict(model,newdata=data_ts,
                  type="survival",se.fit=TRUE)
    mat_pred %<>%
      bind_rows(
        data.frame(
          t = t,
          pred = mean(pred$fit),
          se_pred = 0.5*sd(pred$fit)+0.5*mean(pred$se.fit)
        )
      )
  }
  # return results
  mat_pred %<>% mutate(obs = surfit_obs[["surv"]])
  return(mat_pred)
}

get_stab_summ<-function(
  rank_lst,
  var_colnm="var",
  rank_colnm="rk",
  metric=c("kuncheva",
          "wci",
          "wci_rel"),
  f=NULL,
  d=NULL
){
  K<-length(rank_lst)
  varlst<-c()
  rk_stack<-c()
  
  for(i in 1:K){
    #subset and reorder columns so that var_colnm is the 1st column
    varlst_i<-rank_lst[[i]][,c(var_colnm,rank_colnm)]
    #rename column for easy referring
    varlst_i<-setNames(varlst_i,c("var","rk"))
    #collect all distinct features
    varlst<-unique(c(varlst,varlst_i$var))
    #stack feature rankings
    rk_stack %<>% 
      bind_rows(varlst_i %>% 
                  mutate(mod_idx=i) %>%
                  mutate(p=nrow(varlst_i))) #add model index
  }
  
  if(metric=="kuncheva"){
    if(is.null(f)){
      f<-length(varlst)
    }
    
    if(is.null(d)){
      d<-round(f/2)
    }
    
    #kalousis and kuncheva index
    pairwise_fset<-c()
    for(i in 1:(K-1)){
      for(j in (i+1):K){
        ki_new<-rk_stack %>%
          filter(mod_idx %in% c(i,j)) %>%
          filter(rk<=d) %>%
          dplyr::select(var,mod_idx) %>%
          mutate(pair_idx=paste0(i,"_",j)) %>%
          mutate(mod_idx=ifelse(mod_idx==i,"ki","kj"),ref=1) %>%
          unique %>% spread(mod_idx,ref,fill=0) %>%
          mutate(inter=ki*kj,union=((ki+kj)>=1)*1)
        
        pairwise_fset %<>%
          bind_rows(ki_new)
      }
    }
    
    stb_idx<-pairwise_fset %>%
      group_by(pair_idx) %>%
      summarize(inter_cnt=sum(inter),
                union_cnt=sum(union),
                .groups="drop") %>%
      summarize(ki_unscale=sum((inter_cnt*f-d^2)/(d*(f-d)))*(2/(K*(K-1))),
                .groups="drop") %>%
      mutate(ki=(ki_unscale-(-1))/(1-(-1))) #scale index to 0,1 range
  }
  
  else if(metric=="wci"){
    #weighted consistency index
    stb_idx<-rk_stack %>%
      group_by(var) %>% 
      dplyr::summarize(phi_f=n(),.groups="drop") %>%
      dplyr::mutate(N=sum(phi_f)) %>%
      dplyr::mutate(phi_f_wt=(phi_f/N)*((phi_f-1)/(K-1))) %>%
      group_by(N) %>%
      dplyr::summarize(cw=sum(phi_f_wt),
                       .groups="drop") %>%
      select(cw) 
  }
  
  else if(metric=="wci_rel"){
    #(relative) weighted consistency index
    if(is.null(f)){
      f<-length(varlst)
    }
    
    stb_idx<-rk_stack %>%
      group_by(var) %>% 
      dplyr::summarize(phi_f=n(),
                       .groups="drop") %>%
      dplyr::mutate(N=sum(phi_f)) %>%
      ungroup %>%
      dplyr::mutate(phi_f_wt=(phi_f/N)*((phi_f-1)/(K-1))) %>%
      group_by(N) %>%
      dplyr::summarize(cw=sum(phi_f_wt),
                       .groups="drop") %>%
      mutate(D=N %% f,
             H=N %% K) %>%
      mutate(cw_min=(N^2-f*(N-D)-D^2)/(f*N*(K-1)),
             cw_max=(H^2+N*(K-1)-H*K)/(N*(K-1))) %>%
      mutate(cw_rel=(cw-cw_min)/(cw_max-cw_min)) %>%
      select(cw_rel) 
  }
  
  else{
    stop("stability index calculation method is not supported!")
  }
  
  return(stb_idx)
}

get_parity_summ<-function(
  pred,
  real,
  strata,
  n_bins = 20,
  boots_n = 10,
  verb = TRUE
){
  # requires(broom,rsample)
  N<-length(pred)
  dt<-data.frame(
    pred = pred,
    real = real,
    strata = strata
  ) %>%
    mutate(
      ns = sum(strata),
      ws = ns/N,
      nr = sum(real),
      wr = nr/N
    ) %>%
    arrange(pred)

  # unique pred
  pred_uni<-data.frame(pred = unique(dt$pred)) %>%
    mutate(
      pred_bin = cut(
        pred,
        breaks = quantile(pred,probs = seq(0,1,length.out=n_bins+1)),
        include.lowest = TRUE,
        labels = FALSE
      )
    ) %>%
    group_by(pred_bin) %>%
    mutate(
      lb=min(pred),
      ub=max(pred)
    ) %>%
    ungroup

  dt %<>%
    left_join(pred_uni,by="pred")

  # main procedure
  rslt<-data.frame(
    thresh = as.double(),
    summ_type = as.character(),
    summ_val = as.double()
  )
  for(b in seq_along(1:n_bins)){
    dt_sub<-dt %>%
      mutate(
        pred_ind = as.numeric(pred_bin>=b)
      )
    
    # bootstrap for CI
    dt_sub_boots<-bootstraps(dt_sub,times = boots_n)
    
    # overall accuracy
    dt_sub_boots_rslt1<-map(
      dt_sub_boots$splits,
      function(x){
        dat <- as.data.frame(x) %>%
          group_by(real) %>%
          summarize(
            pr = sum(pred_ind)/n(),
            nr = (n()-sum(pred_ind))/n(),
            .groups = "drop"
          ) %>% 
          pivot_longer(
            cols = c("pr","nr"),
            names_to = "summ_type",
            values_to = "summ_val"
          ) %>%
          unite(summ_type,c("summ_type","real")) %>%
          mutate(
            summ_type = recode(
              summ_type,
              pr_1 = "tpr",
              pr_0 = "fnr",
              nr_1 = "fpr",
              nr_0 = "tnr"
            ),
            thresh = b
          ) 
        }
      )
    
      # overall precision
      dt_sub_boots_rslt2<-map(
        dt_sub_boots$splits,
        function(x){
          dat <- as.data.frame(x) %>%
            group_by(pred_ind) %>%
            summarize(
              tr = sum(real)/n(),
              fr = (n()-sum(real))/n(),
              .groups = "drop"
            ) %>%
            pivot_longer(
              cols = c("tr","fr"),
              names_to = "summ_type",
              values_to = "summ_val"
            ) %>%
            unite(summ_type,c("summ_type","pred_ind")) %>%
            mutate(
              summ_type = recode(
                summ_type,
                tr_1 = "ppv",
                tr_0 = "fdr",
                fr_1 = "for",
                fr_0 = "npv"
              ),
              thresh = b
            )
          }
        )
      
      # accuracy disparity
      dt_sub_boots_rslt3<-map(
        dt_sub_boots$splits,
        function(x){
          dat <- as.data.frame(x) %>%
            group_by(real,strata,ws,wr) %>%
            summarize(
              pr = sum(pred_ind)/n(),
              nr = (n()-sum(pred_ind))/n(),
              .groups = "drop"
            ) %>%
            pivot_longer(
              cols = c("pr","nr"),
              names_to = "summ_type",
              values_to = "summ_val"
            ) %>%
            unite(summ_type,c("summ_type","real")) %>%
            mutate(
              summ_type = recode(
                summ_type,
                pr_1 = "tpr",
                pr_0 = "fnr",
                nr_1 = "fpr",
                nr_0 = "tnr"
              )
            ) %>%
            pivot_wider(
              names_from = "strata",
              values_from = "summ_val"
            ) %>%
            mutate(
              summ_val = `1` - `0`,
              summ_type = paste0('disp_',summ_type),
              thresh = b
            ) %>%
            select(thresh,summ_type,summ_val)
          }
        )
        
        # precision disparity
        dt_sub_boots_rslt4<-map(
          dt_sub_boots$splits,
          function(x){
            dat <- as.data.frame(x) %>%
              group_by(pred_ind,strata,ws,wr) %>%
                summarize(
                  tr = sum(real)/n(),
                  fr = (n()-sum(real))/n(),
                  .groups = "drop"
                ) %>%
                pivot_longer(
                  cols = c("tr","fr"),
                  names_to = "summ_type",
                  values_to = "summ_val"
                ) %>%
                unite(summ_type,c("summ_type","pred_ind")) %>%
                mutate(
                  summ_type = recode(
                    summ_type,
                    tr_1 = "ppv",
                    tr_0 = "fdr",
                    fr_1 = "for",
                    fr_0 = "npv"
                  )
                ) %>%
                pivot_wider(
                  names_from = "strata",
                  values_from = "summ_val"
                ) %>%
                mutate(
                  summ_val = `1` - `0`,
                  summ_type = paste0('disp_',summ_type),
                  thresh = b
                ) %>%
                select(thresh,summ_type,summ_val)
            }
          )
        
        # summarize over bootstrapped samples
        rslt %<>%
          bind_rows(
            bind_rows(dt_sub_boots_rslt1) %>%
              bind_rows(dt_sub_boots_rslt2) %>%
              bind_rows(dt_sub_boots_rslt3) %>%
              bind_rows(dt_sub_boots_rslt4) %>%
              group_by(summ_type,thresh) %>%
              summarise(
                summ_val_m = median(summ_val),
                summ_val_lb = quantile(summ_val,0.025),
                summ_val_ub = quantile(summ_val,0.975),
                .groups = "drop"
              )
          )
        
        # report progress
        if(verb){
          print(paste0("results generated for risk bin:",b))
        }
  }

  return(rslt)
}
