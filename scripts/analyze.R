options(warn=1)
.libPaths(c('software/Rlib',.libPaths()))
library(jsonlite)
source('software/TwoSampleMR_mr.R')
source('software/TwoSampleMR_mode.R')
dir.create('results',showWarnings=FALSE)
normal_result <- function(b,se) list(beta=b,se=se,p=2*pnorm(-abs(b/se)),ci_low=b-qnorm(.975)*se,ci_high=b+qnorm(.975)*se,OR=exp(b),OR_low=exp(b-qnorm(.975)*se),OR_high=exp(b+qnorm(.975)*se))
ivw <- function(d) {
 w<-1/d$se_outcome^2; x<-d$beta; y<-d$beta_outcome; b<-sum(w*x*y)/sum(w*x*x)
 q<-sum(w*(y-b*x)^2); df<-nrow(d)-1; se<-sqrt(max(1,q/df)/sum(w*x*x))
 c(normal_result(b,se),list(Q=q,Q_df=df,Q_p=pchisq(q,df,lower.tail=FALSE),phi=max(1,q/df)))
}
egger <- function(d,orient=TRUE,underdispersion=FALSE) {
 x<-d$beta; y<-d$beta_outcome
 if(orient) {y<-y*sign(x); x<-abs(x)}
 fit<-lm(y~x,weights=1/d$se_outcome^2); s<-summary(fit); cf<-coef(s); correction<-if(underdispersion) 1 else min(1,s$sigma)
 se<-cf[2,2]/correction; sei<-cf[1,2]/correction; df<-nrow(d)-2
 list(beta=cf[2,1],se=se,p=2*pt(-abs(cf[2,1]/se),df),ci_low=cf[2,1]-qt(.975,df)*se,ci_high=cf[2,1]+qt(.975,df)*se,intercept=cf[1,1],intercept_se=sei,intercept_p=2*pt(-abs(cf[1,1]/sei),df),residual_Q=sum(resid(fit)^2/d$se_outcome^2),df=df)
}
median_est <- function(d,second=FALSE,B=10000) {
 x<-d$beta; y<-d$beta_outcome; sx<-d$se; sy<-d$se_outcome
 w<-x*x/(sy*sy); if(second) w<-1/(sy^2/x^2+y^2*sx^2/x^4)
 b<-weighted_median(y/x,w)
 set.seed(20260908); se<-weighted_median_bootstrap(x,y,sx,sy,w,B)
 normal_result(b,se)
}
orient_data <- function(d) {
 d$egger_reoriented<-d$beta<0; flip<-d$egger_reoriented
 d$beta[flip]<- -d$beta[flip]; d$beta_outcome[flip]<- -d$beta_outcome[flip]
 a<-d$effect_allele[flip]; d$effect_allele[flip]<-d$other_allele[flip]; d$other_allele[flip]<-a
 d$eaf[flip]<-1-d$eaf[flip]; d$eaf_outcome_aligned[flip]<-1-d$eaf_outcome_aligned[flip]
 d
}
run_one <- function(path,key,optional=TRUE) {
 d<-read.csv(path); stopifnot(nrow(d)>2,all(d$se>0),all(d$se_outcome>0),!anyDuplicated(d$rsid))
 if('palindromic'%in%names(d)) d$palindromic<-tolower(as.character(d$palindromic))=='true'
 oriented<-orient_data(d); write.csv(oriented,paste0('results/',key,'_egger_oriented.csv'),row.names=FALSE)
 stopifnot(all(oriented$beta>0),all(oriented$se==d$se),all(oriented$se_outcome==d$se_outcome))
 v<-ivw(d); e<-egger(d); wm<-median_est(d); wm2<-median_est(d,TRUE)
 # Independent package-source implementations cross-check core point estimates and uncertainty.
 pkgv<-mr_ivw(d$beta,d$beta_outcome,d$se,d$se_outcome)
 pkge<-mr_egger_regression(d$beta,d$beta_outcome,d$se,d$se_outcome,list())
 stopifnot(abs(pkgv$b-v$beta)<1e-10,abs(pkgv$se-v$se)<1e-10,abs(pkge$b-e$beta)<1e-10,abs(pkge$se-e$se)<1e-10)
 stopifnot(abs(ivw(oriented)$beta-v$beta)<1e-10,abs(median_est(oriented)$beta-wm$beta)<1e-10)
 x<-abs(d$beta); wx<-1/d$se^2; qgx<-sum(wx*(x-weighted.mean(x,wx))^2); i2gx<-max(0,1-(nrow(d)-1)/qgx)
 f<-(d$beta/d$se)^2
 loo<-do.call(rbind,lapply(seq_len(nrow(d)),function(i) data.frame(left_out=d$rsid[i],as.data.frame(ivw(d[-i,])))))
 write.csv(loo,paste0('results/',key,'_loo.csv'),row.names=FALSE)
 w<-d$beta^2/d$se_outcome^2; h<-w/sum(w); residual<-(d$beta_outcome-v$beta*d$beta)/d$se_outcome
 rad<-data.frame(rsid=d$rsid,radial_x=abs(d$beta)/d$se_outcome,radial_y=d$beta_outcome*sign(d$beta)/d$se_outcome,Q_contribution=residual^2,leverage=h,studentized=residual/sqrt(1-h))
 rad$p<-pchisq(rad$studentized^2,1,lower.tail=FALSE); rad$p_bonferroni<-pmin(1,rad$p*nrow(d)); rad$flag<-rad$p_bonferroni<.05
 write.csv(rad,paste0('results/',key,'_radial.csv'),row.names=FALSE)
 # Q/Anderson-Rubin test uses both SNP-exposure and SNP-outcome uncertainty, df=K.
 arq<-function(b) sum((d$beta_outcome-b*d$beta)^2/(d$se_outcome^2+b*b*d$se^2))
 crit<-qchisq(.95,nrow(d)); grid<-sort(unique(c(seq(-100,100,length.out=20001),-exp(seq(log(100),log(1e8),length.out=1000)),exp(seq(log(100),log(1e8),length.out=1000))))); val<-vapply(grid,function(b)arq(b)-crit,0.0)
 crosses<-which(val[-length(val)]*val[-1]<0)
 roots<-vapply(crosses,function(i)uniroot(function(b) arq(b)-crit,grid[c(i,i+1)],tol=1e-10)$root,0.0)
 endpoints<-c(-Inf,roots,Inf); intervals<-list()
 for(i in seq_len(length(endpoints)-1)) {
  lo<-endpoints[i]; hi<-endpoints[i+1]; mid<-if(!is.finite(lo)) hi-1e9 else if(!is.finite(hi)) lo+1e9 else (lo+hi)/2
  if(arq(mid)<=crit) intervals[[length(intervals)+1]]<-c(lo,hi)
 }
 out<-list(n=nrow(d),exposure_N=unique(d$n),F=list(mean=mean(f),median=median(f),min=min(f),max=max(f)),ivw=v,weighted_median=wm,weighted_median_second_order=wm2,egger=e,I2GX=i2gx,QGX=qgx,oriented_snps=sum(d$beta<0),loo=list(beta_min=min(loo$beta),beta_max=max(loo$beta),p_min=min(loo$p),p_max=max(loo$p),n_nominal=sum(loo$p<.05)),radial=list(n_outliers=sum(rad$flag),outliers=rad$rsid[rad$flag],min_bonferroni_p=min(rad$p_bonferroni)),anderson_rubin=list(Q_null=arq(0),df=nrow(d),p_null=pchisq(arq(0),nrow(d),lower.tail=FALSE),confidence_set=intervals),steiger=list(status='Not run: population prevalence of each occupational category and sex-specific PD prevalence are not supplied; binary exposure is on a linear-probability scale. No unlabelled comparison of observed-exposure and liability-outcome R2.'))
 if(any(rad$flag) && sum(!rad$flag)>=3) out$radial_removed_ivw<-ivw(d[!rad$flag,])
 if('palindromic'%in%names(d) && sum(!d$palindromic)>=3) out$all_palindromes_removed<-c(list(n=sum(!d$palindromic)),ivw(d[!d$palindromic,]))
 if(optional) {
  cat(key,'weighted mode started\n'); flush.console()
  set.seed(20260908); mp<-default_parameters(); mp$nboot<-10000; mp$phi<-1
  mo<-tryCatch(mr_weighted_mode(d$beta,d$beta_outcome,d$se,d$se_outcome,mp),error=function(e)list(error=conditionMessage(e)))
  out$weighted_mode<-mo
  for(loss in c('huber','l2')) {
   warnings<-character()
   rr<-withCallingHandlers(tryCatch(mr.raps::mr.raps(data.frame(beta.exposure=d$beta,beta.outcome=d$beta_outcome,se.exposure=d$se,se.outcome=d$se_outcome),diagnostics=FALSE,over.dispersion=TRUE,loss.function=loss),error=function(e) list(error=conditionMessage(e))),warning=function(w){warnings<<-c(warnings,conditionMessage(w)); invokeRestart('muffleWarning')})
   out[[paste0('raps_',loss)]]<-unclass(rr); out[[paste0('raps_',loss,'_warnings')]]<-warnings
  }
  cat(key,'MR-PRESSO started\n'); flush.console()
  warnings<-character()
  pr<-withCallingHandlers(tryCatch(MRPRESSO::mr_presso(BetaOutcome='beta_outcome',BetaExposure='beta',SdOutcome='se_outcome',SdExposure='se',data=d,OUTLIERtest=TRUE,DISTORTIONtest=TRUE,SignifThreshold=.05,NbDistribution=10000,seed=20260908),error=function(e)list(error=conditionMessage(e))),warning=function(w){warnings<<-c(warnings,conditionMessage(w)); invokeRestart('muffleWarning')})
  out$presso<-pr; out$presso_warnings<-warnings
 }
 saveRDS(out,paste0('results/',key,'_analysis.rds'))
 write_json(out,paste0('results/',key,'_analysis.json'),auto_unbox=TRUE,pretty=TRUE,digits=16,na='null',null='null')
 cat(key,'IVW',v$beta,v$se,v$p,'Egger',e$beta,e$p,'intercept p',e$intercept_p,'I2GX',i2gx,'\n')
 out
}
args<-commandArgs(TRUE)
if('saved'%in%args) {
 audit<-list()
 for(sex in c('female','male')) {
  d<-read.csv(paste0(sex,'_harmonized_main.csv'))
  audit[[sex]]<-list(previous_algorithm_recomputed=egger(d,FALSE,TRUE),orientation_only_recomputed=egger(d,TRUE,TRUE),corrected_egger=egger(d),analysis=run_one(paste0(sex,'_harmonized_main.csv'),paste0('saved_primary_',sex),FALSE))
 }
 write_json(audit,'results/saved_primary_correction_audit.json',auto_unbox=TRUE,pretty=TRUE,digits=16)
} else {
 all<-list()
 for(code in c('22614_1','22614_0','22614_2','22610_2')) for(sex in c('female','male')) {
  key<-paste(code,sex,sep='_'); if(length(args)>0 && !(key%in%args)) next
  all[[key]]<-run_one(paste0('results/',key,'_harmonized.csv'),key)
 }
 if(length(args)==0) write_json(all,'results/all_analyses.json',auto_unbox=TRUE,pretty=TRUE,digits=16,na='null',null='null')
}
capture.output(sessionInfo(),file=if(length(args)>0) paste0('results/R_sessionInfo_',args[1],'.txt') else 'results/R_sessionInfo.txt')
