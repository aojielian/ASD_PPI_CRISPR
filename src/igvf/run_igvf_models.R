#!/usr/bin/env Rscript
# Frozen custom Firth logistic implementation, adapted from the R1.1 runner.
# Input rows and covariates must already be the locked complete-case designs.
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=3L) stop("Usage: Rscript run_igvf_models.R INPUT_623.tsv INPUT_610.tsv OUTPUT_DIR")
inputs <- args[1:2]; out <- args[3]
dir.create(out, recursive=TRUE, showWarnings=FALSE)
log1pexp <- function(x) pmax(x,0)+log1p(exp(-abs(x)))
firth_eval <- function(beta,X,y,want_score=TRUE) {
  eta<-as.vector(X%*%beta); pr<-plogis(eta); w<-pmax(pr*(1-pr),1e-12)
  info<-crossprod(X,X*w); deti<-try(determinant(info,logarithm=TRUE),silent=TRUE)
  if(inherits(deti,"try-error")||deti$sign<=0||!is.finite(as.numeric(deti$modulus))) return(list(ll=-Inf,score=rep(NA_real_,ncol(X))))
  ll<-sum(y*eta-log1pexp(eta))+.5*as.numeric(deti$modulus)
  if(!want_score)return(list(ll=ll,score=NULL))
  inv<-try(solve(info),silent=TRUE); if(inherits(inv,"try-error"))return(list(ll=-Inf,score=rep(NA_real_,ncol(X))))
  h<-w*rowSums((X%*%inv)*X)
  list(ll=ll,score=as.vector(crossprod(X,y-pr+h*(.5-pr))))
}
fit_firth <- function(X,y) {
  X<-as.matrix(X); storage.mode(X)<-"double"; y<-as.numeric(y); q<-qr(X)
  if(q$rank<ncol(X))stop("Rank-deficient frozen design")
  op<-optim(rep(0,ncol(X)),fn=function(b)-firth_eval(b,X,y)$ll,gr=function(b)-firth_eval(b,X,y)$score,
             method="BFGS",control=list(maxit=3000,reltol=1e-10))
  if(op$convergence!=0L||!is.finite(op$value))stop("Firth optimization did not converge")
  list(beta=op$par,maxll=-op$value,X=X,y=y)
}
profile_exposure <- function(fit) {
  X<-fit$X; y<-fit$y; idx<-match("exposure",colnames(X)); if(is.na(idx))stop("Missing exposure term")
  pll<-function(v) {
    start<-fit$beta[-idx]
    obj<-function(u){b<-numeric(ncol(X));b[idx]<-v;b[-idx]<-u;-firth_eval(b,X,y,FALSE)$ll}
    gr<-function(u){b<-numeric(ncol(X));b[idx]<-v;b[-idx]<-u;-firth_eval(b,X,y)$score[-idx]}
    z<-optim(start,obj,gr,method="BFGS",control=list(maxit=2000,reltol=1e-9)); if(!is.finite(z$value))-Inf else -z$value
  }
  ll0<-pll(0); p<-pchisq(max(0,2*(fit$maxll-ll0)),1,lower.tail=FALSE); cut<-qchisq(.95,1); center<-fit$beta[idx]
  bound<-function(dir){step<-1;edge<-center+dir*step;v<-2*(fit$maxll-pll(edge));while(is.finite(v)&&v<cut&&step<40){step<-step*1.6;edge<-center+dir*step;v<-2*(fit$maxll-pll(edge))};if(!is.finite(v)||v<cut)return(NA_real_);uniroot(function(z)2*(fit$maxll-pll(z))-cut,min(center,edge),max(center,edge),tol=1e-6)$root}
  c(beta=center,lo=bound(-1),hi=bound(1),p=p)
}
forms<-c("y ~ exposure + source_category_model + log_degree + abundance_293t + bioplex_bait_status + log_promoters",
         "y ~ exposure + source_category_model + log_degree + abundance_293t + bioplex_bait_status + log_promoters + loeuf_z")
results<-list()
for(i in seq_along(inputs)){
 d<-fread(inputs[[i]],na.strings=c("NA",""),check.names=FALSE)
 expected<-if(i==1L)623L else 610L
 if(nrow(d)!=expected||anyDuplicated(d$gene)||!all(c("gene","exposure","y","source_category_model","log_degree","abundance_293t","bioplex_bait_status","log_promoters")%in%names(d)))stop("Frozen IGVF input schema/population mismatch")
 if(i==2L&&!"loeuf_z"%in%names(d))stop("LOEUF model input lacks loeuf_z")
 d[, exposure:=as.integer(exposure)]; d[, y:=as.integer(y)]; d[, source_category_model:=factor(source_category_model)]
 d[, bioplex_bait_status:=as.integer(bioplex_bait_status)]
 if(anyNA(d))stop("Missing data in the supplied complete-case population")
 X<-model.matrix(as.formula(forms[[i]]),as.data.frame(d)); fit<-fit_firth(X,d$y); pr<-profile_exposure(fit)
 results[[i]]<-data.table(model=if(i==1L)"IGVF_623_FIRTH" else "IGVF_610_LOEUF_FIRTH",status="COMPLETE",N=nrow(d),events=sum(d$y),prey_N=sum(d$exposure==1L),comparator_N=sum(d$exposure==0L),beta_exposure=pr["beta"],OR=exp(pr["beta"]),CI_low_OR=exp(pr["lo"]),CI_high_OR=exp(pr["hi"]),P_raw=pr["p"],estimator="Custom Jeffreys-penalized likelihood; profile-likelihood CI and LR P",formula=forms[[i]])
}
fwrite(rbindlist(results),file.path(out,"IGVF_FIRTH_MODEL_RESULTS.tsv"),sep="\t",na="NA",quote=FALSE)
