bestfit<- function(data, dist , order=NULL, start=NULL, conf=0.95){
  if(!is.numeric(data)|| is.character(data) || is.matrix(data))
    stop("data must be a numeric vector")
  distname<- dist
  pdistname<-paste("p", distname, sep="")
  if (!exists(pdistname, mode="function"))
    stop(paste("The ", pdistname, " function must be defined"))
  if (is.element(distname, c("binom", "nbinom", "geom", "hyper", "pois"))){
    stop(paste("The distribution must be continuous"))}
  memp<- function(x, order){
    mean(x^order)
  }
  if (is.element(distname, c("exp", "norm", "lnorm", "gamma", "logis", "weibull", "unif"))){
    fit1<-fitdistrplus::fitdist(data, distname)
    fit2<-fitdistrplus::fitdist(data, distname, method="mme", order=c(1,2), memp=memp)
    parameters1<-coef(fit1)
    parameters2<-coef(fit2)
    ks1<-do.call(ks.test, c(list(data), as.list(pdistname), as.list(parameters1)))
    ks2<-do.call(ks.test, c(list(data), as.list(pdistname), as.list(parameters2)))
  }
  else {
    if(is.null(order)){
      stop("moment matching estimation needs an 'order' argument")}
    else {
      fit1<-fitdistrplus::fitdist(data, distname, start=start)
      fit2<-fitdistrplus::fitdist(data, distname, method="mme", order=order, memp=memp, start=start)
      parameters1<-coef(fit1)
      parameters2<-coef(fit2)
      ks1<-do.call(ks.test, c(list(data), as.list(pdistname), as.list(parameters1)))
      ks2<-do.call(ks.test, c(list(data), as.list(pdistname), as.list(parameters2)))}
  }
  pvalue1<-ks1$p.value
  pvalue2<-ks2$p.value
  if(!is.na(pvalue1) & !is.na(pvalue2)){
    if(pvalue1>=pvalue2 & pvalue1 <(1- conf)){stop(paste( distname , "distribution is not the best"))}
    else if(pvalue1>=pvalue2 & pvalue1>=(1-conf)) {par<-parameters1
    pvalue<-pvalue1
    fit<-fit1}
    else {par<-parameters2
    pvalue<-pvalue2
    fit<-fit2}
  }
  else if(!is.na(pvalue1) & is.na(pvalue2)){
    if(pvalue1>=(1-conf)){par<-parameters1
    pvalue<-pvalue1
    fit<-fit1}
    else {stop(paste( distname , "distribution is not the best"))}
  }
  else if(is.na(pvalue1) & !is.na(pvalue2)){
    if(pvalue2>=(1-conf)){par<-parameters2
    pvalue<-pvalue2
    fit<-fit2}
    else {stop(paste( distname , "distribution is not the best"))}
  }
  else{stop(paste( "The function failed to estimate the parameters"))}
  names(pvalue)<-c("pvalue")
  names(fit)<-c("fit")
  p.value<-pvalue[c("pvalue")]
  fit<-fit[c("fit")]
  re<-c(fit, p.value)
  return(re)
}
