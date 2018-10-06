predI <- function(data, dist, s, n, order=NULL, start=NULL, conf=0.95){
  if(!is.numeric(data)|| is.character(data) || is.matrix(data))
    stop("data must be a numeric vector")
  if(!is.numeric(s)|| is.character(s) || is.matrix(s) )
    stop("s must be a numeric vector")
  if(length(s) != 1 )
    stop("the length of s must be equal to 1")
  if(!is.numeric(n)|| is.character(n) || is.matrix(n) )
    stop("n must be a numeric vector")
  if( length(data)>s) { stop("s must be greater than data length")}
  if( length(data)>n) { stop("n must be greater than data length")}
  if( s>n) { stop("n must be greater than s")}
  distname<- dist
  pdistname<-paste("p", distname, sep="")
  qdistname<-paste("q", distname, sep="")
  if (!exists(pdistname, mode="function"))
    stop(paste("The ", pdistname, " function must be defined"))
  if (!exists(qdistname, mode="function"))
    stop(paste("The ", qdistname, " function must be defined"))
  if (is.element(distname, c("exp", "norm", "lnorm", "gamma", "logis", "weibull", "unif"))){
    fitting<-bestfit(data, distname, conf=conf)
    parameters<-fitting$fit
  }
  else {
    if(is.null(order)){
      stop("moment matching estimation needs an 'order' argument")}
    else {
      fitting<-bestfit(data, distname, order=order, start=start, conf=conf)
      parameters<-fitting$fit}
  }
  comb<-function(l,m){
    return(factorial(l)/(factorial(m)*factorial(l-m)))
  }
  r<-length(data)
  x<-sort(data)
  j<-0:(s-r-1)
  S<-numeric(1)
  a=factorial(n)/(factorial(r-1)*factorial(n-s)*factorial(s-r-1))
  eq=function(f){
    for(i in 0:(r-1)){
      S=S+sum((-1)^(i+j)*comb(r-1,i)*comb(s-r-1,j)*((n-s+j+1)*((n-r+i+1)+(f*(n-s+j+1))))^(-1))}
    return((1-conf)-a*S)}
  sol<-uniroot(eq,c(0,1),extendInt = "yes")
  root<-sol$root
  lower<-x[r]
  pdist<- do.call(pdistname, c(x[r], as.list(parameters)))
  ptot<- (1-(1-pdist)^(root+1))
  upper<- do.call(qdistname, c(ptot, as.list(parameters)))
  if(lower>= upper){stop("Lower bound can not be greater than upper")}
  interval<-c(lower,upper)
  names(interval)<-c("lower","upper")
  interval<-interval[c("lower", "upper")]
  int<-list(interval=interval, lower=lower, upper=upper, distname=distname, r=r, s=s, n=n, parameters=parameters)
  return(structure(int, class = "predI"))
}
print.predI <- function(x, ...) {
  if (!inherits(x, "predI"))
    stop("Use only with 'predI' objects")
  cat("Prediction interval for the next observation based on '", x$distname, "' distribution \n")
  cat("Parameters:\n")
  print(cbind.data.frame("etsimate" = x$parameters), ...)
  cat("Interval:\n")
  print(cbind.data.frame("PCI" = x$interval), ...)
}
