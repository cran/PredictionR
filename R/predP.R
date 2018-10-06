predP<-function(data, dist, no, n, order=NULL, start=NULL, conf=0.95){
  if(!is.numeric(data)|| is.character(data) || is.matrix(data))
    stop("data must be a numeric vector")
  if(!is.numeric(no)|| is.character(no) || is.matrix(no) )
    stop("no must be a numeric vector")
  if(!is.numeric(n)|| is.character(n) || is.matrix(n) )
    stop("n must be a numeric vector")
  distname<- dist
  ld<-length(data)
  ls<-no+ld
  if(ls > n ){stop("number of observations exceeds n")}
  if( ld > n) { stop("n must be greater than data length")}
  if (is.element(distname, c("exp", "norm", "lnorm", "gamma", "logis", "weibull", "unif"))){
    for(k in (ld+1):ls){
      s<-k
      PI<-predI(data, distname, s, n, conf=conf)
      l<-PI$lower
      u<-PI$upper
      xnew<- (l+u)/2
      data<-c(data,xnew)}
  }
  else{
    if(is.null(order)){
      stop("moment matching estimation needs an 'order' argument")}
    else {
      for(k in (ld+1):ls){
        s<-k
        PI<-predI(data, distname, s, n, order=order, start=start, conf=conf)
        l<-PI$lower
        u<-PI$upper
        xnew<- (l+u)/2
        data<-c(data,xnew)}
    }
  }
  ns<- (ld+1):ls
  newobs<- data[ns]
  results <-list(data=data,  newobs=newobs , ns=ns , no=no , distname=distname , ld=ld  ,  n=n)
  return(structure(results, class = "predP"))
}
print.predP <- function(x,...) {
  if (!inherits(x, "predP"))
    stop("Use only with 'predP' objects")
  cat("Prediction points for xs based on '", x$distname, "' distribution \n")
  cat("Prediction points:\n")
  print(cbind.data.frame("s" = x$ns, "xs" = x$newobs), ...)
}
