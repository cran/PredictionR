predmix <- function(data, s, n, parameters,  conf=0.95){
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
  if(!is.numeric(parameters)|| is.character(parameters) || is.matrix(parameters))
    stop("parameters must be a numeric vector")
  if(length(parameters) == 0 )
    stop("the length of parameters must be greater than or equal to 1")
  distname<- "mixexp2"
  pdistname<-paste("p", distname, sep="")
  qdistname<-paste("q", distname, sep="")
  usolve=function(r,s){
    comb=function(l,m){
      return( factorial(l) / (factorial(m) * factorial(l-m)) )
    }
    n<-n
    r<-r
    s<-s
    j<-0:(s-r-1)
    S<-numeric(1)
    a<-factorial(n)/(factorial(r-1)*factorial(n-s)*factorial(s-r-1))
    eq<-function(x){
      for(i in 0:(r-1)){
        S<-S+sum((-1)^(i+j)*comb(r-1,i)*comb(s-r-1,j)*((n-s+j+1)*((n-r+i+1)+(x*(n-s+j+1))))^(-1))}
      return((1-conf)-a*S)}
    uniroot(eq,c(0,1),extendInt = "yes")
  }

  r<-length(data)
  x<-sort(data)
  q1<-(s)/(n+1)
  q2<-(r)/(n+1)
  q4<-(s-3)/(n+1)
  root1<-usolve(r,s)$root
  if((s-r)>3){
    root2<-usolve(r,s-3)$root
  }
  pdist<- do.call(pdistname, c(x[r], as.list(parameters)))
  ptot1<- (1-(1-pdist)^(root1+1))
  if((s-r)>3){
    ptot2<- (1-(1-pdist)^(root2+1))
  }
  upper1<- do.call(qdistname, c(ptot1, as.list(parameters)))
  if((s-r)>3){
    upper2<- do.call(qdistname, c(ptot2, as.list(parameters)))
    c2=(do.call(qdistname, c(q4, as.list(parameters)))-do.call(qdistname, c(q2, as.list(parameters))))/(upper2-do.call(qdistname, c(q2, as.list(parameters))))
    pointest2=x[r]+c2*(upper2-x[r])
  }
  c=(do.call(qdistname, c(q1, as.list(parameters)))-do.call(qdistname, c(q2, as.list(parameters))))/(upper1-do.call(qdistname, c(q2, as.list(parameters))))
  pointest=x[r]+c*(upper1-x[r])
  if((s-r)<=3){lower<-x[r]
  }
  else{
    lower<-pointest2
  }
  if(lower>= upper1){stop("Lower bound can not be greater than upper bound")}
  interval<-c(lower,upper1)
  point<- pointest
  names(interval)<-c("lower","upper")
  names(point)<-c("predicted point")
  interval<-interval[c("lower", "upper")]
  point<-point[c("predicted point")]
  int<-list(interval=interval, point=point, lower=lower, upper=upper1)
  return(structure(int, class = "predmix"))
}
print.predmix <- function(x, ...) {
  if (!inherits(x, "predmix"))
    stop("Use only with 'predmix' objects")
  cat("Prediction  for the future observation based on mixture of two exponential distribution \n")
  cat("Point:\n")
  print(cbind.data.frame("predicted point" = x$point), ...)
  cat("Interval:\n")
  print(cbind.data.frame("PI" = x$interval), ...)
}
