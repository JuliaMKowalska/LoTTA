
#' normalize continuous score function
#'
#' @param X - score data
#' @param trimmed - NULL by default, provide vector of score values to trim the data
#' @param normalize - whether data is to be normalized
#'
#' @return list with X - normalized score, d - parameter used to shift the score


normalize_cont_X <- function(X,trimmed=NULL,normalize=TRUE) {
  d=ifelse(normalize==TRUE,(min(X)+max(X))*0.5,0)
  s=ifelse(normalize==TRUE,max(X)-min(X),1)
  if(is.null(trimmed)==FALSE){
    a=min(trimmed)
    b=max(trimmed)
    X=X[X >= a & X<= b]
  }
  X_norm=(X-d)/s
  return(list(X=X_norm,d=d,s=s))
}


#' normalize continuous outcome function
#'
#' @param Y - outcome data
#' @param X - score data
#' @param trimmed - NULL by default, provide vector of values to trim the data
#' @param normalize - whether data is to be normalized
#'
#' @return list with Y - normalized outcome, mu - shift parameter, sd - scaling parameter


normalize_cont_Y <- function(Y,X,trimmed=NULL,normalize=TRUE) {
  mu=ifelse(normalize==TRUE,mean(Y),0)
  sd=ifelse(normalize==TRUE,sd(Y),1)
  Y_norm=(Y-mu)/sd
  if(is.null(trimmed)==FALSE){
    a=min(trimmed)
    b=max(trimmed)

    Y_norm=Y_norm[X >= a & X <= b]
  }
  return(list(Y=Y_norm,mu=mu,sd=sd))
}

#' normalize discrete score function
#'
#' @param X score data
#' @param trimmed - NULL by default, provide vector of values to trim the data
#' @param normalize - TRUE by default, whether data is to be normalized
#'
#' @return list with X - normalized score, s - parameter used to scale the score


normalize_dis_X <- function(X,trimmed=NULL,normalize=TRUE) {
  M=as.integer(max(abs(X)))
  if(normalize==TRUE){
    sc1=10^(-nchar(as.character(M))+1)
    sc2=10^(-nchar(as.character(M)))
    if(abs(max(abs(X*sc1))-1)<abs(max(abs(X*sc2))-1)){
      sc=sc1
    }
    else{
      sc=sc2
    }
  }
  else{
    sc=1
  }
  if(is.null(trimmed)==FALSE){
    a=trimmed[1]
    b=trimmed[2]
    X=X[X >= a & X<= b]
  }
  X_norm=X*sc

  return(list(X=X_norm,d=0,s=sc^{-1}))
}
#' Binary outcomes for trimmed score
#'
#' @param Y - outcome data
#' @param X - score data
#' @param trimmed - NULL by default, provide vector of values to trim the data
#'
#' @return list with Y - outcomes for trimmed score

trim_dis_Y<-function(Y,X,trimmed=NULL){
  Y_trim=Y
  if(is.null(trimmed)==FALSE){
    a=trimmed[1]
    b=trimmed[2]
    Y_trim=Y[X >= a & X<= b]
  }
  return(list(Y=Y_trim))
}
