
#' normalize continuous score function
#'
#' @param x - score data
#' @param trimmed - NULL by default, provide vector of score values to trim the data
#' @param normalize - whether data is to be normalized
#'
#' @return list with x - normalized score, d - parameter used to shift the score


normalize_cont_x <- function(x,trimmed=NULL,normalize=TRUE) {
  d=ifelse(normalize==TRUE,(min(x)+max(x))*0.5,0)
  s=ifelse(normalize==TRUE,max(x)-min(x),1)
  if(is.null(trimmed)==FALSE){
    a=min(trimmed)
    b=max(trimmed)
    x=x[x >= a & x<= b]
  }
  x_norm=(x-d)/s
  return(list(x=x_norm,d=d,s=s))
}


#' normalize continuous outcome function
#'
#' @param y - outcome data
#' @param x - score data
#' @param trimmed - NULL by default, provide vector of values to trim the data
#' @param normalize - whether data is to be normalized
#'
#' @return list with y - normalized outcome, mu - shift parameter, sd - scaling parameter


normalize_cont_y <- function(y,x,trimmed=NULL,normalize=TRUE) {
  mu=ifelse(normalize==TRUE,mean(y),0)
  sd=ifelse(normalize==TRUE,sd(y),1)
  y_norm=(y-mu)/sd
  if(is.null(trimmed)==FALSE){
    a=min(trimmed)
    b=max(trimmed)

    y_norm=y_norm[x >= a & x <= b]
  }
  return(list(y=y_norm,mu=mu,sd=sd))
}

#' normalize discrete score function
#'
#' @param x score data
#' @param trimmed - NULL by default, provide vector of values to trim the data
#' @param normalize - TRUE by default, whether data is to be normalized
#'
#' @return list with x - normalized score, s - parameter used to scale the score


normalize_dis_x <- function(x,trimmed=NULL,normalize=TRUE) {
  M=as.integer(max(abs(x)))
  if(normalize==TRUE){
    sc1=10^(-nchar(as.character(M))+1)
    sc2=10^(-nchar(as.character(M)))
    if(abs(max(abs(x*sc1))-1)<abs(max(abs(x*sc2))-1)){
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
    x=x[x >= a & x<= b]
  }
  x_norm=x*sc

  return(list(x=x_norm,d=0,s=sc^{-1}))
}
#' Binary outcomes for trimmed score
#'
#' @param y - outcome data
#' @param x - score data
#' @param trimmed - NULL by default, provide vector of values to trim the data
#'
#' @return list with y - outcomes for trimmed score

trim_dis_y<-function(y,x,trimmed=NULL){
  y_trim=y
  if(is.null(trimmed)==FALSE){
    a=trimmed[1]
    b=trimmed[2]
    y_trim=y[x >= a & x<= b]
  }
  return(list(y=y_trim))
}
