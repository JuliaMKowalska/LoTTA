#' logit function
#' @param x - score data
#' @return value of logit function at x
logit<-function(x){
  return(log(x/(1-x)))
}

#' inverse logit function
#' @param x - score data
#' @return value of inverse logit function at x
invlogit<-function(x){
  return(1/(1+exp(-x)))
}

#' function that finds maximum widow size to searxh for a cutoff
#' @param x - score data
#' @param ns - minimum number of data points on each side of the cutoff to which cubic
#' parts are fitted
#' @return list with ubl - minimum value of the window's left boundary point, ubr - maximum value of the window's right boundary point
# lb - minimum window size (grid size in case of discrete score)


bounds<-function(x,ns=25){
  xu=sort(unique(x))
  ql=ns/length(x)
  q=as.numeric(quantile(x,c(ql,1-ql)))
  ubr=q[2]
  ubl=q[1]
  N=length(xu)
  diff=xu[2:N]-xu[1:(N-1)]
  diff1=diff[1:(N-2)]
  diff2=diff[2:(N-1)]
  Diff=ifelse(diff1>diff2,diff1,diff2)
  qd=quantile(Diff,c(0.75))
  lb=qd[[1]]
  return(list(ubl=ubl,ubr=ubr,lb=lb))
}


#' function that samples initial values for LoTTA with a discerete prior and binary ourcomes
#' @param x - score data
#' @param t - treatment data
#' @param y - outcome data
#' @param Ct_start - posterior samples of cutoff location (categorized by natural numbers)
#' obtained through "cutoff_initial_DIS.txt"
#' @param cstart - the first point with a positive prior mass
#' @param grid - distance between two consecutive points with nonzero prior mass
#' @param lb - minimum window size (grid size in case of discrete score)
#' @param ubr - maximum value of the window's right boundary point
#' @param ubl - minimum value of the window's left boundary point
#' @param s - seed
#' @param jlb - minimum jump size
#' @return list with initial parameters values for LoTTA model with discrete score
#'  and binary outcomes, and .RNG.seed value
#' @import stats

Initial_DIS_BIN<-function(x,t,y,Ct_start,cstart,grid,lb,ubr,ubl,s,jlb=0.2){
  set.seed(s)
  ct=sample(Ct_start,1)
  c=(ct*grid+cstart-grid)
  tl=mean(t[x<c])
  tr=mean(t[x>=c])
  K=Initial_SHARP_BIN(x,y,c,lb,ubr,ubl,s)
  kl=K[['kl']]
  kr=K[['kr']]
  klt=(kl-ubl)/(c-lb-ubl)
  krt=(kr-c-lb)/(ubr-c-lb)
  kld=kl-c
  krd=kr-c
  a0l=K[['a0l']]
  a1l=K[['a1l']]
  a2l=K[['a2l']]
  a3l=K[['a3l']]

  a0r=K[['a0r']]
  a1r=K[['a1r']]
  a2r=K[['a2r']]
  a3r=K[['a3r']]
  j=ifelse(tr-tl>=abs(a0r-a0l),max(jlb+0.0001,tr-tl),(1+max(jlb,abs(a0r-a0l)))*0.5)

  k1t=runif(1,lb,c-ubl)
  k2t=runif(1,lb,ubr-c)
  a2lt=0
  b2lt=ifelse(j+tl>=0.99,(1-j)*0.8,tl)
  a1lt=0
  b1lt=b2lt
  a1rt=0
  b1rt=b1lt+j
  a2rt=0
  b2rt=b1lt+j

  return(list(ct=ct,j=j,a0l=a0l,a1l=a1l,a2l=a2l,a3l=a3l,a0r=a0r,a1r=a1r,a2r=a2r,a3r=a3r,kl=kl,kr=kr,k1t=k1t,k2t=k2t,a1lt=a1lt,a2lt=a2lt,b2lt=b2lt,a1rt=a1rt,a2rt=a2rt,.RNG.seed=s))
}

#' function that samples initial values for the treatment model with a continuous prior
#' @param x - score data
#' @param t - treatment data
#' @param C_start - posterior samples of cutoff location
#' obtained through "cutoff_initial_CONT.txt"
#' @param clb - left end of an interval on which the prior is supported
#' @param cub - right end of an interval on which the prior is supported
#' @param lb - minimum window size
#' @param ubr - maximum value of the window's right boundary point
#' @param ubl - minimum value of the window's left boundary point
#' @param s - seed
#' @param jlb - minimum jump size
#' @return list with initial parameters values for LoTTA model with continuous score
#' and .RNG.seed value
#' @import stats

Initial_treatment_CONT<-function(x,t,C_start,clb,cub,lb,ubr,ubl,s,jlb=0.2){
  set.seed(s)
  c0=sample(C_start,1)
  c=c0*(cub-clb)+clb
  tl=mean(t[x<c])
  tr=mean(t[x>=c])

  j=ifelse(tr-tl>=jlb+0.0001,tr-tl,(1+jlb)*0.5)

  k1t=runif(1,lb,c-ubl)
  k2t=runif(1,lb,ubr-c)
  a2lt=0
  b2lt=ifelse(j+tl>=0.99,(1-j)*0.8,tl)
  a1lt=0
  b1lt=b2lt
  a1rt=0
  b1rt=b1lt+j
  a2rt=0
  b2rt=b1lt+j

  return(list(c0=c0,j=j,k1t=k1t,k2t=k2t,a1lt=a1lt,a2lt=a2lt,b2lt=b2lt,a1rt=a1rt,a2rt=a2rt,.RNG.seed=s))
}

#' function that samples initial values for the treatment model with a discrete prior
#' @param x - score data
#' @param t - treatment data
#' @param Ct_start - posterior samples of cutoff location (categorized by natural numbers)
#' obtained through "cutoff_initial_dis.txt"c
#' @param cstart - the first point with a positive prior mass
#' @param grid - distance between two consecutive points with nonzero prior mass
#' @param lb - minimum window size (grid size in case of discrete score)
#' @param ubr - maximum value of the window's right boundary point
#' @param ubl - minimum value of the window's left boundary point
#' @param s - seed
#' @param jlb - minimum jump size
#' @return list with initial parameters values for treatment model with discrete score
#'  and .RNG.seed value
#' @import stats

Initial_treatment_DIS<-function(x,t,Ct_start,cstart,grid,lb,ubr,ubl,s,jlb=0.2){
  set.seed(s)
  ct=sample(Ct_start,1)
  c=(ct*grid+cstart-grid)
  tl=mean(t[x<c])
  tr=mean(t[x>=c])

  j=ifelse(tr-tl>=jlb+0.0001,tr-tl,(1+jlb)*0.5)

  k1t=runif(1,lb,c-ubl)
  k2t=runif(1,lb,ubr-c)
  a2lt=0
  b2lt=ifelse(j+tl>=0.99,(1-j)*0.8,tl)
  a1lt=0
  b1lt=b2lt
  a1rt=0
  b1rt=b1lt+j
  a2rt=0
  b2rt=b1lt+j

  return(list(ct=ct,j=j,k1t=k1t,k2t=k2t,a1lt=a1lt,a2lt=a2lt,b2lt=b2lt,a1rt=a1rt,a2rt=a2rt,.RNG.seed=s))
}

#' function that samples initial values for the treatment model with a known cutoff
#' @param x - score data
#' @param t - treatment data
#' @param c - cutoff location
#' @param lb - minimum window size
#' @param ubr - maximum value of the window's right boundary point
#' @param ubl - minimum value of the window's left boundary point
#' @param s - seed
#' @param jlb - minimum jump size
#' @return list with initial parameters values for LoTTA model with continuous score
#' and .RNG.seed value
#' @import stats

Initial_treatment_c<-function(x,t,c,lb,ubr,ubl,s,jlb=0.2){
  set.seed(s)

  tl=mean(t[x<c])
  tr=mean(t[x>=c])

  j=ifelse(tr-tl>=jlb+0.0001,tr-tl,(1+jlb)*0.5)

  k1t=runif(1,lb,c-ubl)
  k2t=runif(1,lb,ubr-c)
  a2lt=0
  b2lt=ifelse(j+tl>=0.99,(1-j)*0.8,tl)
  a1lt=0
  b1lt=b2lt
  a1rt=0
  b1rt=b1lt+j
  a2rt=0
  b2rt=b1lt+j

  return(list(j=j,k1t=k1t,k2t=k2t,a1lt=a1lt,a2lt=a2lt,b2lt=b2lt,a1rt=a1rt,a2rt=a2rt,.RNG.seed=s))
}

#' function that samples initial values for fuzzy LoTTA model with a continuous prior and binary outcomes
#' @param x - score data
#' @param t - treatment data
#' @param y - outcome data
#' @param C_start - posterior samples of cutoff location
#' obtained through "cutoff_initial_CONT.txt"c
#' @param clb - left end of an interval on which the prior is supported
#' @param cub - right end of an interval on which the prior is supported
#' @param lb - minimum window size
#' @param ubr - maximum value of the window's right boundary point
#' @param ubl - minimum value of the window's left boundary point
#' @param s - seed
#' @param jlb - minimum jump size
#' @return list with initial parameters values for LoTTA model with continuous score
#' and binary outcomes, and .RNG.seed value
#' @import stats

Initial_CONT_BIN<-function(x,t,y,C_start,clb,cub,lb,ubr,ubl,s,jlb=0.2){
  set.seed(s)
  MIN=min(x)
  MAx=max(x)
  c0=sample(C_start,1)
  c=c0*(cub-clb)+clb
  tl=mean(t[x<c])
  tr=mean(t[x>=c])
  K=Initial_SHARP_BIN(x,y,c,lb,ubr,ubl,s)
  kl=K[['kl']]
  kr=K[['kr']]
  klt=(kl-ubl)/(c-lb-ubl)
  krt=(kr-c-lb)/(ubr-c-lb)
  kld=kl-c
  krd=kr-c
  a0l=K[['a0l']]
  a1l=K[['a1l']]
  a2l=K[['a2l']]
  a3l=K[['a3l']]

  a0r=K[['a0r']]
  a1r=K[['a1r']]
  a2r=K[['a2r']]
  a3r=K[['a3r']]
  j=ifelse(tr-tl>=abs(a0r-a0l),max(jlb+0.0001,tr-tl),(1+max(jlb,abs(a0r-a0l)))*0.5)

  k1t=runif(1,lb,c-ubl)
  k2t=runif(1,lb,ubr-c)
  a2lt=0
  b2lt=ifelse(j+tl>=0.99,(1-j)*0.8,tl)
  a1lt=0
  b1lt=b2lt
  a1rt=0
  b1rt=b1lt+j
  a2rt=0
  b2rt=b1lt+j

  return(list(c0=c0,j=j,a0l=a0l,a1l=a1l,a2l=a2l,a3l=a3l,a0r=a0r,a1r=a1r,a2r=a2r,a3r=a3r,kl=kl,kr=kr,k1t=k1t,k2t=k2t,a1lt=a1lt,a2lt=a2lt,b2lt=b2lt,a1rt=a1rt,a2rt=a2rt,.RNG.seed=s))
}

#' function that samples initial values for fuzzy LoTTA model with a discrete prior and binary outcomes
#' @param x - score data
#' @param t - treatment data
#' @param y - outcome data
#' @param Ct_start - posterior samples of cutoff location (categorized by natural numbers)
#' obtained through "cutoff_initial_dis.txt"c
#' @param cstart - the first point with a positive prior mass
#' @param grid - distance between two consecutive points with nonzero prior mass
#' @param lb - minimum window size (grid size in case of discrete score)
#' @param ubr - maximum value of the window's right boundary point
#' @param ubl - minimum value of the window's left boundary point
#' @param s - seed
#' @param jlb - minimum jump size
#' @return list with initial parameters values for LoTTA model with discrete score
#'  and continuous outcomes, and .RNG.seed value
#' @import stats

Initial_DIS_CONT<-function(x,t,y,Ct_start,cstart,grid,lb,ubr,ubl,s,jlb=0.2){
  set.seed(s)
  ct=sample(Ct_start,1)
  c=(ct*grid+cstart-grid)
  tl=mean(t[x<c])
  tr=mean(t[x>=c])
  K=Initial_SHARP_CONT(x,y,c,lb,ubr,ubl,s)
  kl=K[['kl']]
  kr=K[['kr']]
  a0l=K[['a0l']]
  a1l=K[['a1l']]
  a2l=K[['a2l']]
  a3l=K[['a3l']]

  a0r=K[['a0r']]
  a1r=K[['a1r']]
  a2r=K[['a2r']]
  a3r=K[['a3r']]
  tau1r=K[['tau1r']]
  tau2pr=K[['tau2pr']]
  tau1l=K[['tau1l']]
  tau2pl=K[['tau2pl']]
  j=ifelse(tr-tl>=jlb+0.0001,tr-tl,(1+jlb)*0.5)

  k1t=runif(1,lb,c-ubl)
  k2t=runif(1,lb,ubr-c)
  a2lt=0
  b2lt=ifelse(j+tl>=0.99,(1-j)*0.8,tl)
  a1lt=0
  b1lt=b2lt
  a1rt=0
  b1rt=b1lt+j
  a2rt=0
  b2rt=b1lt+j
  return(list(ct=ct,j=j,a0l=a0l,a1l=a1l,a2l=a2l,a3l=a3l,a0r=a0r,a1r=a1r,a2r=a2r,a3r=a3r,tau1r=tau1r,tau2pr=tau2pr,tau1l=tau1l,tau2pl=tau2pl,kl=kl,kr=kr,k1t=k1t,k2t=k2t,a1lt=a1lt,a2lt=a2lt,b2lt=b2lt,a1rt=a1rt,a2rt=a2rt,.RNG.seed=s))
}

#' function that samples initial values for fuzzy LoTTA model with a ontinuous prior and continuous outcomes
#' @param x - score data
#' @param t - treatment data
#' @param y - outcome data
#' @param C_start - posterior samples of cutoff location
#' obtained through "cutoff_initial_CONT.txt"
#' @param clb - left end of an interval on which the prior is supported
#' @param cub - right end of an interval on which the prior is supported
#' @param lb - minimum window size
#' @param ubr - maximum value of the window's right boundary point
#' @param ubl - minimum value of the window's left boundary point
#' @param s - seed
#' @param jlb - minimum jump size
#' @return list with initial parameters values for LoTTA model with continuous score
#' and continuous outcomes, and .RNG.seed value
#' @import stats
#'
Initial_CONT_CONT<-function(x,t,y,C_start,clb,cub,lb,ubr,ubl,s,jlb=0.2){
  set.seed(s)
  c0=sample(C_start,1)
  c=c0*(cub-clb)+clb
  tl=mean(t[x<c])
  tr=mean(t[x>=c])
  K=Initial_SHARP_CONT(x,y,c,lb,ubr,ubl,s)
  kl=K[['kl']]
  kr=K[['kr']]
  a0l=K[['a0l']]
  a1l=K[['a1l']]
  a2l=K[['a2l']]
  a3l=K[['a3l']]

  a0r=K[['a0r']]
  a1r=K[['a1r']]
  a2r=K[['a2r']]
  a3r=K[['a3r']]
  tau1r=K[['tau1r']]
  tau2pr=K[['tau2pr']]
  tau1l=K[['tau1l']]
  tau2pl=K[['tau2pl']]
  j=ifelse(tr-tl>=jlb+0.0001,tr-tl,(1+jlb)*0.5)

  k1t=runif(1,lb,c-ubl)
  k2t=runif(1,lb,ubr-c)
  a2lt=0
  b2lt=ifelse(j+tl>=0.99,(1-j)*0.8,tl)
  a1lt=0
  b1lt=b2lt
  a1rt=0
  b1rt=b1lt+j
  a2rt=0
  b2rt=b1lt+j
  return(list(c0=c0,j=j,a0l=a0l,a1l=a1l,a2l=a2l,a3l=a3l,a0r=a0r,a1r=a1r,a2r=a2r,a3r=a3r,tau1r=tau1r,tau2pr=tau2pr,tau1l=tau1l,tau2pl=tau2pl,kl=kl,kr=kr,k1t=k1t,k2t=k2t,a1lt=a1lt,a2lt=a2lt,b2lt=b2lt,a1rt=a1rt,a2rt=a2rt,.RNG.seed=s))
}

#' function that samples initial values for fuzzy LoTTA model with a known cutoff and continuous outcomes
#' @param x - score data
#' @param t - treatment data
#' @param y - outcome data
#' @param c - cutoff value
#' @param lb - minimum window size
#' @param ubr - maximum value of the window's right boundary point
#' @param ubl - minimum value of the window's left boundary point
#' @param s - seed
#' @param jlb - minimum jump size
#' @return list with initial parameters values for LoTTA model with continuous score
#' and continuous outcomes, and .RNG.seed value
#' @import stats
#'

Initial_FUZZy_CONT<-function(x,t,y,c,lb,ubr,ubl,s,jlb=0.2){
  set.seed(s)
  tl=mean(t[x<c])
  tr=mean(t[x>=c])
  K=Initial_SHARP_CONT(x,y,c,lb,ubr,ubl,s)
  kl=K[['kl']]
  kr=K[['kr']]
  a0l=K[['a0l']]
  a1l=K[['a1l']]
  a2l=K[['a2l']]
  a3l=K[['a3l']]

  a0r=K[['a0r']]
  a1r=K[['a1r']]
  a2r=K[['a2r']]
  a3r=K[['a3r']]
  tau1r=K[['tau1r']]
  tau2pr=K[['tau2pr']]
  tau1l=K[['tau1l']]
  tau2pl=K[['tau2pl']]
  j=ifelse(tr-tl>=jlb+0.0001,tr-tl,(1+jlb)*0.5)

  k1t=runif(1,lb,c-ubl)
  k2t=runif(1,lb,ubr-c)
  a2lt=0
  b2lt=ifelse(j+tl>=0.99,(1-j)*0.8,tl)
  a1lt=0
  b1lt=b2lt
  a1rt=0
  b1rt=b1lt+j
  a2rt=0
  b2rt=b1lt+j

  return(list(j=j,a0l=a0l,a1l=a1l,a2l=a2l,a3l=a3l,a0r=a0r,a1r=a1r,a2r=a2r,a3r=a3r,tau1r=tau1r,tau2pr=tau2pr,tau1l=tau1l,tau2pl=tau2pl,kl=kl,kr=kr,k1t=k1t,k2t=k2t,a1lt=a1lt,a2lt=a2lt,b2lt=b2lt,a1rt=a1rt,a2rt=a2rt,.RNG.seed=s))
}

#' function that samples initial values for fuzzy LoTTA model with a known cutoff and binary outcomes
#' @param x - score data
#' @param t - treatment data
#' @param y - outcome data
#' @param c - cutoff location
#' @param lb - minimum window size
#' @param ubr - maximum value of the window's right boundary point
#' @param ubl - minimum value of the window's left boundary point
#' @param s - seed
#' @param jlb - minimum jump size
#' @return list with initial parameters values for fuzzy LoTTA model with a known cutoff and
#' binary outcomes, and .RNG.seed value
#' @import stats
#'

Initial_FUZZy_BIN<-function(x,t,y,c,lb,ubr,ubl,s,jlb=0.2){
  set.seed(s)
  MIN=min(x)
  MAx=max(x)
  tl=mean(t[x<c])
  tr=mean(t[x>=c])
  K=Initial_SHARP_BIN(x,y,c,lb,ubr,ubl,s)
  kl=K[['kl']]
  kr=K[['kr']]
  klt=(kl-ubl)/(c-lb-ubl)
  krt=(kr-c-lb)/(ubr-c-lb)
  kld=kl-c
  krd=kr-c
  a0l=K[['a0l']]
  a1l=K[['a1l']]
  a2l=K[['a2l']]
  a3l=K[['a3l']]

  a0r=K[['a0r']]
  a1r=K[['a1r']]
  a2r=K[['a2r']]
  a3r=K[['a3r']]
  j=ifelse(tr-tl>=abs(a0r-a0l),max(jlb+0.0001,tr-tl),(1+max(jlb,abs(a0r-a0l)))*0.5)

  k1t=runif(1,lb,c-ubl)
  k2t=runif(1,lb,ubr-c)
  a2lt=0
  b2lt=ifelse(j+tl>=0.99,(1-j)*0.8,tl)
  a1lt=0
  b1lt=b2lt
  a1rt=0
  b1rt=b1lt+j
  a2rt=0
  b2rt=b1lt+j

  return(list(j=j,a0l=a0l,a1l=a1l,a2l=a2l,a3l=a3l,a0r=a0r,a1r=a1r,a2r=a2r,a3r=a3r,kl=kl,kr=kr,k1t=k1t,k2t=k2t,a1lt=a1lt,a2lt=a2lt,b2lt=b2lt,a1rt=a1rt,a2rt=a2rt,.RNG.seed=s))
}

#' function that samples initial values for sharp LoTTA model with continuous outcomes
#' @param x - score data
#' @param y - outcome data
#' @param c - cutoff location
#' @param lb - minimum window size
#' @param ubr - maximum value of the window's right boundary point
#' @param ubl - minimum value of the window's left boundary point
#' @param s - seed
#' @return list with initial parameters values for sharp LoTTA model with
#' continuous outcomes, and .RNG.seed value
#' @import stats
#'
Initial_SHARP_CONT<-function(x,y,c,lb,ubr,ubl,s){
  l=as.numeric(length(x))
  set.seed(s)
  I=sample(seq(1,l),l,replace = TRUE)

  xn=x[I]
  yn=y[I]
  (K=optimal_k(xn,yn,c,lb,ubr,ubl))

  kl=K[['kl']]
  kr=K[['kr']]
  kld=kl-c
  krd=kr-c
  a0l=K[['a0l']]
  a1l=K[['a1l']]
  a2l=K[['a2l']]
  a3l=K[['a3l']]

  a0r=K[['a0r']]
  a1r=K[['a1r']]
  a2r=K[['a2r']]
  a3r=K[['a3r']]
  tau1r=1/var(y)
  tau2pr=rbeta(1,1,1)
  tau1l=1/var(y)
  tau2pl=rbeta(1,1,1)

  return(list(a0l=a0l,a1l=a1l,a2l=a2l,a3l=a3l,a0r=a0r,a1r=a1r,a2r=a2r,a3r=a3r,tau1r=tau1r,tau2pr=tau2pr,tau1l=tau1l,tau2pl=tau2pl,kl=kl,kr=kr,.RNG.seed=s))
}

#' function that samples initial values for sharp LoTTA model with binary outcomes
#' @param x - score data
#' @param y - outcome data
#' @param c - cutoff location
#' @param lb - minimum window size
#' @param ubr - maximum value of the window's right boundary point
#' @param ubl - minimum value of the window's left boundary point
#' @param s - seed
#' @return list with initial parameters values for sharp LoTTA model with
#' binary outcomes, and .RNG.seed value
#' @import stats
#'
Initial_SHARP_BIN<-function(x,y,c,lb,ubr,ubl,s){
  l=as.numeric(length(x))
  set.seed(s)
  I=sample(seq(1,l),l,replace = TRUE)

  xn=x[I]
  yn=y[I]
  (K=optimal_k_bin(xn,yn,c,lb,ubr,ubl))

  kl=K[['kl']]
  kr=K[['kr']]
  klt=(kl-ubl)/(c-lb-ubl)
  krt=(kr-c-lb)/(ubr-c-lb)
  kld=kl-c
  krd=kr-c
  a0l=K[['a0l']]
  a1l=K[['a1l']]
  a2l=K[['a2l']]
  a3l=K[['a3l']]

  a0r=K[['a0r']]
  a1r=K[['a1r']]
  a2r=K[['a2r']]
  a3r=K[['a3r']]

  a0l=min(max(a0l,0.01),0.99)
  a0r=min(max(a0r,0.01),0.99)
  a1l=min(max(a1l,(0.98-a0l)/kld),(0.02-a0l)/kld)
  a1r=min(max(a1r,(0.02-a0r)/krd),(0.98-a0r)/krd)
  b0l=logit(a0l)
  b1l=a1l*(invlogit(-b0l)*(1-invlogit(-b0l)))^(-1)
  b0r=logit(a0r)
  b1r=a1r*(invlogit(-b0r)*(1-invlogit(-b0r)))^(-1)

  return(list(a0l=a0l,a1l=a1l,a2l=a2l,a3l=a3l,a0r=a0r,a1r=a1r,a2r=a2r,a3r=a3r,kl=kl,kr=kr,.RNG.seed=s))
}

#' function that checks the type of a prior and whether it is correct
#' @param prior - list with prior parameters or a cutoff value
#' @param minx - minimum value of the score
#' @param maxx - maximum value of the score
#' @param ubl - minimum value of the window's left boundary point
#' @param ubr - maximum value of the window's right boundary point
#' @param lb - minimum window size
#' @return list with type of prior and prior parameters
#'
read_prior<-function(prior,minx,maxx,ubl,ubr,lb){
  if(is.numeric(prior)==TRUE){
    cMIN=ubl+lb
    cMAx=ubr-lb
    if(prior<minx||prior>maxx){
      stop('c outside the range of x')
    }
    if(prior<cMIN || prior>cMAx){
      stop('not enough observations with scores less than c or with scores higher than c. Decrease n_min')
    }
    return(list(type='F',c=prior))
  }
  continous_mandatory=c('clb','cub')
  continous_opt=c('alpha','beta')
  discrete_mandatory=c('cstart','cend')
  discrete_opt=c('weights','grid')
  keys=names(prior)
  check_cont=all(continous_mandatory%in%keys)
  check_dis=all(discrete_mandatory%in%keys)
  check_cont_opt=any(continous_opt%in%keys)
  check_dis_opt=any(discrete_opt%in%keys)
  check_CONT=any(continous_mandatory%in%keys)||check_cont_opt
  check_DIS=any(discrete_mandatory%in%keys)||check_dis_opt
  if(check_CONT==TRUE & check_DIS==TRUE)
  {
    stop('parameters for both discrete and continuous priors were provided.\n For continous prior provide list with clb,cub,alpha,beta.\n For discrete prior provide list with cstart,cend,weights.')
  }
  else{
    if(check_cont==TRUE){
      clb=prior$clb
      cub=prior$cub
      clbMIN=ubl+lb
      cubMAx=ubr-lb
      if(clb<clbMIN || cub>cubMAx){
        stop(paste('not enough observations with scores less than clb or with scores higher than cub. Set clb >',as.character(round(clbMIN,3)),' cub < ',as.character(round(cubMAx,3))))
      }
      alpha=ifelse(is.null(prior$alpha)==TRUE,1,prior$alpha)
      beta=ifelse(is.null(prior$beta)==TRUE,1,prior$beta)
      return(list(type='C',clb=clb,cub=cub,alpha=alpha,beta=beta))
    }
    if(check_dis==TRUE){
      cstart=prior$cstart
      cend=prior$cend
      cstartMIN=ubl+lb
      cendMAx=ubr-lb
      if(cstart<cstartMIN || cend>cendMAx){
        stop(paste('not enough observations with scores less than cstart or with scores higher than cend. Set cstart >',as.character(round(cstartMIN,3)),' cend < ',as.character(round(cendMAx,3))))
      }
      grid=ifelse(is.null(prior$grid)==TRUE,lb,prior$grid)
      if((cend-cstart)%%grid!=0){
        stop('cstart-cend is not multiply of the grid length. ')
      }
      if(is.null(prior$weights)==FALSE){
        weightsl=length(prior$weights)
        if((floor((cend-cstart)/grid)+1)!=weightsl){
          stop(paste('length of weights not compatible with the choice of cstart,cend and grid. Length of weights should equal ',as.character(floor((cend-cstart)/grid)+1)))
        }
      }
      else{
        weights=rep(1,(floor((cend-cstart)/grid)+1))
      }


      return(list(type='D',cstart=cstart,cend=cend,grid=grid,weights=weights))
    }
    if(check_CONT==TRUE){
      if(is.null(prior$clb)==TRUE)
      {
        stop('continous prior: argument "clb" is missing, with no default')
      }
      else{
        stop('continous prior: argument "cub" is missing, with no default')
      }
    }
    if(check_DIS==TRUE){
      if(is.null(prior$cstart)==TRUE){
        stop('discrete prior: argument "cstart" is missing, with no default')
      }
      else{
        stop('discrete prior: argument "cend" is missing, with no default')
      }


    }
  }

}
#' function that searches for initial parameters of outcome function to initiate the sampler
#' @param x - score data
#' @param y - outcome data
#' @param c - cutoff
#' @param ubl - minimum value of the window's left boundary point
#' @param ubr - maximum value of the window's right boundary point
#' @param lb - minimum window size
#' @return - a list with kl and kr and cubic functions parameters
#'
optimal_k<-function(x,y,c,lb,ubr,ubl){

  y1=y[x<c]
  x1=x[x<c]
  y2=y[x>=c]
  x2=x[x>=c]
  A0=c()
  A1=c()
  A2=c()
  A3=c()
  Res=c()
  KR=seq(c+1.0001*lb,(ubr+c+1.0001*lb)*0.5,length.out = 50)
  for(kr in KR){
    d3=data.frame(x=(x2-c),x2=(x2-c)^2*(x2>=kr),x3=(x2-c)^3*(x2>=kr),y=y2)
    (l<-lm(y~.,data=d3))
    Res=append(Res,mean(l$residuals^2))
    A0=append(A0,as.numeric(l$coefficients[1]))
    A1=append(A1,as.numeric(l$coefficients[2]))
    A2=append(A2,ifelse(is.na(as.numeric(l$coefficients[3]))==TRUE,0,as.numeric(l$coefficients[3]) ))
    A3=append(A3,ifelse(is.na(as.numeric(l$coefficients[4]))==TRUE,0,as.numeric(l$coefficients[4]) ))

  }
  kr=KR[which.min(Res)]
  a0r=A0[which.min(Res)]
  a1r=A1[which.min(Res)]
  a2r=A2[which.min(Res)]
  a3r=A3[which.min(Res)]

  A0=c()
  A1=c()
  A2=c()
  A3=c()
  Res=c()
  KL=seq((c-1.0001*lb+ubl)*0.5,c-1.0001*lb,length.out = 50)


  for(kl in KL){
    d3=data.frame(x=(x1-c),x2=(x1-c)^2*(x1<=kl),x3=(x1-c)^3*(x1<=kl),y=y1)
    (l<-lm(y~.,data=d3))
    Res=append(Res,mean(l$residuals^2))
    A0=append(A0,as.numeric(l$coefficients[1]))
    A1=append(A1,as.numeric(l$coefficients[2]))
    A2=append(A2,ifelse(is.na(as.numeric(l$coefficients[3]))==TRUE,0,as.numeric(l$coefficients[3]) ))
    A3=append(A3,ifelse(is.na(as.numeric(l$coefficients[4]))==TRUE,0,as.numeric(l$coefficients[4]) ))

  }
  kl=KL[which.min(Res)]
  a0l=A0[which.min(Res)]
  a1l=A1[which.min(Res)]
  a2l=A2[which.min(Res)]
  a3l=A3[which.min(Res)]


  return(list(kl=kl,kr=kr,a0l=a0l,a1l=a1l,a2l=a2l,a3l=a3l,a0r=a0r,a1r=a1r,a2r=a2r,a3r=a3r))
}

#' function that searches for initial parameters of binary outcome function to initiate the sampler
#' @param x - score data
#' @param y - outcome data
#' @param c - cutoff
#' @param ubl - minimum value of the window's left boundary point
#' @param ubr - maximum value of the window's right boundary point
#' @param lb - minimum window size
#' @return - a list with kl and kr and cubic functions parameters
#'

optimal_k_bin<-function(x,y,c,lb,ubr,ubl){

  y1=y[x<c]
  x1=x[x<c]
  y2=y[x>=c]
  x2=x[x>=c]
  A0=c()
  A1=c()
  A2=c()
  A3=c()
  Res=c()
  KR=seq(c+1.0001*lb,(ubr+c+1.0001*lb)*0.5,length.out = 50)
  for(kr in KR){
    d3=data.frame(x=(x2-c),x2=(x2-c)^2*(x2>=kr),x3=(x2-c)^3*(x2>=kr),y=y2)
    suppressWarnings(l<-glm(y~.,data=d3,family=binomial(link='logit')))
    Res=append(Res,mean(l$residuals^2))
    A0=append(A0,as.numeric(l$coefficients[1]))
    A1=append(A1,as.numeric(l$coefficients[2]))
    A2=append(A2,ifelse(is.na(as.numeric(l$coefficients[3]))==TRUE,0,as.numeric(l$coefficients[3]) ))
    A3=append(A3,ifelse(is.na(as.numeric(l$coefficients[4]))==TRUE,0,as.numeric(l$coefficients[4]) ))
  }
  kr=KR[which.min(Res)]
  krd=kr-c
  a0r=max(min(invlogit(as.numeric(A0[which.min(Res)])),0.95),0.05)
  a1r=max(min(A1[which.min(Res)]*(invlogit(-A0[which.min(Res)])*(1-invlogit(-A0[which.min(Res)]))),(0.99-a0r)/krd),(0.01-a0r)/krd)
  a2r=A2[which.min(Res)]
  a3r=A3[which.min(Res)]

  A0=c()
  A1=c()
  A2=c()
  A3=c()
  Res=c()
  KL=seq((c-1.0001*lb+ubl)*0.5,c-1.0001*lb,length.out = 50)
  for(kl in KL){
    d3=data.frame(x=(x1-c),x2=(x1-c)^2*(x1<=kl),x3=(x1-c)^3*(x1<=kl),y=y1)
    suppressWarnings(l<-glm(y~.,data=d3,family=binomial(link='logit')))
    Res=append(Res,mean(l$residuals^2))
    A0=append(A0,as.numeric(l$coefficients[1]))
    A1=append(A1,as.numeric(l$coefficients[2]))
    A2=append(A2,ifelse(is.na(as.numeric(l$coefficients[3]))==TRUE,0,as.numeric(l$coefficients[3]) ))
    A3=append(A3,ifelse(is.na(as.numeric(l$coefficients[4]))==TRUE,0,as.numeric(l$coefficients[4]) ))
  }
  kl=KL[which.min(Res)]
  kld=kl-c
  a0l=max(min(invlogit(as.numeric(A0[which.min(Res)])),0.95),0.05)
  a1l=min(max(A1[which.min(Res)]*(invlogit(-A0[which.min(Res)])*(1-invlogit(-A0[which.min(Res)]))),(0.99-a0l)/kld),(0.01-a0l)/kld)
  a2l=A2[which.min(Res)]
  a3l=A3[which.min(Res)]

  xc=x-c
  b0l=logit(a0l)
  b1l=a1l*(invlogit(-b0l)*(1-invlogit(-b0l)))^(-1)
  b0r=logit(a0r)
  b1r=a1r*(invlogit(-b0r)*(1-invlogit(-b0r)))^(-1)
  fun=ifelse(x<c,a1l*(xc)+a0l+invlogit(100*(kl-x))*(invlogit((a3l)*(xc)^3+(a2l)*(xc)^2+b1l*(xc)+b0l)-a0l-a1l*xc),a1r*(xc)+a0r+invlogit(100*(x-kr))*(invlogit((a3r)*(xc)^3+(a2r)*(xc)^2+b1r*(xc)+b0r)-a0r-a1r*xc))
  if(max(fun[x<c])>0.999||min(fun[x<c])<0.001){
    a1l=0
    a2l=0
    a3l=0
  }
  if(max(fun[x>=c])>0.999||min(fun[x>=c])<0.001){
    a1r=0
    a2r=0
    a3r=0
  }

  return(list(kl=kl,kr=kr,a0l=a0l,a1l=a1l,a2l=a2l,a3l=a3l,a0r=a0r,a1r=a1r,a2r=a2r,a3r=a3r))
}



