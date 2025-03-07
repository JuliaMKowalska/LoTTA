
#' LoTTA_sharp_CONT
#'
#' Function that fits LoTTA model to the sharp RD data with continuous outcomes. The data does NOT have to be normalized
#' beforehand. We recommend NOT to transform the data before imputing it into the function, except for initial trimming which should be done beforehand.
#' The further trimming for the sensitivity analysis can be done through the function, which ensures that the data is normalized before the trimming.
#'
#' @param X - is the score data
#' @param Y - is the continuous outcome data
#' @param c - specifies the cutoff point
#' @param ci - specifies the probability level 1-alpha for the highest posterior density intervals; default is ci = 0.95
#' @param trimmed - takes as a value NULL or a vector of two values. It specifies potential trimming of the data.
#' If set to NULL no trimming is applied to the data. If a list of two values is provided the data is trimmed
#' to include data points with the score X in between those values; deafult is trimmed=NULL
#' @param outcome_prior - takes as a value a list with elements 'pr' and 'df'. 'pr' specifies precision in the normal priors on the coefficients in the
#' outcome function. 'df' specifies degrees of freedom in the chisquare prior on the precision of the error terms; default is list('pr'=0.0001,'df'=7)
#' @param n_min - specifies the minimum number of data points to which a cubic part of the outcome function is fit to ensure stability of the sampling procedure; default is n_min=25
#' @param param - takes as a value a vector with names of the parameters that are to be sampled; default is the list of all parameters
#' @param normalize - specifies if the data is to be normalized. The data is normalized as follows. X_normalized=(X-d)/s, where d=(min(X)+max(X))*0.5 and s=max(X)-min(X). Y_normalized=(Y-mu)/sd, where mu=mean(Y) and sd=sd(Y).
#' The priors inside the model are specified for the normalized data, in extreme cases not normalizing the data may lead to unreliable results; default is normalize=TRUE
#' @param n.chains - specifies the number of chains in the MCMC sampler; default is n.chains=4
#' @param burnin - specifies the number of burnin iterations without the adaptive iterations; default is burnin=5000
#' @param sample - specifies the number of samples per chain; default is samples=5000
#' @param adapt - specifies the number of adaptive iterations in the MCMC sampler; default is adapt=1000
#' @param inits - initial values for the sampler. By default the initial values are sampled inside the function. If user wants to provide
#' its own values please refer to run.jags manual; default inits=NULL
#' @param method - set to default as 'parallel', which allows to sample the chains in parallel reducing computation time.
#' To read more about possible method values type ?run.jags; default method='parallel'
#' @param seed - specifies the seed for replicability purposes; default is seed=NULL
#' @param ... - other arguments of run.jags function. For more details type ?run.jags
#' @return  The function returns the list with the elements:
#' \itemize{
#'  \item Effect_estimate: contains a list with MAP estimate and HDI of the treatment effect on the original, unnormalized scale;
#'  \item JAGS_output: contains output of the run.jags function for the normalized data if normalize=TRUE, based on this output mixing of the chains can be assessed;
#'  \item Samples: contains posterior samples of the treatment effect (eff);
#'  \item Normalized_data: contains a list of the normalized data (if normalized=TRUE) and the parameters used to normalize the data (see arg normalize);
#'  \item Priors: contains a list of the outcome prior parameters ;
#'  \item Inits contains the list of initial values  and .RNG.seed value
#' }
#' @import stats
#' @import runjags
#' @import bayestestR
#' @export
#' @examples
#' # code to generate the data
#'
#' invlogit <- function(x) {
#'   return(1 / (1 + exp(-x)))
#' }
#'
#' funB <- function(X) {
#'   Y = X
#'   X2 = X[X >= 0]
#'   X1 = X[X < 0]
#'   Y[X < 0] = 1 / (1 + exp(-2 * X1)) - 0.5 + 0.4
#'   Y[X >= 0] = (log(X2 * 2 + 1) - 0.15 * X2^2) * 0.6 - 0.20 + 0.4
#'   return(Y)
#' }
#'
#' funB_sample <- function(X) {
#'   Y = funB(X)+ rnorm(length(X), 0, 0.025)
#'   return(Y)
#' }
#'
#' # data generation
#' N=150 # try different dataset size
#' set.seed(1234)
#' X = sort(runif(N, -1, 1))
#' Y = funB_sample(X)
#' c = 0
#'
#' # running LoTTA function on sharp RDD with continuous outcomes;
#' # cutoff = 0, treatment effect = -0.2
#' # remember to check convergence and adjust burnin, sample and adapt if needed
#' out = LoTTA_sharp_CONT(X, Y, c, burnin = 3000, sample = 2000, adapt = 1000, seed = 1234)
#'

LoTTA_sharp_CONT<-function(X,Y,c,ci=0.95,trimmed=NULL,outcome_prior=list('pr'=0.0001,'df'=7),n_min=25,param=c('eff','a0l','a1l','a2l','a3l','a0r','a1r','a2r','a3r','tau1r','tau2r','tau1l','tau2l','kl','kr')
                           ,normalize=TRUE,n.chains=4,burnin=10000,sample=5000,adapt=1000,inits=NULL,method='parallel',seed=NULL,...){
  normalized_parameters=list(X=X,Y=Y,c_normalized=c,d_X=0,s_X=1,mu_Y=0,sd_Y=1)
  if(all(unique(Y)%in%c(0,1))==TRUE){
    warning('Binary outcome detected. Using LoTTA_sharp_BIN is recommended. ')
  }
  if(is.null(trimmed)==TRUE){
    b=bounds(X,n_min)
    ubr=b$ubr
    ubl=b$ubl
    lb=b$lb
    read_prior(c,min(X),max(X),ubl,ubr,lb)
  }
  else{
    if(is.vector(trimmed)==TRUE){
      if(length(unique(trimmed))==2){
        Xtr=X[X>=min(trimmed)&X<=max(trimmed)]
        if(length(Xtr)==0){
          stop('there is no data in the trimmed dataset. Change the values of trimmed')
        }
        else{
          if(length(Xtr)<100){
            warning('there is less than 100 datapoints in the trimmed dataset')
          }
          b=bounds(Xtr,n_min)
          ubr=b$ubr
          ubl=b$ubl
          lb=b$lb
          read_prior(c,min(Xtr),max(Xtr),ubl,ubr,lb)
        }
      }
      else{
        stop('provided trimmed variable is neither NULL nor a vector with precisely two unique values')
      }
    }
    else{
      stop('wrong variable type. trimmed must be either NULL or a vector with precisely two unique values')
    }
  }
  if(normalize==TRUE || is.null(trimmed)==FALSE){
    x_list=normalize_cont_X(X,trimmed,normalize)
    y_list=normalize_cont_Y(Y,X,trimmed,normalize)
    X=x_list$X
    Y=y_list$Y
    c=(c-x_list$d)/x_list$s
    normalized_parameters=list(X=X,Y=Y,c=c,d_X=x_list$d,s_X=x_list$s,mu_Y=y_list$mu,sd_Y=y_list$sd)
  }

  seed=ifelse(is.null(seed)==TRUE,sample(1:1000000,1),seed)
  b=bounds(X,n_min)
  ubr=b$ubr
  ubl=b$ubl
  lb=b$lb
  if(n_min<25){
    warning('setting n_min<25 may cause convergence problems in MCMC and lead to unstable results')
  }
  pr=outcome_prior$pr
  df=outcome_prior$df
  dat=list(N=length(X),x=X,y=Y,c=c,ubr=ubr,ubl=ubl,lb=lb,pr=pr,df=df)
  if(is.null(inits)==TRUE){
    inits=list()
    for (i in 1:n.chains) {
      (inits[[i]]=Initial_SHARP_CONT(X,Y,c,lb,ubr,ubl,seed+i))
    }
  }
  model_path <- system.file("models", 'LoTTA_SHARP_CONT.txt', package = "LoTTA")
  posterior<- tryCatch({run.jags(model_path, data=dat,inits = inits,monitor=param,burnin = burnin,sample=sample,adapt = adapt,n.chains = n.chains,method = method,...)},error = function(e) {
    # Handle the error and return a specific value
    message("An error occurred, only initial values returned: ", e$message)
    return(NA)
    # Return a specific value on error
  })
  if (identical(posterior, NA)) {
    return(inits)  # Exit the entire function
  }

  Samples=combine.mcmc(posterior)
  Eff=as.numeric(Samples[,'eff'])*normalized_parameters$sd_Y
  mapeff=as.numeric(map_estimate(Eff))
  hdieff=as.numeric(ci(Eff,ci=ci,method='HDI'))[2:3]
  Res=list(Effect_estimate=list(MAP=mapeff,CI=hdieff),JAGS_output=posterior,Samples=data.frame('eff'=Eff),Normalized_data=normalized_parameters,Priors=list('Outcome prior'=outcome_prior),Inits=inits)
  return(Res)
}


#' LoTTA_sharp_BIN
#'
#' Function that fits LoTTA model to the sharp RD data with binary outcomes. The score does NOT have to be normalized
#' beforehand. We recommend NOT to transform the data before imputing it into the function, except for initial trimming of the score which should be done beforehand.
#' The further trimming for the sensitivity analysis can be done through the function, which ensures that the data is normalized before the trimming.
#'
#' @param X - is the score data
#' @param Y - is the binary outcome data
#' @param c - specifies the cutoff point
#' @param ci - specifies the probability level 1-alpha for the highest posterior density intervals; default is ci = 0.95
#' @param trimmed - takes as a value NULL or a vector of two values. It specifies potential trimming of the data.
#' If set to NULL no trimming is applied to the data. If a list of two values is provided the data is trimmed
#' to include data points with the score X in between those values; deafult is trimmed=NULL
#' @param outcome_prior - takes as a value a list with elements 'pr'. 'pr' specifies precision in the normal priors on the coefficients in the
#' outcome function; default is list('pr'=0.0001)
#' @param n_min - specifies the minimum number of data points to which a cubic part of the outcome function is fit to ensure stability of the sampling procedure; default is n_min=25
#' @param param - takes as a value a vector with names of the parameters that are to be sampled; default is the list of all parameters
#' @param normalize - specifies if the data is to be normalized. The data is normalized as follows. X_normalized=(X-d)/s, where d=(min(X)+max(X))*0.5 and s=max(X)-min(X).
#' The priors inside the model are specified for the normalized data, in extreme cases not normalizing the data may lead to unreliable results; default is normalize=TRUE
#' @param n.chains - specifies the number of chains in the MCMC sampler; default is n.chains=4
#' @param burnin - specifies the number of burnin iterations without the adaptive iterations; default is burnin=5000
#' @param sample - specifies the number of samples per chain; default is samples=5000
#' @param adapt - specifies the number of adaptive iterations in the MCMC sampler; default is adapt=1000
#' @param inits - initial values for the sampler. By default the initial values are sampled inside the function. If user wants to provide
#' its own values please refer to run.jags manual; default inits=NULL
#' @param method - set to default as 'parallel', which allows to sample the chains in parallel reducing computation time.
#' To read more about possible method values type ?run.jags; default method='parallel'
#' @param seed - specifies the seed for replicability purposes; default is seed=NULL
#' @param ... - other arguments of run.jags function. For more details type ?run.jags
#' @return  The function returns the list with the elements:
#' \itemize{
#'  \item Effect_estimate: contains a list with MAP estimate and HDI of the treatment effect on the original, unnormalized scale;
#'  \item JAGS_output: contains output of the run.jags function for the normalized data if normalize=TRUE, based on this output mixing of the chains can be assessed;
#'  \item Samples: contains posterior samples of the treatment effect (eff);
#'  \item Normalized_data: contains a list of the normalized data (if normalized=TRUE) and the parameters used to normalize the data (see arg normalize);
#'  \item Priors: contains a list of the outcome prior parameters ;
#'  \item Inits contains the list of initial values  and .RNG.seed value
#' }
#' @import stats
#' @import runjags
#' @import bayestestR
#' @export
#' @examples
#' # code to generate the data
#'
#' invlogit <- function(x) {
#'   return(1 / (1 + exp(-x)))
#' }
#'
#' fun_prob55 <- function(X) {
#'   P = rep(0, length(X))
#'   P[X >= 0.] = invlogit((8.5 * X[X >= 0.] - 1.5)) / 10.5 + 0.65 - 0.0007072
#'   P[X < 0.] = (X[X < 0.] + 1)^4 / 15 + 0.05
#'   return(P)
#' }
#'
#' sample_prob55 <- function(X) {
#'   P = rep(0, length(X))
#'   P[X >= 0.] = invlogit((8.5 * X[X >= 0.] - 1.5)) / 10.5 + 0.65 - 0.0007072
#'   P[X < 0.] = (X[X < 0.] + 1)^4 / 15 + 0.05
#'   T = rep(0, length(X))
#'   for (j in 1:length(X)) {
#'     T[j] = sample(c(1, 0), 1, prob = c(P[j], 1 - P[j]))
#'   }
#'   return(T)
#' }
#'
#' # data generation
#' N=150 # try different dataset sizes
#' set.seed(1234)
#' X = sort(runif(N, -1, 1))
#' Y = sample_prob55(X)
#' c = 0
#'
#' # running LoTTA function on sharp RDD with a binary outcome
#' # cutoff = 0, treatment effect = 0.55
#' # remember to check convergence and adjust burnin, sample and adapt if needed
#' out = LoTTA_sharp_BIN(X, Y, c, burnin = 3000, sample = 2000, adapt = 1000, seed = 1234)
#'

LoTTA_sharp_BIN<-function(X,Y,c,ci=0.95,trimmed=NULL,outcome_prior=list('pr'=0.0001),n_min=25,param=c('eff','a0l','a1l','a2l','a3l','a0r','a1r','a2r','a3r','kl','kr')
                          ,normalize=TRUE,n.chains=4,burnin=5000,sample=5000,adapt=1000,inits=NULL,method='parallel',seed=NULL,...){
  if(all(unique(Y)%in%c(0,1))==FALSE){
    stop('Y is not a binary variable')
  }
  if(is.null(trimmed)==TRUE){
    b=bounds(X,n_min)
    ubr=b$ubr
    ubl=b$ubl
    lb=b$lb
    read_prior(c,min(X),max(X),ubl,ubr,lb)
  }
  else{
    if(is.vector(trimmed)==TRUE){
      if(length(unique(trimmed))==2){
        Xtr=X[X>=min(trimmed)&X<=max(trimmed)]
        if(length(Xtr)==0){
          stop('there is no data in the trimmed dataset. Change the values of trimmed')
        }
        else{
          if(length(Xtr)<100){
            warning('there is less than 100 datapoints in the trimmed dataset')
          }
          b=bounds(Xtr,n_min)
          ubr=b$ubr
          ubl=b$ubl
          lb=b$lb
          read_prior(c,min(Xtr),max(Xtr),ubl,ubr,lb)
        }
      }
      else{
        stop('provided trimmed variable is neither NULL nor a vector with precisely two unique values')
      }
    }
    else{
      stop('wrong variable type. trimmed must be either NULL or a vector with precisely two unique values')
    }
  }
  normalized_parameters=list(X=X,Y=Y,c_normalized=c,d_X=0,s_X=1,mu_Y=0,sd_Y=1)
  if(normalize==TRUE || is.null(trimmed)==FALSE){
    x_list=normalize_cont_X(X,trimmed,normalize)
    y_list=normalize_cont_Y(Y,X,trimmed,normalize=FALSE)
    X=x_list$X
    Y=y_list$Y
    c=(c-x_list$d)/x_list$s
    normalized_parameters=list(X=X,Y=Y,c=c,d_X=x_list$d,s_X=x_list$s,mu_Y=y_list$mu,sd_Y=y_list$sd)
  }
  seed=ifelse(is.null(seed)==TRUE,sample(1:1000000,1),seed)
  b=bounds(X,n_min)
  ubr=b$ubr
  ubl=b$ubl
  lb=b$lb

  if(n_min<25){
    warning('setting n_min<25 may cause convergence problems in MCMC and lead to unstable results')
  }
  pr=outcome_prior$pr
  dat=list(N=length(X),x=X,y=Y,c=c,ubr=ubr,ubl=ubl,lb=lb,pr=pr)
  if(is.null(inits)==TRUE){
    inits=list()
    for (i in 1:n.chains) {
      (inits[[i]]=Initial_SHARP_BIN(X,Y,c,lb,ubr,ubl,(seed+i)))

    }
  }
  model_path <- system.file("models", 'LoTTA_SHARP_BIN.txt', package = "LoTTA")
  posterior<- tryCatch({run.jags(model_path, data=dat,inits = inits,monitor=param,burnin = burnin,sample=sample,adapt = adapt,n.chains = n.chains,method = method,...)},error = function(e) {
    # Handle the error and return a specific value
    message("An error occurred, only initial values returned: ", e$message)
    return(NA)
    # Return a specific value on error
  })
  if (identical(posterior, NA)) {
    return(inits)  # Exit the entire function
  }
  Samples=combine.mcmc(posterior)
  Eff=as.numeric(Samples[,'eff'])*normalized_parameters$sd_Y
  mapeff=as.numeric(map_estimate(Eff))
  hdieff=as.numeric(ci(Eff,ci=ci,method='HDI'))[2:3]
  Res=list(Effect_estimate=list(MAP=mapeff,CI=hdieff),JAGS_output=posterior,Samples=data.frame('eff'=Eff),Normalized_data=normalized_parameters,Priors=list('Outcome prior'=outcome_prior),Inits=inits)
  return(Res)
}

#' LoTTA_fuzzy_BIN
#'
#' Function that fits LoTTA model to the fuzzy RD data with binary outcomes with an either known or unknown/suspected cutoff.
#' It supports two types of priors on the cutoff location: a scaled beta distribution of the form beta(alpha,beta)(cub-clb)+clb and a discrete distribution with the support of the form cstart+grid i for i=0,...,(cend-cstart)/grid.
#' The score does NOT have to be normalized beforehand. We recommend NOT to transform the data before imputing it into the function, except for initial trimming of the score which should be done beforehand.
#' The further trimming for the sensitivity analysis can be done through the function, which ensures that the data is normalized before the trimming.
#'
#' @param X - is the score data
#' @param T - is the treatment allocation data
#' @param Y - is the binary outcome data
#' @param c_prior - specifies the cutoff prior in case of the unknown cutoff or the cutoff point
#' if the cutoff is known. Takes as value a number if the cutoff is known or a list of values otherwise.
#' For a continuous prior the list requires the following elements: clb - left end of the interval cub - right end of the interval
#' in which the scaled and translated beta distribution is defined, alpha (optional) - shape parameter, default value = 1, beta (optional) - shape parameter, default value = 1.
#' For a discrete prior the list requires the following elements: cstart - first point with positive prior mass, cend - last point with positive prior mass, grid - distance between the consecutive points in the support
#' weights (optional) - vector of weights assigned to each point in the support, default is vector of 1's (uniform distribution)
#' @param jlb - minimum jump size
#' @param ci - specifies the probability level 1-alpha for the highest posterior density intervals; default is ci = 0.95
#' @param trimmed - takes as a value NULL or a vector of two values. It specifies potential trimming of the data.
#' If set to NULL no trimming is applied to the data. If a list of two values is provided the data is trimmed
#' to include data points with the score X in between those values; default is trimmed=NULL
#' @param outcome_prior - takes as a value a list with elements 'pr'. 'pr' specifies precision in the normal priors on the coefficients in the
#' outcome function; default is list('pr'=0.0001)
#' @param n_min - specifies the minimum number of data points to which a cubic part of the outcome function is fit to ensure stability of the sampling procedure; default is n_min=25
#' @param param - takes as a value a vector with names of the parameters that are to be sampled; default is the list of all parameters
#' @param normalize - specifies if the data is to be normalized. The data is normalized as follows. If the prior is continuous: X_normalized=(X-d)/s, where d=(min(X)+max(X))*0.5 and s=max(X)-min(X).
#' If the prior is discrete: X_normalized=X/s, where s=10^m, where m is chosen so that |max(abs(X))-1| is minimal.
#' The priors inside the model are specified for the normalized data, in extreme cases not normalizing the data may lead to unreliable results; default is normalize=TRUE
#' @param n.chains - specifies the number of chains in the MCMC sampler; default is n.chains=4
#' @param burnin - specifies the number of burnin iterations without the adaptive iterations; default is burnin=5000
#' @param sample - specifies the number of samples per chain; default is samples=5000
#' @param adapt - specifies the number of adaptive iterations in the MCMC sampler; default is adapt=1000
#' @param inits - initial values for the sampler. By default the initial values are sampled inside the function. If user wants to provide
#' its own values please refer to run.jags manual; default inits=NULL
#' @param method - set to default as 'parallel', which allows to sample the chains in parallel reducing computation time.
#' To read more about possible method values type ?run.jags; default method='parallel'
#' @param seed - specifies the seed for replicability purposes; default is seed=NULL
#' @param ... - other arguments of run.jags function. For more details type ?run.jags
#' @return  The function returns the list with the elements:
#' \itemize{
#'  \item Effect_estimate: contains a list with MAP estimate and HDI of the treatment effect, the cutoff location (if unknown) and the discontinuity size in the treatment probability function (compliance rate at c) on the original, unnormalized scale;
#'  \item JAGS_output: contains output of the run.jags function for the normalized data if normalize=TRUE, based on this output mixing of the chains can be assessed;
#'  \item Samples: contains posterior samples of the treatment effect (eff), cutoff location (c) if unknown, and compliance rate (j);
#'  \item Normalized_data: contains a list of the normalized data (if normalized=TRUE) and the parameters used to normalize the data (see arg normalize);
#'  \item Priors: contains a list of the priors' parameters ;
#'  \item Inits contains the list of initial values  and .RNG.seed value
#' }
#' @import stats
#' @import runjags
#' @import bayestestR
#' @export
#' @examples
#' # code to generate the data
#'
#' invlogit <- function(x) {
#'   return(1 / (1 + exp(-x)))
#' }
#'
#' fun_prob55 <- function(X) {
#'   P = rep(0, length(X))
#'   P[X >= 0.] = invlogit((8.5 * X[X >= 0.] - 1.5)) / 10.5 + 0.65 - 0.0007072
#'   P[X < 0.] = (X[X < 0.] + 1)^4 / 15 + 0.05
#'   return(P)
#' }
#'
#' sample_prob55 <- function(X) {
#'   P = rep(0, length(X))
#'   P[X >= 0.] = invlogit((8.5 * X[X >= 0.] - 1.5)) / 10.5 + 0.65 - 0.0007072
#'   P[X < 0.] = (X[X < 0.] + 1)^4 / 15 + 0.05
#'   T = rep(0, length(X))
#'   for (j in 1:length(X)) {
#'     T[j] = sample(c(1, 0), 1, prob = c(P[j], 1 - P[j]))
#'   }
#'   return(T)
#' }
#'
#' funB <- function(X) {
#'   Y = X
#'   X2 = X[X >= 0]
#'   X1 = X[X < 0]
#'   Y[X < 0] = 1 / (1 + exp(-2 * X1)) - 0.5 + 0.4
#'   Y[X >= 0] = (log(X2 * 2 + 1) - 0.15 * X2^2) * 0.6 - 0.20 + 0.4
#'   return(Y)
#' }
#'
#' funB_sample_binary <- function(X) {
#'   Y = X
#'   P = funB(X)
#'   for (j in 1:length(X)) {
#'     Y[j] = sample(c(1, 0), 1, prob = c(P[j], 1 - P[j]))
#'   }
#'   return(Y)
#' }
#'
#' N=150 # try different dataset size
#' set.seed(1234)
#' X = sort(runif(N, -1, 1))
#' T = sample_prob55(X)
#' Y = funB_sample_binary(X)
#'
#' # comment out to try different priors:
#' c_prior=0 # known cutoff c=0
#' # c_prior=list(clb=-0.25,cub=0.25) # uniform prior on the interval [-0.25,0.25]
#' # c_prior=list(cstart=-0.25,cend=0.25,grid=0.05) # uniform discrete prior on -0.25, -0.2, ..., 0.25
#'
#' # running LoTTA model on fuzzy RDD with binary outcomes;
#' # cutoff = 0, compliance rate = 0.55, treatment effect = -0.3636364
#' # remember to check convergence and adjust burnin, sample and adapt if needed
#' out = LoTTA_fuzzy_BIN(X,T,Y,c_prior,burnin = 2000,sample = 1000,adapt=1000,seed=1234)

LoTTA_fuzzy_BIN<-function(X,T,Y,c_prior,jlb=0.2,ci=0.95,trimmed=NULL,outcome_prior=list('pr'=0.0001),n_min=25,param=c('c','j','kl','kr','eff','a0l','a1l','a2l','a3l','a0r','a1r','a2r','a3r','b1lt','a1lt','a2lt','b2lt','b1rt','a1rt','a2rt','b2rt','k1t','k2t')
                          ,normalize=TRUE,n.chains=4,burnin=10000,sample=1000,adapt=500,inits=NULL,method='parallel',seed=NULL,...){

  if(all(unique(Y)%in%c(0,1))==FALSE){
    stop('Y is not a binary variable')
  }
  if(is.null(trimmed)==TRUE){
    b=bounds(X,n_min)
    ubr=b$ubr
    ubl=b$ubl
    lb=b$lb
    prior=read_prior(c_prior,min(X),max(X),ubl,ubr,lb)
  }
  else{
    if(is.vector(trimmed)==TRUE){
      if(length(unique(trimmed))==2){
        Xtr=X[X>=min(trimmed)&X<=max(trimmed)]
        if(length(Xtr)==0){
          stop('there is no data in the trimmed dataset. Change the values of trimmed')
        }
        else{
          if(length(Xtr)<100){
            warning('there is less than 100 datapoints in the trimmed dataset')
          }
          b=bounds(Xtr,n_min)
          ubr=b$ubr
          ubl=b$ubl
          lb=b$lb
          prior=read_prior(c_prior,min(Xtr),max(Xtr),ubl,ubr,lb)
        }
      }
      else{
        stop('provided trimmed variable is neither NULL nor a vector with precisely two unique values')
      }
    }
    else{
      stop('wrong variable type. trimmed must be either NULL or a vector with precisely two unique values')
    }
  }
  normalized_parameters=list(X=X,Y=Y,T=T,d_X=0,s_X=1)
  if(normalize==TRUE || is.null(trimmed)==FALSE){
    if(prior$type=='D'){
      x_list=normalize_dis_X(X,trimmed,normalize)
    }
    else{
      x_list=normalize_cont_X(X,trimmed,normalize)
    }

    y_list=normalize_cont_Y(Y,X,trimmed,normalize=FALSE)
    t_list=normalize_cont_Y(T,X,trimmed,normalize=FALSE)
    X=x_list$X
    Y=y_list$Y
    T=t_list$Y
    normalized_parameters=list(X=X,Y=Y,T=T,d_X=x_list$d,s_X=x_list$s)

  }

  pr=outcome_prior$pr
  seed=ifelse(is.null(seed)==TRUE,sample(1:1000000,1),seed)


  if(prior$type=='C'){

    alpha=prior$alpha
    beta=prior$beta
    clb=prior$clb
    cub=prior$cub
    normalized_parameters$clb=clb
    normalized_parameters$cub=cub
    if(normalize==TRUE){
      clb=(clb-normalized_parameters$d_X)/normalized_parameters$s_X
      cub=(cub-normalized_parameters$d_X)/normalized_parameters$s_X
      normalized_parameters$clb=clb
      normalized_parameters$cub=cub
    }
    b=bounds(X,n_min)
    ubr=b$ubr
    ubl=b$ubl
    ubrt=b$ubr
    ublt=b$ubl
    lb=b$lb

    datT=list(N=length(X),x=X,t=T,jlb=jlb,clb=clb,cub=cub,alpha=alpha,beta=beta)
    param_c=c('c0')
    c0=ifelse(alpha>1 & beta>1,(alpha-1)/(alpha+beta-2),0.5)
    initc1=list(c0=c0,al=0.5-0.5*((1+jlb)*0.5),j=(1+jlb)*0.5,.RNG.seed=seed,.RNG.name="base::Mersenne-Twister")
    model_path1 <- system.file("models", 'cutoff_initial_CONT.txt', package = "LoTTA")
    suppressWarnings(dat1_c<- run.jags(model_path1,inits = list(initc1) ,data=datT,monitor=param_c,burnin = 1000,sample=500,adapt = 100,n.chains = 1,method = 'simple'))
    C_start=as.numeric(combine.mcmc(dat1_c$mcmc))
    if(is.null(inits)==TRUE){
      inits=list()
      for (i in 1:n.chains) {

        suppressWarnings(inits[[i]]<-Initial_CONT_BIN(X,T,Y,C_start,clb,cub,lb,ubr,ubl,seed+i,jlb))
      }
    }
    dat=list(N=length(X),x=X,t=T,y=Y,ubr=ubr,ubl=ubl,ubrt=ubrt,ublt=ublt,lb=lb,jlb=jlb,clb=clb,cub=cub,alpha=alpha,beta=beta,pr=pr,seed=seed)
    model_path2 <- system.file("models", 'LoTTA_CONT_BIN.txt', package = "LoTTA")
    posterior<- tryCatch({run.jags(model_path2, data=dat,inits = inits,monitor=param,burnin = burnin,sample=sample,adapt = adapt,n.chains = n.chains,method = method,...)},error = function(e) {
      # Handle the error and return a specific value
      message("An error occurred, only initial values are returned: ", e$message)
      return(NA)
      # Return a specific value on error
    })
    if (identical(posterior, NA)) {
      return(inits)  # Exit the entire function
    }
    Samples=combine.mcmc(posterior)
    C=as.numeric(Samples[,'c'])*normalized_parameters$s_X+normalized_parameters$d_X
    mapc=as.numeric(map_estimate(C))
    hdic=as.numeric(ci(C,ci=ci,method='HDI'))[2:3]
  }

  if(prior$type=='D'){

    cstart=prior$cstart
    cend=prior$cend
    grid=prior$grid
    weights=prior$weights
    normalized_parameters$cstart=cstart
    normalized_parameters$cend=cend
    normalized_parameters$grid=grid

    if(normalize==TRUE){
      cstart=(cstart-normalized_parameters$d_X)/normalized_parameters$s_X
      cend=(cend-normalized_parameters$d_X)/normalized_parameters$s_X
      grid=grid/normalized_parameters$s_X
      normalized_parameters$cstart=cstart
      normalized_parameters$cend=cend
      normalized_parameters$grid=grid
    }
    b=bounds(X,n_min)
    ubr=b$ubr
    ubl=b$ubl
    ubrt=b$ubr
    ublt=b$ubl
    lb=b$lb

    datT=list(N=length(X),x=X,t=T,jlb=jlb,cstart=cstart,grid=grid,weights=weights,lb=lb)
    param_c=c('ct')
    ct=ifelse(length(unique(weights))==1,ceiling(length(weights)*0.5),which.max(weights))
    initc1=list(ct=ct,al=0.5-0.5*((1+jlb)*0.5),j=(1+jlb)*0.5,.RNG.seed=seed,.RNG.name="base::Mersenne-Twister")
    model_path1 <- system.file("models", 'cutoff_initial_DIS.txt', package = "LoTTA")
    suppressWarnings(dat1_c<- run.jags(model_path1,inits = list(initc1) ,data=datT,monitor=param_c,burnin = 1000,sample=500,adapt = 100,n.chains = 1,method = 'simple'))
    C_start=as.numeric(combine.mcmc(dat1_c$mcmc))
    if(is.null(inits)==TRUE){
      inits=list()
      for (i in 1:n.chains) {

        suppressWarnings(inits[[i]]<-Initial_DIS_BIN(X,T,Y,C_start,cstart,grid,lb,ubr,ubl,seed+i,jlb))
      }
    }
    dat=list(N=length(X),x=X,t=T,y=Y,ubr=ubr,ubl=ubl,ubrt=ubrt,ublt=ublt,lb=lb,jlb=jlb,cstart=cstart,grid=grid,weights=weights,pr=pr,seed=seed)
    model_path2 <- system.file("models", 'LoTTA_DIS_BIN.txt', package = "LoTTA")
    posterior<- tryCatch({run.jags(model_path2, data=dat,inits = inits,monitor=param,burnin = burnin,sample=sample,adapt = adapt,n.chains = n.chains,method = method,...)},error = function(e) {
      # Handle the error and return a specific value
      message("An error occurred, only initial values are returned: ", e$message)
      return(NA)
      # Return a specific value on error
    })
    if (identical(posterior, NA)) {
      return(inits)  # Exit the entire function
    }
    Samples=combine.mcmc(posterior)
    C=as.numeric(Samples[,'c'])*normalized_parameters$s_X+normalized_parameters$d_X
    mapc=which.max(table(C))
    mapc=as.numeric(names(mapc))
    hdic=as.numeric(ci(C,ci=ci,method='HDI'))[2:3]
  }


  if(prior$type=='F'){
    c=c_prior
    normalized_parameters$c=c
    if(normalize==TRUE || is.null(trimmed)==FALSE){
      c=(c-x_list$d)/x_list$s
      normalized_parameters$c=c
    }
    b=bounds(X,25)
    ubr=b$ubr
    ubl=b$ubl
    ubrt=b$ubr
    ublt=b$ubl
    lb=b$lb
    param=param[!(param %in% c('c'))]
    if(is.null(inits)==TRUE){
      inits=list()
      for (i in 1:n.chains) {

        inits[[i]]=Initial_FUZZY_BIN(X,T,Y,c,lb,ubr,ubl,seed+i,jlb)
      }
    }
    dat=list(N=length(X),x=X,t=T,y=Y,c=c,ubr=ubr,ubl=ubl,ubrt=ubrt,ublt=ublt,lb=lb,jlb=jlb,pr=pr,seed=seed)
    model_path <- system.file("models", 'LoTTA_BIN_c.txt', package = "LoTTA")
    posterior<- tryCatch({run.jags(model_path, data=dat,inits = inits,monitor=param,burnin = burnin,sample=sample,adapt = adapt,n.chains = n.chains,method = method,...)},error = function(e) {
      # Handle the error and return a specific value
      message("An error occurred, only initial values are returned: ", e$message)
      return(NA)
      # Return a specific value on error
    })
    if (identical(posterior, NA)) {
      return(inits)  # Exit the entire function
    }
    Samples=combine.mcmc(posterior)
  }

  Eff=as.numeric(Samples[,'eff'])
  mapeff=as.numeric(map_estimate(Eff))
  hdieff=as.numeric(ci(Eff,ci=ci,method='HDI'))[2:3]
  J=as.numeric(Samples[,'j'])
  mapj=as.numeric(map_estimate(J))
  hdij=as.numeric(ci(J,ci=ci,method='HDI'))[2:3]


  if(is.numeric(c_prior)==FALSE){
    Samples=data.frame('eff'=Eff,'c'=C,'j'=J)
    Res=list(Effect_estimate=list('MAP treatment effect'=mapeff,'CI treatment effect'=hdieff,'MAP cutoff'=mapc,'CI cutoff'=hdic,'MAP compliance rate'=mapj,'CI compliance rate'=hdij),JAGS_output=posterior,Samples=Samples,Normalized_data=normalized_parameters,Priors=list('Cutoff prior'=prior[-1],'Outcome prior'=list(pr=pr,df=df)),Inits=inits)
  }
  else{
    Samples=data.frame('eff'=Eff,'j'=J)
    Res=list(Effect_estimate=list('MAP treatment effect'=mapeff,'CI treatment effect'=hdieff,'MAP compliance rate'=mapj,'CI compliance rate'=hdij),JAGS_output=posterior,Samples=Samples,Normalized_data=normalized_parameters,Priors=list('Cutoff'=c_prior,'Outcome prior'=list(pr=pr,df=df)),Inits=inits)
  }

  return(Res)


}

#' LoTTA_fuzzy_CONT
#'
#' Function that fits LoTTA model to the fuzzy RD data with continuous outcomes with an either known or unknown/suspected cutoff.
#' It supports two types of priors on the cutoff location: a scaled beta distribution of the form beta(alpha,beta)(cub-clb)+clb and a discrete distribution with the support of the form cstart+grid i for i=0,...,(cend-cstart)/grid.
#' The score does NOT have to be normalized beforehand. We recommend NOT to transform the data before imputing it into the function, except for initial trimming of the score which should be done beforehand.
#' The further trimming for the sensitivity analysis can be done through the function, which ensures that the data is normalized before the trimming.
#'
#' @param X - is the score data
#' @param T - is the treatment allocation data
#' @param Y - is the binary outcome data
#' @param c_prior - specifies the cutoff prior in case of the unknown cutoff or the cutoff point
#' if the cutoff is known. Takes as value a number if the cutoff is known or a list of values otherwise.
#' For a continuous prior the list requires the following elements: clb - left end of the interval cub - right end of the interval
#' in which the scaled and translated beta distribution is defined, alpha (optional) - shape parameter, default value = 1, beta (optional) - shape parameter, default value = 1.
#' For a discrete prior the list requires the following elements: cstart - first point with positive prior mass, cend - last point with positive prior mass, grid - distance between the consecutive points in the support
#' weights (optional) - vector of weights assigned to each point in the support, default is vector of 1's (uniform distribution)
#' @param jlb - minimum jump size
#' @param ci - specifies the probability level 1-alpha for the highest posterior density intervals; default is ci = 0.95
#' @param trimmed - takes as a value NULL or a vector of two values. It specifies potential trimming of the data.
#' If set to NULL no trimming is applied to the data. If a list of two values is provided the data is trimmed
#' to include data points with the score X in between those values; default is trimmed=NULL
#' @param outcome_prior - takes as a value a list with elements 'pr' and 'df'. 'pr' specifies precision in the normal priors on the coefficients in the
#' outcome function. 'df' specifies degrees of freedom in the chisquare prior on the precision of the error terms; default is list('pr'=0.0001,'df'=7)
#' @param n_min - specifies the minimum number of data points to which a cubic part of the outcome function is fit to ensure stability of the sampling procedure; default is n_min=25
#' @param param - takes as a value a vector with names of the parameters that are to be sampled; default is the list of all parameters
#' @param normalize - specifies if the data is to be normalized. The data is normalized as follows. If the prior is continuous: X_normalized=(X-d)/s, where d=(min(X)+max(X))*0.5 and s=max(X)-min(X),
#' If the prior is discrete: X_normalized=X/s, where s=10^m, where m is chosen so that |max(abs(X))-1| is minimal. The outcome data is normalized as follows: Y_normalized=(Y-mu)/sd, where mu=mean(Y) and sd=sd(Y).
#' The priors inside the model are specified for the normalized data, in extreme cases not normalizing the data may lead to unreliable results; default is normalize=TRUE
#' @param n.chains - specifies the number of chains in the MCMC sampler; default is n.chains=4
#' @param burnin - specifies the number of burnin iterations without the adaptive iterations; default is burnin=5000
#' @param sample - specifies the number of samples per chain; default is samples=5000
#' @param adapt - specifies the number of adaptive iterations in the MCMC sampler; default is adapt=1000
#' @param inits - initial values for the sampler. By default the initial values are sampled inside the function. If user wants to provide
#' its own values please refer to run.jags manual; default inits=NULL
#' @param method - set to default as 'parallel', which allows to sample the chains in parallel reducing computation time.
#' To read more about possible method values type ?run.jags; default method='parallel'
#' @param seed - specifies the seed for replicability purposes; default is seed=NULL
#' @param ... - other arguments of run.jags function. For more details type ?run.jags
#' @return  The function returns the list with the elements:
#' \itemize{
#'  \item Effect_estimate: contains a list with MAP estimate and HDI of the treatment effect, the cutoff location (if unknown) and the discontinuity size in the treatment probability function (compliance rate at c) on the original, unnormalized scale;
#'  \item JAGS_output: contains output of the run.jags function for the normalized data if normalize=TRUE, based on this output mixing of the chains can be assessed;
#'  \item Samples: contains posterior samples of the treatment effect (eff), cutoff location (c) if unknown, and compliance rate (j);
#'  \item Normalized_data: contains a list of the normalized data (if normalized=TRUE) and the parameters used to normalize the data (see arg normalize);
#'  \item Priors: contains a list of the priors' parameters ;
#'  \item Inits contains the list of initial values  and .RNG.seed value
#' }
#' @import stats
#' @import runjags
#' @import bayestestR
#' @export
#' @examples
#' # code to generate the data
#'
#' invlogit <- function(x) {
#'   return(1 / (1 + exp(-x)))
#' }
#'
#' fun_prob55 <- function(X) {
#'   P = rep(0, length(X))
#'   P[X >= 0.] = invlogit((8.5 * X[X >= 0.] - 1.5)) / 10.5 + 0.65 - 0.0007072
#'   P[X < 0.] = (X[X < 0.] + 1)^4 / 15 + 0.05
#'   return(P)
#' }
#'
#' sample_prob55 <- function(X) {
#'   P = rep(0, length(X))
#'   P[X >= 0.] = invlogit((8.5 * X[X >= 0.] - 1.5)) / 10.5 + 0.65 - 0.0007072
#'   P[X < 0.] = (X[X < 0.] + 1)^4 / 15 + 0.05
#'   T = rep(0, length(X))
#'   for (j in 1:length(X)) {
#'     T[j] = sample(c(1, 0), 1, prob = c(P[j], 1 - P[j]))
#'   }
#'   return(T)
#' }
#'
#' funB <- function(X) {
#'   Y = X
#'   X2 = X[X >= 0]
#'   X1 = X[X < 0]
#'   Y[X < 0] = 1 / (1 + exp(-2 * X1)) - 0.5 + 0.4
#'   Y[X >= 0] = (log(X2 * 2 + 1) - 0.15 * X2^2) * 0.6 - 0.20 + 0.4
#'   return(Y)
#' }
#'
#' funB_sample <- function(X) {
#'   Y = funB(X)+ rnorm(length(X), 0, 0.025)
#'   return(Y)
#' }
#'
#' # data generation
#' N=150 # try different dataset size
#' set.seed(1234)
#' X = sort(runif(N, -1, 1))
#' T = sample_prob55(X)
#' Y = funB_sample(X)
#'
#' # comment out to try different priors:
#' c_prior=0 # known cutoff c=0
#' # c_prior=list(clb=-0.25,cub=0.25) # uniform prior on the interval [-0.25,0.25]
#' # c_prior=list(cstart=-0.25,cend=0.25,grid=0.05) # uniform discrete prior on -0.25, -0.2, ..., 0.25
#'
#' # running LoTTA model on fuzzy RDD with continuous outcomes;
#' # cutoff = 0, compliance rate = 0.55, treatment effect = -0.3636364
#' # remember to check convergence and adjust burnin, sample and adapt if needed
#' out = LoTTA_fuzzy_CONT(X,T,Y,c_prior,burnin = 3000,sample = 1000,adapt=1000,seed=1234)



LoTTA_fuzzy_CONT<-function(X,T,Y,c_prior,jlb=0.2,ci=0.95,trimmed=NULL,outcome_prior=list('pr'=0.0001,'df'=7),n_min=25,param=c('c','j','kl','kr','eff','a0l','a1l','a2l','a3l','a0r','a1r','a2r','a3r','b1lt','a1lt','a2lt','b2lt','b1rt','a1rt','a2rt','b2rt','k1t','k2t','tau1r','tau2r','tau1l','tau2l')
                           ,normalize=TRUE,n.chains=4,burnin=10000,sample=1000,adapt=500,inits=NULL,method='parallel',seed=NULL,...){

  if(all(unique(Y)%in%c(0,1))==TRUE){
    warning('Binary outcome detected. Using LoTTA_fuzzy_BIN is recommended. ')
  }
  if(is.null(trimmed)==TRUE){
    b=bounds(X,n_min)
    ubr=b$ubr
    ubl=b$ubl
    lb=b$lb
    prior=read_prior(c_prior,min(X),max(X),ubl,ubr,lb)
  }
  else{
    if(is.vector(trimmed)==TRUE){
      if(length(unique(trimmed))==2){
        Xtr=X[X>=min(trimmed)&X<=max(trimmed)]
        if(length(Xtr)==0){
          stop('there is no data in the trimmed dataset. Change the values of trimmed')
        }
        else{
          if(length(Xtr)<100){
            warning('there is less than 100 datapoints in the trimmed dataset')
          }
          b=bounds(Xtr,n_min)
          ubr=b$ubr
          ubl=b$ubl
          lb=b$lb
          prior=read_prior(c_prior,min(Xtr),max(Xtr),ubl,ubr,lb)
        }
      }
      else{
        stop('provided trimmed variable is neither NULL nor a vector with precisely two unique values')
      }
    }
    else{
      stop('wrong variable type. trimmed must be either NULL or a vector with precisely two unique values')
    }
  }
  normalized_parameters=list(X=X,Y=Y,T=T,d_X=0,s_X=1,mu_Y=0,sd_Y=1)
  if(normalize==TRUE || is.null(trimmed)==FALSE){
    if(prior$type=='D'){
      x_list=normalize_dis_X(X,trimmed,normalize)

    }
    else{
      x_list=normalize_cont_X(X,trimmed,normalize)

    }

    y_list=normalize_cont_Y(Y,X,trimmed,normalize)
    t_list=normalize_cont_Y(T,X,trimmed,normalize=FALSE)
    X=x_list$X
    Y=y_list$Y
    T=t_list$Y
    normalized_parameters=list(X=X,Y=Y,T=T,d_X=x_list$d,s_X=x_list$s,mu_Y=y_list$mu,sd_Y=y_list$sd)

  }

  pr=outcome_prior$pr
  df=outcome_prior$df
  seed=ifelse(is.null(seed)==TRUE,sample(1:1000000,1),seed)


  if(prior$type=='C'){

    alpha=prior$alpha
    beta=prior$beta
    clb=prior$clb
    cub=prior$cub
    normalized_parameters$clb=clb
    normalized_parameters$cub=cub
    if(normalize==TRUE){
      clb=(clb-normalized_parameters$d_X)/normalized_parameters$s_X
      cub=(cub-normalized_parameters$d_X)/normalized_parameters$s_X
      normalized_parameters$clb=clb
      normalized_parameters$cub=cub
    }
    b=bounds(X,n_min)
    ubr=b$ubr
    ubl=b$ubl
    ubrt=b$ubr
    ublt=b$ubl
    lb=b$lb

    datT=list(N=length(X),x=X,t=T,jlb=jlb,clb=clb,cub=cub,alpha=alpha,beta=beta)
    param_c=c('c0')
    c0=ifelse(alpha>1 & beta>1,(alpha-1)/(alpha+beta-2),0.5)
    initc1=list(c0=c0,al=0.5-0.5*((1+jlb)*0.5),j=(1+jlb)*0.5,.RNG.seed=seed,.RNG.name="base::Mersenne-Twister")
    model_path1 <- system.file("models", 'cutoff_initial_CONT.txt', package = "LoTTA")
    suppressWarnings(dat1_c<- run.jags(model_path1,inits = list(initc1) ,data=datT,monitor=param_c,burnin = 1000,sample=500,adapt = 100,n.chains = 1,method = 'simple'))
    C_start=as.numeric(combine.mcmc(dat1_c$mcmc))
    if(is.null(inits)==TRUE){
      inits=list()
      for (i in 1:n.chains) {

        suppressWarnings(inits[[i]]<-Initial_CONT_CONT(X,T,Y,C_start,clb,cub,lb,ubr,ubl,seed+i,jlb))
      }
    }
    dat=list(N=length(X),x=X,t=T,y=Y,ubr=ubr,ubl=ubl,ubrt=ubrt,ublt=ublt,lb=lb,jlb=jlb,clb=clb,cub=cub,alpha=alpha,beta=beta,pr=pr,df=df,seed=seed)
    model_path2 <- system.file("models", 'LoTTA_CONT_CONT.txt', package = "LoTTA")
    posterior<- tryCatch({run.jags(model_path2, data=dat,inits = inits,monitor=param,burnin = burnin,sample=sample,adapt = adapt,n.chains = n.chains,method = method,...)},error = function(e) {
      # Handle the error and return a specific value
      message("An error occurred, only initial values are returned: ", e$message)
      return(NA)
      # Return a specific value on error
    })
    if (identical(posterior, NA)) {
      return(inits)  # Exit the entire function
    }
    Samples=combine.mcmc(posterior)
    C=as.numeric(Samples[,'c'])*normalized_parameters$s_X+normalized_parameters$d_X
    mapc=as.numeric(map_estimate(C))
    hdic=as.numeric(ci(C,ci=ci,method='HDI'))[2:3]
  }

  if(prior$type=='D'){

    cstart=prior$cstart
    cend=prior$cend
    grid=prior$grid
    weights=prior$weights
    normalized_parameters$cstart=cstart
    normalized_parameters$cend=cend
    normalized_parameters$grid=grid
    if(normalize==TRUE){
      cstart=(cstart-normalized_parameters$d_X)/normalized_parameters$s_X
      cend=(cend-normalized_parameters$d_X)/normalized_parameters$s_X
      grid=grid/normalized_parameters$s_X

      normalized_parameters$cstart=cstart
      normalized_parameters$cend=cend
      normalized_parameters$grid=grid
    }
    b=bounds(X,n_min)
    ubr=b$ubr
    ubl=b$ubl
    ubrt=b$ubr
    ublt=b$ubl
    lb=b$lb

    datT=list(N=length(X),x=X,t=T,jlb=jlb,cstart=cstart,grid=grid,weights=weights,lb=lb)
    param_c=c('ct')
    ct=ifelse(length(unique(weights))==1,ceiling(length(weights)*0.5),which.max(weights))
    initc1=list(ct=ct,al=0.5-0.5*((1+jlb)*0.5),j=(1+jlb)*0.5,.RNG.seed=seed,.RNG.name="base::Mersenne-Twister")
    model_path1 <- system.file("models", 'cutoff_initial_DIS.txt', package = "LoTTA")
    suppressWarnings(dat1_c<- run.jags(model_path1,inits = list(initc1) ,data=datT,monitor=param_c,burnin = 1000,sample=500,adapt = 100,n.chains = 1,method = 'simple'))
    C_start=as.numeric(combine.mcmc(dat1_c$mcmc))

    if(is.null(inits)==TRUE){
      inits=list()
      for (i in 1:n.chains) {

        suppressWarnings(inits[[i]]<-Initial_DIS_CONT(X,T,Y,C_start,cstart,grid,lb,ubr,ubl,seed+i,jlb))
      }
    }

    dat=list(N=length(X),x=X,t=T,y=Y,ubr=ubr,ubl=ubl,ubrt=ubrt,ublt=ublt,lb=lb,jlb=jlb,cstart=cstart,grid=grid,weights=weights,pr=pr,df=df,seed=seed)
    model_path2 <- system.file("models", 'LoTTA_DIS_CONT.txt', package = "LoTTA")
    posterior<- tryCatch({run.jags(model_path2, data=dat,inits = inits,monitor=param,burnin = burnin,sample=sample,adapt = adapt,n.chains = n.chains,method = method,...)},error = function(e) {
      # Handle the error and return a specific value
      message("An error occurred, only initial values are returned: ", e$message)
      return(NA)
      # Return a specific value on error
    })
    if (identical(posterior, NA)) {
      return(inits)  # Exit the entire function
    }
    Samples=combine.mcmc(posterior)
    C=as.numeric(Samples[,'c'])*normalized_parameters$s_X+normalized_parameters$d_X
    mapc=which.max(table(C))
    mapc=as.numeric(names(mapc))
    hdic=as.numeric(ci(C,ci=ci,method='HDI'))[2:3]
  }


  if(prior$type=='F'){
    c=c_prior
    normalized_parameters$c=c
    if(normalize==TRUE || is.null(trimmed)==FALSE){
      c=(c-x_list$d)/x_list$s
      normalized_parameters$c=c
    }
    b=bounds(X,25)
    ubr=b$ubr
    ubl=b$ubl
    ubrt=b$ubr
    ublt=b$ubl
    lb=b$lb
    param=param[!(param %in% c('c'))]
    if(is.null(inits)==TRUE){
      inits=list()
      for (i in 1:n.chains) {

        inits[[i]]=Initial_FUZZY_CONT(X,T,Y,c,lb,ubr,ubl,seed+i,jlb)
      }
    }
    dat=list(N=length(X),x=X,t=T,y=Y,c=c,ubr=ubr,ubl=ubl,ubrt=ubrt,ublt=ublt,lb=lb,jlb=jlb,pr=pr,df=df,seed=seed)
    model_path <- system.file("models", 'LoTTA_CONT_c.txt', package = "LoTTA")
    posterior<- tryCatch({run.jags(model_path, data=dat,inits = inits,monitor=param,burnin = burnin,sample=sample,adapt = adapt,n.chains = n.chains,method = method,...)},error = function(e) {
      # Handle the error and return a specific value
      message("An error occurred, only initial values are returned: ", e$message)
      return(NA)
      # Return a specific value on error
    })
    if (identical(posterior, NA)) {
      return(inits)  # Exit the entire function
    }
    Samples=combine.mcmc(posterior)
  }

  Eff=as.numeric(Samples[,'eff'])*normalized_parameters$sd_Y
  mapeff=as.numeric(map_estimate(Eff))
  hdieff=as.numeric(ci(Eff,ci=ci,method='HDI'))[2:3]
  J=as.numeric(Samples[,'j'])
  mapj=as.numeric(map_estimate(J))
  hdij=as.numeric(ci(J,ci=ci,method='HDI'))[2:3]


  if(is.numeric(c_prior)==FALSE){
    Samples=data.frame('eff'=Eff,'c'=C,'j'=J)
    Res=list(Effect_estimate=list('MAP treatment effect'=mapeff,'CI treatment effect'=hdieff,'MAP cutoff'=mapc,'CI cutoff'=hdic,'MAP compliance rate'=mapj,'CI compliance rate'=hdij),JAGS_output=posterior,Samples=Samples,Normalized_data=normalized_parameters,Priors=list('Cutoff prior'= prior[-1],'Outcome prior'=list(pr=pr,df=df)),Inits=inits)
  }
  else{
    Samples=data.frame('eff'=Eff,'j'=J)
    Res=list(Effect_estimate=list('MAP treatment effect'=mapeff,'CI treatment effect'=hdieff,'MAP compliance rate'=mapj,'CI compliance rate'=hdij),JAGS_output=posterior,Samples=Samples,Normalized_data=normalized_parameters,Priors=list('Cutoff'=c_prior,'Outcome prior'=list(pr=pr,df=df)),Inits=inits)
  }

  return(Res)


}


#' LoTTA_treatment
#'
#' Function that fits LoTTA treatment model to the fuzzy RD treatment data with an either known or unknown/suspected cutoff.
#' It supports two types of priors on the cutoff location: a scaled beta distribution of the form beta(alpha,beta)(cub-clb)+clb and a discrete distribution with the support of the form cstart+grid i for i=0,...,(cend-cstart)/grid.
#' The score does NOT have to be normalized beforehand. We recommend NOT to transform the data before imputing it into the function, except for initial trimming of the score which should be done beforehand.
#' The further trimming for the sensitivity analysis can be done through the function, which ensures that the data is normalized before the trimming.
#'
#' @param X - is the score data
#' @param T - is the treatment allocation data
#' @param c_prior - specifies the cutoff prior in case of the unknown cutoff or the cutoff point
#' if the cutoff is known. Takes as value a number if the cutoff is known or a list of values otherwise.
#' For a continuous prior the list requires the following elements: clb - left end of the interval cub - right end of the interval
#' in which the scaled and translated beta distribution is defined, alpha (optional) - shape parameter, default value = 1, beta (optional) - shape parameter, default value = 1.
#' For a discrete prior the list requires the following elements: cstart - first point with positive prior mass, cend - last point with positive prior mass, grid - distance between the consecutive points in the support
#' weights (optional) - vector of weights assigned to each point in the support, default is vector of 1's (uniform distribution)
#' @param jlb - minimum jump size
#' @param ci - specifies the probability level 1-alpha for the highest posterior density intervals; default is ci = 0.95
#' @param trimmed - takes as a value NULL or a vector of two values. It specifies potential trimming of the data.
#' If set to NULL no trimming is applied to the data. If a list of two values is provided the data is trimmed
#' to include data points with the score X in between those values; default is trimmed=NULL
#' @param n_min - specifies the minimum number of data points to which a cubic part of the outcome function is fit to ensure stability of the sampling procedure; default is n_min=25
#' @param param - takes as a value a vector with names of the parameters that are to be sampled; default is the list of all parameters
#' @param normalize - specifies if the data is to be normalized. The data is normalized as follows. If the prior is continuous: X_normalized=(X-d)/s, where d=(min(X)+max(X))*0.5 and s=max(X)-min(X),
#' If the prior is discrete: X_normalized=X/s, where s=10^m, where m is chosen so that |max(abs(X))-1| is minimal. The outcome data is normalized as follows: Y_normalized=(Y-mu)/sd, where mu=mean(Y) and sd=sd(Y).
#' The priors inside the model are specified for the normalized data, in extreme cases not normalizing the data may lead to unreliable results; default is normalize=TRUE
#' @param n.chains - specifies the number of chains in the MCMC sampler; default is n.chains=4
#' @param burnin - specifies the number of burnin iterations without the adaptive iterations; default is burnin=5000
#' @param sample - specifies the number of samples per chain; default is samples=5000
#' @param adapt - specifies the number of adaptive iterations in the MCMC sampler; default is adapt=1000
#' @param inits - initial values for the sampler. By default the initial values are sampled inside the function. If user wants to provide
#' its own values please refer to run.jags manual; default inits=NULL
#' @param method - set to default as 'parallel', which allows to sample the chains in parallel reducing computation time.
#' To read more about possible method values type ?run.jags; default method='parallel'
#' @param seed - specifies the seed for replicability purposes; default is seed=NULL
#' @param ... - other arguments of run.jags function. For more details type ?run.jags
#' @return  The function returns the list with the elements:
#' \itemize{
#'  \item Effect_estimate: contains a list with MAP estimate and HDI of the cutoff location (if unknown) and the discontinuity size in the treatment probability function (compliance rate at c) on the original, unnormalized scale;
#'  \item JAGS_output: contains output of the run.jags function for the normalized data if normalize=TRUE, based on this output mixing of the chains can be assessed;
#'  \item Samples: contains posterior samples of the cutoff location (c) if unknown, and compliance rate (j);
#'  \item Normalized_data: contains a list of the normalized data (if normalized=TRUE) and the parameters used to normalize the data (see arg normalize);
#'  \item Priors: contains a list of the priors' parameters ;
#'  \item Inits contains the list of initial values  and .RNG.seed value
#' }

#' @import stats
#' @import runjags
#' @import bayestestR
#' @export
#' @examples
#' # code to generate the data
#'
#' invlogit <- function(x) {
#'   return(1 / (1 + exp(-x)))
#' }
#'
#' fun_prob55 <- function(X) {
#'   P = rep(0, length(X))
#'   P[X >= 0.] = invlogit((8.5 * X[X >= 0.] - 1.5)) / 10.5 + 0.65 - 0.0007072
#'   P[X < 0.] = (X[X < 0.] + 1)^4 / 15 + 0.05
#'   return(P)
#' }
#'
#' sample_prob55 <- function(X) {
#'   P = rep(0, length(X))
#'   P[X >= 0.] = invlogit((8.5 * X[X >= 0.] - 1.5)) / 10.5 + 0.65 - 0.0007072
#'   P[X < 0.] = (X[X < 0.] + 1)^4 / 15 + 0.05
#'   T = rep(0, length(X))
#'   for (j in 1:length(X)) {
#'     T[j] = sample(c(1, 0), 1, prob = c(P[j], 1 - P[j]))
#'   }
#'   return(T)
#' }
#'
#' # data generation
#' N=150 # try different dataset size
#' set.seed(1234)
#' X = sort(runif(N, -1, 1))
#' T = sample_prob55(X)
#'
#' # comment out to try different priors:
#' c_prior=0 # known cutoff c=0
#' # c_prior=list(clb=-0.25,cub=0.25) # uniform prior on the interval [-0.25,0.25]
#' # c_prior=list(cstart=-0.25,cend=0.25,grid=0.05) # uniform discrete prior on -0.25, -0.2, ..., 0.25
#'
#' # running LoTTA model on the treatment allocation data;
#' # cutoff = 0, compliance rate = 0.55
#' # remember to check convergence and adjust burnin, sample and adapt if needed
#' out = LoTTA_treatment(X, T, c_prior, burnin = 3000, sample = 2000, adapt = 1000, seed = 1234)
#'

LoTTA_treatment<-function(X,T,c_prior,jlb=0.2,ci=0.95,trimmed=NULL,n_min=25,param=c('c','j','b1lt','a1lt','a2lt','b2lt','b1rt','a1rt','a2rt','b2rt','k1t','k2t')
                          ,normalize=TRUE,n.chains=4,burnin=5000,sample=1000,adapt=500,inits=NULL,method='parallel',seed=NULL,...){

  if(is.null(trimmed)==TRUE){
    b=bounds(X,n_min)
    ubr=b$ubr
    ubl=b$ubl
    lb=b$lb
    prior=read_prior(c_prior,min(X),max(X),ubl,ubr,lb)
  }
  else{
    if(is.vector(trimmed)==TRUE){
      if(length(unique(trimmed))==2){
        Xtr=X[X>=min(trimmed)&X<=max(trimmed)]
        if(length(Xtr)==0){
          stop('there is no data in the trimmed dataset. Change the values of trimmed')
        }
        else{
          if(length(Xtr)<100){
            warning('there is less than 100 datapoints in the trimmed dataset')
          }
          b=bounds(Xtr,n_min)
          ubr=b$ubr
          ubl=b$ubl
          lb=b$lb
          prior=read_prior(c_prior,min(Xtr),max(Xtr),ubl,ubr,lb)
        }
      }
      else{
        stop('provided trimmed variable is neither NULL nor a vector with precisely two unique values')
      }
    }
    else{
      stop('wrong variable type. trimmed must be either NULL or a vector with precisely two unique values')
    }
  }
  seed=ifelse(is.null(seed)==TRUE,sample(1:1000000,1),seed)
  normalized_parameters=list(X=X,T=T,d_X=0,s_X=1)
  if(normalize==TRUE || is.null(trimmed)==FALSE){
    if(prior$type=='D'){
      x_list=normalize_dis_X(X,trimmed,normalize)

    }
    else{
      x_list=normalize_cont_X(X,trimmed,normalize)

    }

    t_list=normalize_cont_Y(T,X,trimmed,normalize=FALSE)
    X=x_list$X
    T=t_list$Y
    normalized_parameters=list(X=X,T=T,d_X=x_list$d,s_X=x_list$s)

  }
  if(prior$type=='C'){
    alpha=prior$alpha
    beta=prior$beta
    clb=prior$clb
    cub=prior$cub
    normalized_parameters$clb=clb
    normalized_parameters$cub=cub
    if(normalize==TRUE){
      clb=(clb-normalized_parameters$d_X)/normalized_parameters$s_X
      cub=(cub-normalized_parameters$d_X)/normalized_parameters$s_X
      normalized_parameters$clb=clb
      normalized_parameters$cub=cub
    }
    b=bounds(X,n_min)
    ubrt=b$ubr
    ublt=b$ubl
    lb=b$lb

    datT=list(N=length(X),x=X,t=T,jlb=jlb,ubrt=ubrt,ublt=ublt,clb=clb,cub=cub,alpha=alpha,beta=beta)
    param_c=c('c0')
    c0=ifelse(alpha>1 & beta>1,(alpha-1)/(alpha+beta-2),0.5)
    initc1=list(c0=c0,al=0.5-0.5*((1+jlb)*0.5),j=(1+jlb)*0.5,.RNG.seed=seed,.RNG.name="base::Mersenne-Twister")
    model_path1 <- system.file("models", 'cutoff_initial_CONT.txt', package = "LoTTA")
    suppressWarnings(dat1_c<- run.jags(model_path1,inits = list(initc1) ,data=datT,monitor=param_c,burnin = 1000,sample=500,adapt = 100,n.chains = 1,method = 'simple'))
    C_start=as.numeric(combine.mcmc(dat1_c$mcmc))
    if(is.null(inits)==TRUE){
      inits=list()
      for (i in 1:n.chains) {

        suppressWarnings(inits[[i]]<-Initial_treatment_CONT(X,T,C_start,clb,cub,lb,ubrt,ublt,seed+i,jlb))
      }
    }
    dat=list(N=length(X),x=X,t=T,ubrt=ubrt,ublt=ublt,lb=lb,jlb=jlb,clb=clb,cub=cub,alpha=alpha,beta=beta,seed=seed)
    model_path2 <- system.file("models", 'LoTTA_treatment_CONT.txt', package = "LoTTA")
    posterior<- tryCatch({run.jags(model_path2, data=dat,inits = inits,monitor=param,burnin = burnin,sample=sample,adapt = adapt,n.chains = n.chains,method = method,...)},error = function(e) {
      # Handle the error and return a specific value
      message("An error occurred, only initial values are returned: ", e$message)
      return(NA)
      # Return a specific value on error
    })
    if (identical(posterior, NA)) {
      return(inits)  # Exit the entire function
    }
    Samples=combine.mcmc(posterior)
    C=as.numeric(Samples[,'c'])*normalized_parameters$s_X+normalized_parameters$d_X
    mapc=as.numeric(map_estimate(C))
    hdic=as.numeric(ci(C,ci=ci,method='HDI'))[2:3]
  }
  if(prior$type=='D'){
    cstart=prior$cstart
    cend=prior$cend
    grid=prior$grid
    weights=prior$weights
    normalized_parameters$cstart=cstart
    normalized_parameters$cend=cend
    normalized_parameters$grid=grid
    if(normalize==TRUE){
      cstart=(cstart-normalized_parameters$d_X)/normalized_parameters$s_X
      cend=(cend-normalized_parameters$d_X)/normalized_parameters$s_X
      grid=grid/normalized_parameters$s_X

      normalized_parameters$cstart=cstart
      normalized_parameters$cend=cend
      normalized_parameters$grid=grid
    }
    b=bounds(X,n_min)
    ubrt=b$ubr
    ublt=b$ubl
    lb=b$lb

    datT=list(N=length(X),x=X,t=T,jlb=jlb,cstart=cstart,grid=grid,weights=weights,lb=lb)
    param_c=c('ct')
    ct=ifelse(length(unique(weights))==1,ceiling(length(weights)*0.5),which.max(weights))
    initc1=list(ct=ct,al=0.5-0.5*((1+jlb)*0.5),j=(1+jlb)*0.5,.RNG.seed=seed,.RNG.name="base::Mersenne-Twister")
    model_path1 <- system.file("models", 'cutoff_initial_DIS.txt', package = "LoTTA")
    suppressWarnings(dat1_c<- run.jags(model_path1,inits = list(initc1) ,data=datT,monitor=param_c,burnin = 1000,sample=500,adapt = 100,n.chains = 1,method = 'simple'))
    C_start=as.numeric(combine.mcmc(dat1_c$mcmc))
    if(is.null(inits)==TRUE){
      inits=list()
      for (i in 1:n.chains) {

        suppressWarnings(inits[[i]]<-Initial_treatment_DIS(X,T,C_start,cstart,grid,lb,ubrt,ublt,seed+i,jlb))
      }
    }
    dat=list(N=length(X),x=X,t=T,ubrt=ubrt,ublt=ublt,lb=lb,jlb=jlb,cstart=cstart,grid=grid,weights=weights,seed=seed)
    model_path2 <- system.file("models", 'LoTTA_treatment_DIS.txt', package = "LoTTA")
    posterior<- tryCatch({run.jags(model_path2, data=dat,inits = inits,monitor=param,burnin = burnin,sample=sample,adapt = adapt,n.chains = n.chains,method = method,...)},error = function(e) {
      # Handle the error and return a specific value
      message("An error occurred, only initial values are returned: ", e$message)
      return(NA)
      # Return a specific value on error
    })
    if (identical(posterior, NA)) {
      return(inits)  # Exit the entire function
    }
    Samples=combine.mcmc(posterior)
    C=as.numeric(Samples[,'c'])*normalized_parameters$s_X+normalized_parameters$d_X
    mapc=which.max(table(C))
    mapc=as.numeric(names(mapc))
    hdic=as.numeric(ci(C,ci=ci,method='HDI'))[2:3]
  }
  if(prior$type=='F'){
    c=c_prior
    normalized_parameters$c=c
    if(normalize==TRUE || is.null(trimmed)==FALSE){
      c=(c-x_list$d)/x_list$s
      normalized_parameters$c=c
    }
    b=bounds(X,25)
    ubrt=b$ubr
    ublt=b$ubl
    lb=b$lb
    param=param[!(param %in% c('c'))]
    if(is.null(inits)==TRUE){
      inits=list()
      for (i in 1:n.chains) {

        inits[[i]]=Initial_treatment_c(X,T,c,lb,ubrt,ublt,seed+i,jlb)
      }
    }
    dat=list(N=length(X),x=X,t=T,c=c,ubrt=ubrt,ublt=ublt,lb=lb,jlb=jlb,seed=seed)
    model_path <- system.file("models", 'LoTTA_treatment_c.txt', package = "LoTTA")
    posterior<- tryCatch({run.jags(model_path, data=dat,inits = inits,monitor=param,burnin = burnin,sample=sample,adapt = adapt,n.chains = n.chains,method = method,...)},error = function(e) {
      # Handle the error and return a specific value
      message("An error occurred, only initial values are returned: ", e$message)
      return(NA)
      # Return a specific value on error
    })
    if (identical(posterior, NA)) {
      return(inits)  # Exit the entire function
    }
    Samples=combine.mcmc(posterior)
  }
  J=as.numeric(Samples[,'j'])
  mapj=as.numeric(map_estimate(J))
  hdij=as.numeric(ci(J,ci=ci,method='HDI'))[2:3]
  if(is.numeric(c_prior)==FALSE){
    Samples=data.frame('c'=C,'j'=J)
    Res=list(Effect_estimate=list('MAP cutoff'=mapc,'CI cutoff'=hdic,'MAP compliance rate'=mapj,'CI compliance rate'=hdij),JAGS_output=posterior,Samples=Samples,Normalized_data=normalized_parameters,Priors=list('Cutoff prior'=prior[-1]),Inits=inits)
  }
  else{
    Samples=data.frame('j'=J)
    Res=list(Effect_estimate=list('MAP compliance rate'=mapj,'CI compliance rate'=hdij),JAGS_output=posterior,Samples=Samples,Normalized_data=normalized_parameters,Priors=list('Cutoff'=c_prior),Inits=inits)
  }

  return(Res)
}
