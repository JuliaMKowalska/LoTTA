#' Function that splits the data into bins and computes the average in each bin
#' @param Y - is the outcome data
#' @param X - is the score data
#' @param c - specifies the cutoff point, set to NULL if the binning of the data should be independent from c
#' @param binsize - length of the interval that is one bin
#' @return list with elements Y_bin: list of average values of Y in each bin, X_bin: list of middle points for each bin

Bin_data<-function(Y,X,c=NULL,binsize=0.2){
  if(is.null(c)==1){
    L=seq(min(X),max(X),binsize)
  }
  else{
    L=seq(min(X),c,binsize)
    L=append(L,seq(c,max(X),binsize))
  }

  L=append(L,c(max(X)+0.0000001))
  L=unique(L)
  groups <- cut(X, L,right=FALSE,labels = L[-length(L)])
  Group<-as.numeric(groups)
  Binned_X=as.numeric(lapply(split(X, groups), mean))
  Binned=lapply(split(Y, groups), mean)
  Binned=as.numeric(Binned)
  return(list('Y_binned'=Binned,'X_binned'= Binned_X))
}

#' Function that evaluates the binary outcome function in a domain x, given the coefficients
#' @param coef_s - list with the function coefficients
#' @param x - points at which the function is evaluated
#' @param d_x - shifting parameter that was used to normalize X
#' @param s_x - scaling parameter that was used to normalize X
#' @return list with the values of the function

BIN_outcome_function_sample<-function(coef_s,x,d_x,s_x){
  x=x
  c=coef_s[['c']]
  a0l=coef_s[['a0l']]
  a1l=coef_s[['a1l']]
  a2l=coef_s[['a2l']]
  a3l=coef_s[['a3l']]

  a0r=coef_s[['a0r']]
  a1r=coef_s[['a1r']]
  a2r=coef_s[['a2r']]
  a3r=coef_s[['a3r']]


  b0l=logit(a0l)
  b1l=a1l*(invlogit(-b0l)*(1-invlogit(-b0l)))^(-1)

  b0r=logit(a0r)
  b1r=a1r*(invlogit(-b0r)*(1-invlogit(-b0r)))^(-1)

  kl=coef_s[['kl']]
  kr=coef_s[['kr']]
  xc=x-c
  param <-ifelse(x<c,a1l*(xc)+a0l+invlogit(100*(kl-x))*(invlogit((a3l)*(xc)^3+(a2l)*(xc)^2+b1l*(xc)+b0l)-a0l-a1l*xc),a1r*(xc)+a0r+invlogit(100*(x-kr))*(invlogit((a3r)*(xc)^3+(a2r)*(xc)^2+b1r*(xc)+b0r)-a0r-a1r*xc))
  return(param)
}

#' Function that evaluates the continuous outcome function in a domain x, given the coefficients
#' @param coef_s - list with the function coefficients
#' @param x - points at which the function is evaluated
#' @param d_x - shifting parameter that was used to normalize X
#' @param s_x - scaling parameter that was used to normalize X
#' @param mu_y - shifting parameter that was used to normalize Y
#' @param sd_y - scaling parameter that was used to normalize Y
#' @return list with the values of the function

CONT_outcome_function_sample<-function(coef_s,x,d_x,s_x,mu_y,sd_y){
  x=(x-d_x)/s_x
  c=coef_s[['c']]
  a0l=coef_s[['a0l']]
  a1l=coef_s[['a1l']]
  a2l=coef_s[['a2l']]
  a3l=coef_s[['a3l']]

  a0r=coef_s[['a0r']]
  a1r=coef_s[['a1r']]
  a2r=coef_s[['a2r']]
  a3r=coef_s[['a3r']]


  kl=coef_s[['kl']]
  kr=coef_s[['kr']]
  xc=x-c
  param <-ifelse(x<c,a1l*(xc)+a0l+invlogit(100*(kl-x))*((a3l)*(xc)^3+(a2l)*(xc)^2),a1r*(xc)+a0r+invlogit(100*(x-kr))*((a3r)*(xc)^3+(a2r)*(xc)^2))
  param=param*sd_y+mu_y
  return(param)
}

#' Function that evaluates the treatment probability function in a domain x, given the coefficients
#' @param coef_s - list with the function coefficients
#' @param x - points at which the function is evaluated
#' @param d_x - shifting parameter that was used to normalize X
#' @param s_x - scaling parameter that was used to normalize X
#' @return list with the values of the function

treatment_function_sample<-function(coef_s,x,d_x,s_x){
  x=(x-d_x)/s_x
  c=coef_s[['c']]
  a1lt=coef_s[['a1lt']]
  a2lt=coef_s[['a2lt']]

  a1rt=coef_s[['a1rt']]
  a2rt=coef_s[['a2rt']]

  b1lt=coef_s[['b1lt']]
  b2lt=coef_s[['b2lt']]

  b1rt=coef_s[['b1rt']]
  b2rt=coef_s[['b2rt']]

  k1t=coef_s[['k1t']]
  k2t=coef_s[['k2t']]

  paramt<- ifelse(x<c,ifelse(x>=c-k1t,a1lt*x+b1lt,a2lt*x+b2lt),ifelse(x<=c+k2t,a1rt*x+b1rt,a2rt*x+b2rt))

  return(paramt)
}

#' Function that plots the median (or another quantile) of the posterior function
#' with binary outcome along with the quanatile-based credible interval.
#' The function is plotted on top of the binned input data. To bin the data, the score data
#' is divided into bins of fixed length, then the average outcome in each bin is calculated.
#' The average outcomes are plotted against the average values of the score in the corresponding bins.
#' @param LoTTA_posterior - output of one of the LoTTA functions (LoTTA_sharp, LoTTA_fuzzy)
#' @param nbins - number of bins to aggregate the input data
#' @param probs - list of three quantiles, the first and the last one define the quanatile-based credible interval,
#' the middle value defines the quantile of the posterior function to plot;
#' by default the quantiles correspond to the median posterior function and 95% credible interval
#' probs=c(0.025,0.5,0.975)
#' @param n_eval - n_eval*range(X) is the number of points at which each posterior function is evaluated,
#' the higher number means slower computing time and a smoother plot; default n_eval=200
#' @param col_line - the color of the line and the band
#' @param size_line - thickness of the line
#' @param alpha_interval - alpha value of the band, lower values correspond to a more transparent color
#' @param col_dots - color of the dots that correspond to the binned data
#' @param size_dots - size of the dots that correspond to the binned data
#' @param alpha_dots - transparency of the dots that correspond to the binned data,
#'  lower values correspond to a more transparent color
#' @param col_cutoff - color of the dotted line at the cutoff
#' @param title - title of the plot
#' @param subtitle - subtitle of the plot
#' @param y_lab - label of the y-axis
#' @param x_lab - label of the x-axis
#' @param plot.theme - ggplot2 plot theme (see https://ggplot2.tidyverse.org/reference/ggtheme.html)
#' possibly with additional arguments, it takes the default value plot.theme=theme_classic(base_size = 14),
#' @param legend.position - position of the legend, refer to ggplot2 manual for the possible values;
#' by default legend is not printed legend.position='none'
#' @param name_legend - title of the legend
#' @param labels_legend - the label of the plotted function in the legend
#' @param text - can be any value that is accepted in the argument _text_ in the _theme_ function
#' of ggplot2 package,refer to ggplot2 manual for the possible values;
#' by default is changes font to a serif one text=element_text(family='serif')
#' @param legend.text - can be any value that is accepted in the argument _legend.text_ in the _theme_ function
#' of ggplot2 package,refer to ggplot2 manual for the possible values;
#' by default is changes font size to the legend to 14 legend.text=element_text(size = 14)
#' @param plot.title - can be any value that is accepted in the argument _plot.title_ in the _theme_ function
#' of ggplot2 package,refer to ggplot2 manual for the possible values;
#' by default it centers the plot title plot.title = element_text(hjust = 0.5)
#' @param plot.subtitle - can be any value that is accepted in the argument _plot.subtitle_ in the _theme_ function
#' of ggplot2 package,refer to ggplot2 manual for the possible values;
#' by default it centers the plot subtitle plot.title = element_text(hjust = 0.5)
#' @param ... -  other arguments of the _theme_ function, refer to ggplot2 manual
#' @return ggplot2 object
#' @import ggplot2
#' @import ggpubr

plot_outcome_BIN<-function(LoTTA_posterior,nbins=100,probs=c(0.025,0.5,0.975),n_eval=200,col_line='#E69F00',size_line=0.1,alpha_interval=0.35,col_dots='gray',size_dots=3,alpha_dots=0.6,col_cutoff='black',
                           title="Outcome function",subtitle=NULL,y_lab="",x_lab=expression(paste(italic('X')," - score")) ,plot.theme=theme_classic(base_size = 14),
                           legend.position='none', name_legend='Legend', labels_legend = 'median outcome fun.',
                           text=element_text(family='serif'),legend.text=element_text(size = 14),
                           plot.title = element_text(hjust = 0.5),plot.subtitle=element_text(hjust = 0.5),...){
  LoTTA_jags=LoTTA_posterior$JAGS_output
  Samples_FULL=data.frame(data.frame(combine.mcmc(LoTTA_jags)))
  d_x=LoTTA_posterior$Normalized_data$d_X
  s_x=LoTTA_posterior$Normalized_data$s_X
  X=LoTTA_posterior$Normalized_data$'X'*s_x+d_x
  Y=LoTTA_posterior$Normalized_data$'Y'
  if(is.null(LoTTA_posterior$Effect_estimate$`MAP cutoff`)==TRUE){
    c=LoTTA_posterior$Normalized_data$c
    Samples_FULL[,'c']=rep(c,nrow(Samples_FULL))
  }
  else{
    c=LoTTA_posterior$Effect_estimate$`MAP cutoff`
  }
  rangeX=max(X)-min(X)
  binsize=rangeX/nbins
  binned_data=Bin_data(Y,X,c,binsize)
  X_plot=binned_data[[2]]
  Y_plot=binned_data[[1]]

  x=seq(min(X),max(X),length.out = n_eval*rangeX)
  fun_sampleF=apply(Samples_FULL,1,BIN_outcome_function_sample,x=x,d_x=d_x,s_x=s_x)
  q_funF=apply(fun_sampleF,1,quantile,probs=probs)
  plot.df=data.frame('X'=X_plot,'Y'=Y_plot)
  lower95=q_funF[1,]
  upper95=q_funF[3,]
  lineF.df=data.frame(X=x,Y=q_funF[2,],lower95=lower95,upper95=upper95)
  theme_def=plot.theme
  ggscatter(plot.df, x = "X", y = "Y",
            color = col_dots,
            size = size_dots, alpha = alpha_dots)+geom_vline(xintercept = c,linetype='dotted',color=col_cutoff)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    geom_ribbon(aes(x=X,y=Y,ymin = lower95, ymax = upper95,fill='Full'),data=lineF.df,alpha=alpha_interval)+geom_point(aes(x=X,y=Y),data=lineF.df,size=size_line,col=col_line)+
    labs( y=y_lab, x=x_lab,title=title,subtitle = subtitle)+theme_def+
    theme(text = text,legend.text =legend.text,plot.title = plot.title,plot.subtitle = plot.subtitle,legend.position = legend.position,...)+ scale_fill_manual(name=name_legend,values=c('Full'=col_line),labels = labels_legend)

}


#' Function that plots the median (or another quantile) of the posterior function of a continous outcome along with the quanatile-based credible interval.
#' The function is plotted on top of the binned input data. To bin the data, the score data
#' is divided into bins of fixed length, then the average outcome in each bin is calculated.
#' The average outcomes are plotted against the average values of the score in the corresponding bins.
#' @param LoTTA_posterior - output of one of the LoTTA functions (LoTTA_sharp, LoTTA_fuzzy)
#' @param nbins - number of bins to aggregate the input data
#' @param probs - list of three quantiles, the first and the last one define the quanatile-based credible interval,
#' the middle value defines the quantile of the posterior function to plot;
#' by default the quantiles correspond to the median posterior function and 95% credible interval
#' probs=c(0.025,0.5,0.975)
#' @param n_eval - n_eval*range(X) is the number of points at which each posterior function is evaluated,
#' the higher number means slower computing time and a smoother plot; default n_eval=200
#' @param col_line - the color of the line and the band
#' @param size_line - thickness of the line
#' @param alpha_interval - alpha value of the band, lower values correspond to a more transparent color
#' @param col_dots - color of the dots that correspond to the binned data
#' @param size_dots - size of the dots that correspond to the binned data
#' @param alpha_dots - transparency of the dots that correspond to the binned data,
#'  lower values correspond to a more transparent color
#' @param col_cutoff - color of the dotted line at the cutoff
#' @param title - title of the plot
#' @param subtitle - subtitle of the plot
#' @param y_lab - label of the y-axis
#' @param x_lab - label of the x-axis
#' @param plot.theme - ggplot2 plot theme (see https://ggplot2.tidyverse.org/reference/ggtheme.html)
#' possibly with additional arguments, it takes the default value plot.theme=theme_classic(base_size = 14),
#' @param legend.position - position of the legend, refer to ggplot2 manual for the possible values;
#' by default legend is not printed legend.position='none'
#' @param name_legend - title of the legend
#' @param labels_legend - the label of the plotted function in the legend
#' @param text - can be any value that is accepted in the argument _text_ in the _theme_ function
#' of ggplot2 package,refer to ggplot2 manual for the possible values;
#' by default is changes font to a serif one text=element_text(family='serif')
#' @param legend.text - can be any value that is accepted in the argument _legend.text_ in the _theme_ function
#' of ggplot2 package,refer to ggplot2 manual for the possible values;
#' by default is changes font size to the legend to 14 legend.text=element_text(size = 14)
#' @param plot.title - can be any value that is accepted in the argument _plot.title_ in the _theme_ function
#' of ggplot2 package,refer to ggplot2 manual for the possible values;
#' by default it centers the plot title plot.title = element_text(hjust = 0.5)
#' @param plot.subtitle - can be any value that is accepted in the argument _plot.subtitle_ in the _theme_ function
#' of ggplot2 package,refer to ggplot2 manual for the possible values;
#' by default it centers the plot subtitle plot.title = element_text(hjust = 0.5)
#' @param ... -  other arguments of the _theme_ function, refer to ggplot2 manual
#' @return ggplot2 object
#' @import ggplot2
#' @import ggpubr

plot_outcome_CONT<-function(LoTTA_posterior,nbins=100,probs=c(0.025,0.5,0.975),n_eval=200,col_line='#E69F00',size_line=0.1,alpha_interval=0.35,col_dots='gray',size_dots=3,alpha_dots=0.6,col_cutoff='black',
                            title="Outcome function",subtitle=NULL,y_lab="",x_lab=expression(paste(italic('X')," - score")) ,plot.theme=theme_classic(base_size = 14),
                            legend.position='none', name_legend='Legend', labels_legend = 'median outcome fun.',
                            text=element_text(family='serif'),legend.text=element_text(size = 14),
                            plot.title = element_text(hjust = 0.5),plot.subtitle=element_text(hjust = 0.5),...){
  LoTTA_jags=LoTTA_posterior$JAGS_output
  Samples_FULL=data.frame(data.frame(combine.mcmc(LoTTA_jags)))
  d_x=LoTTA_posterior$Normalized_data$d_X
  s_x=LoTTA_posterior$Normalized_data$s_X
  mu_y=LoTTA_posterior$Normalized_data$mu_Y
  sd_y=LoTTA_posterior$Normalized_data$sd_Y
  X=LoTTA_posterior$Normalized_data$'X'*s_x+d_x
  Y=LoTTA_posterior$Normalized_data$'Y'*sd_y+mu_y
  if(is.null(LoTTA_posterior$Effect_estimate$`MAP cutoff`)==TRUE){
    c=LoTTA_posterior$Normalized_data$c
    Samples_FULL[,'c']=rep(c,nrow(Samples_FULL))
  }
  else{
    c=LoTTA_posterior$Effect_estimate$`MAP cutoff`
  }
  rangeX=max(X)-min(X)
  binsize=rangeX/nbins
  binned_data=Bin_data(Y,X,c,binsize)
  X_plot=binned_data[[2]]
  Y_plot=binned_data[[1]]

  x=seq(min(X),max(X),length.out = n_eval*rangeX)
  fun_sampleF=apply(Samples_FULL,1,CONT_outcome_function_sample,x=x,d_x=d_x,s_x=s_x,mu_y=mu_y,sd_y=sd_y)
  q_funF=apply(fun_sampleF,1,quantile,probs=probs)
  plot.df=data.frame('X'=X_plot,'Y'=Y_plot)
  lower95=q_funF[1,]
  upper95=q_funF[3,]
  lineF.df=data.frame(X=x,Y=q_funF[2,],lower95=lower95,upper95=upper95)
  theme_def=plot.theme
  ggscatter(plot.df, x = "X", y = "Y",
            color = col_dots,
            size = size_dots, alpha = alpha_dots)+geom_vline(xintercept = c,linetype='dotted',color=col_cutoff)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    geom_ribbon(aes(x=X,y=Y,ymin = lower95, ymax = upper95,fill='Full'),data=lineF.df,alpha=alpha_interval)+geom_point(aes(x=X,y=Y),data=lineF.df,size=size_line,col=col_line)+
    labs( y=y_lab, x=x_lab,title=title,subtitle = subtitle)+theme_def+
    theme(text = text,legend.text =legend.text,plot.title = plot.title,plot.subtitle = plot.subtitle,legend.position = legend.position,...)+ scale_fill_manual(name=name_legend,values=c('Full'=col_line),labels = labels_legend)

}
#' LoTTA_plot_outcome
#'
#' Function that plots the median (or another quantile) of the LoTTA posterior outcome function along with the quanatile-based credible interval.
#' The function is plotted on top of the binned input data. To bin the data, the score data
#' is divided into bins of fixed length, then the average outcome in each bin is calculated.
#' The average outcomes are plotted against the average values of the score in the corresponding bins.
#' The data is binned separately on each side of the cutoff, the cutoff is marked on the plot
#' with a dotted line. In case of an unknown cutoff, the MAP estimate is used.
#' @param LoTTA_posterior - output of one of the LoTTA functions (LoTTA_sharp_CONT, LoTTA_fuzzy_CONT, LoTTA_sharp_BIN, LoTTA_fuzzy_BIN)
#' with all the parameters sampled (the default option in those functions)
#' @param nbins - number of bins to aggregate the input data
#' @param probs - list of three quantiles, the first and the last one define the quanatile-based credible interval,
#' the middle value defines the quantile of the posterior function to plot;
#' by default the quantiles correspond to the median posterior function and 95% credible interval
#' probs=c(0.025,0.5,0.975)
#' @param n_eval - n_eval*range(X) is the number of points at which each posterior function is evaluated,
#' the higher number means slower computing time and a smoother plot; default n_eval=200
#' @param col_line - the color of the line and the band
#' @param size_line - thickness of the line
#' @param alpha_interval - alpha value of the band, lower values correspond to a more transparent color
#' @param col_dots - color of the dots that correspond to the binned data
#' @param size_dots - size of the dots that correspond to the binned data
#' @param alpha_dots - transparency of the dots that correspond to the binned data,
#'  lower values correspond to a more transparent color
#' @param col_cutoff - color of the dotted line at the cutoff
#' @param title - title of the plot
#' @param subtitle - subtitle of the plot
#' @param y_lab - label of the y-axis
#' @param x_lab - label of the x-axis
#' @param plot.theme - ggplot2 plot theme (see https://ggplot2.tidyverse.org/reference/ggtheme.html)
#' possibly with additional arguments, it takes the default value plot.theme=theme_classic(base_size = 14),
#' @param legend.position - position of the legend, refer to ggplot2 manual for the possible values;
#' by default legend is not printed legend.position='none'
#' @param name_legend - title of the legend
#' @param labels_legend - the label of the plotted function in the legend
#' @param text - can be any value that is accepted in the argument _text_ in the _theme_ function
#' of ggplot2 package,refer to ggplot2 manual for the possible values;
#' by default is changes font to a serif one text=element_text(family='serif')
#' @param legend.text - can be any value that is accepted in the argument _legend.text_ in the _theme_ function
#' of ggplot2 package,refer to ggplot2 manual for the possible values;
#' by default is changes font size to the legend to 14 legend.text=element_text(size = 14)
#' @param plot.title - can be any value that is accepted in the argument _plot.title_ in the _theme_ function
#' of ggplot2 package,refer to ggplot2 manual for the possible values;
#' by default it centers the plot title plot.title = element_text(hjust = 0.5)
#' @param plot.subtitle - can be any value that is accepted in the argument _plot.subtitle_ in the _theme_ function
#' of ggplot2 package,refer to ggplot2 manual for the possible values;
#' by default it centers the plot subtitle plot.title = element_text(hjust = 0.5)
#' @param ... -  other arguments of the _theme_ function, refer to ggplot2 manual
#' @return ggplot2 object
#' @import ggplot2
#' @import ggpubr
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
#' LoTTA_plot_outcome(out,n_eval=100)

LoTTA_plot_outcome<-function(LoTTA_posterior,nbins=100,probs=c(0.025,0.5,0.975),n_eval=200,col_line='#E69F00',size_line=0.1,alpha_interval=0.35,col_dots='gray',size_dots=3,alpha_dots=0.6,col_cutoff='black',
                             title="Outcome function",subtitle=NULL,y_lab="",x_lab=expression(paste(italic('X')," - score")) ,plot.theme=theme_classic(base_size = 14),
                             legend.position='none', name_legend='Legend', labels_legend = 'median outcome fun.',
                             text=element_text(family='serif'),legend.text=element_text(size = 14),
                             plot.title = element_text(hjust = 0.5),plot.subtitle=element_text(hjust = 0.5),...){
  if(length(LoTTA_posterior$Priors$`Outcome prior`)==2){
    plot_outcome_CONT(LoTTA_posterior,nbins,probs,n_eval,col_line,size_line,alpha_interval,col_dots,size_dots,alpha_dots,col_cutoff,
                      title,subtitle,y_lab,x_lab ,plot.theme,legend.position, name_legend, labels_legend,text,legend.text,plot.title,plot.subtitle,...)

  }
  else{
    plot_outcome_BIN(LoTTA_posterior,nbins,probs,n_eval,col_line,size_line,alpha_interval,col_dots,size_dots,alpha_dots,col_cutoff,
                     title,subtitle,y_lab,x_lab ,plot.theme,legend.position, name_legend, labels_legend,text,legend.text,plot.title,plot.subtitle,...)

  }


}


#' Function that plots the median (or another quantile) of the LoTTA posterior treatment probability function along with the quanatile-based credible interval.
#' The function is plotted on top of the binned input data. To bin the data, the score data
#' is divided into bins of fixed length, then the proportion of treated is calculated in each bin.
#' The proportions are plotted against the average values of the score in the corresponding bins.
#' The data is binned separately on each side of the cutoff, the cutoff is marked on the plot
#' with a dotted line. In case of an unknown cutoff, the MAP estimate is used.
#' @param LoTTA_posterior - output of one of the LoTTA functions (LoTTA_fuzzy_CONT, LoTTA_fuzzy_BIN, LoTTA_treatment)
#' with all the parameters sampled (the default option in those functions)
#' @param nbins - number of bins to aggregate the input data
#' @param probs - list of three quantiles, the first and the last one define the quanatile-based credible interval,
#' the middle value defines the quantile of the posterior function to plot;
#' by default the quantiles correspond to the median posterior function and 95% credible interval
#' probs=c(0.025,0.5,0.975)
#' @param n_eval - n_eval*range(X) is the number of points at which each posterior function is evaluated,
#' the higher number means slower computing time and a smoother plot; default n_eval=200
#' @param col_line - the color of the line and the band
#' @param size_line - thickness of the line
#' @param alpha_interval - alpha value of the band, lower values correspond to a more transparent color
#' @param col_dots - color of the dots that correspond to the binned data
#' @param size_dots - size of the dots that correspond to the binned data
#' @param alpha_dots - transparency of the dots that correspond to the binned data,
#'  lower values correspond to a more transparent color
#' @param col_cutoff - color of the dotted line at the cutoff
#' @param title - title of the plot
#' @param subtitle - subtitle of the plot
#' @param y_lab - label of the y-axis
#' @param x_lab - label of the x-axis
#' @param plot.theme - ggplot2 plot theme (see https://ggplot2.tidyverse.org/reference/ggtheme.html)
#' possibly with additional arguments, it takes the default value plot.theme=theme_classic(base_size = 14),
#' @param legend.position - position of the legend, refer to ggplot2 manual for the possible values;
#' by default legend is not printed legend.position='none'
#' @param name_legend - title of the legend
#' @param labels_legend - the label of the plotted function in the legend
#' @param text - can be any value that is accepted in the argument _text_ in the _theme_ function
#' of ggplot2 package,refer to ggplot2 manual for the possible values;
#' by default is changes font to a serif one text=element_text(family='serif')
#' @param legend.text - can be any value that is accepted in the argument _legend.text_ in the _theme_ function
#' of ggplot2 package,refer to ggplot2 manual for the possible values;
#' by default is changes font size to the legend to 14 legend.text=element_text(size = 14)
#' @param plot.title - can be any value that is accepted in the argument _plot.title_ in the _theme_ function
#' of ggplot2 package,refer to ggplot2 manual for the possible values;
#' by default it centers the plot title plot.title = element_text(hjust = 0.5)
#' @param plot.subtitle - can be any value that is accepted in the argument _plot.subtitle_ in the _theme_ function
#' of ggplot2 package,refer to ggplot2 manual for the possible values;
#' by default it centers the plot subtitle plot.title = element_text(hjust = 0.5)
#' @param ... -  other arguments of the _theme_ function, refer to ggplot2 manual
#' @return ggplot2 object
#' @import ggplot2
#' @import ggpubr
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
#' LoTTA_plot_treatment(out,nbins =15,n_eval = 100)

LoTTA_plot_treatment<-function(LoTTA_posterior,nbins=100,probs=c(0.025,0.5,0.975),n_eval=200,col_line='#E69F00',size_line=0.1,alpha_interval=0.35,col_dots='gray',size_dots=3,alpha_dots=0.6,col_cutoff='black',
                               title="Treatment probability function",subtitle=NULL,y_lab="",x_lab=expression(paste(italic('X')," - score")) ,plot.theme=theme_classic(base_size = 14),
                               legend.position='none', name_legend='Legend', labels_legend = 'median treatment prob. fun.',
                               text=element_text(family='serif'),legend.text=element_text(size = 14),
                               plot.title = element_text(hjust = 0.5),plot.subtitle=element_text(hjust = 0.5),...){
  LoTTA_jags=LoTTA_posterior$JAGS_output
  Samples_FULL=data.frame(data.frame(combine.mcmc(LoTTA_jags)))
  d_x=LoTTA_posterior$Normalized_data$d_X
  s_x=LoTTA_posterior$Normalized_data$s_X
  X=LoTTA_posterior$Normalized_data$'X'*s_x+d_x
  T=LoTTA_posterior$Normalized_data$'T'
  if(is.null(LoTTA_posterior$Effect_estimate$`MAP cutoff`)==TRUE){
    c=LoTTA_posterior$Normalized_data$c
    Samples_FULL[,'c']=rep(c,nrow(Samples_FULL))
  }
  else{
    c=LoTTA_posterior$Effect_estimate$`MAP cutoff`
  }
  rangeX=max(X)-min(X)
  binsize=rangeX/nbins
  binned_data=Bin_data(T,X,c,binsize)
  X_plot=binned_data[[2]]
  T_plot=binned_data[[1]]

  x=seq(min(X),max(X),length.out = n_eval*rangeX)
  fun_sampleF=apply(Samples_FULL,1,treatment_function_sample,x=x,d_x=d_x,s_x=s_x)
  q_funF=apply(fun_sampleF,1,quantile,probs=probs)
  plot.df=data.frame('X'=X_plot,'Y'=T_plot)
  lower95=q_funF[1,]
  upper95=q_funF[3,]
  Y=q_funF[2,]
  lineF.df=data.frame(X=x,Y=Y,lower95=lower95,upper95=upper95)
  theme_def=plot.theme
  ggscatter(plot.df, x = "X", y = "Y",
            color = col_dots,
            size = size_dots, alpha = alpha_dots)+geom_vline(xintercept = c,linetype='dotted',color=col_cutoff)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    geom_ribbon(aes(x=X,y=Y,ymin = lower95, ymax = upper95,fill='Full'),data=lineF.df,alpha=alpha_interval)+geom_point(aes(x=X,y=Y),data=lineF.df,size=size_line,col=col_line)+
    labs( y=y_lab, x=x_lab,title=title,subtitle = subtitle)+theme_def+
    theme(text = text,legend.text =legend.text,plot.title = plot.title,plot.subtitle = plot.subtitle,legend.position = legend.position,...)+ scale_fill_manual(name=name_legend,values=c('Full'=col_line),labels = labels_legend)

}

#' Function that visualizes the impact of the cutoff location on the treatment effect estimate.
#' It plots too figures. The bottom figure depicts the posterior density of the cutoff location.
#' The top figure depicts the box plot of the treatment effect given the cutoff point.
#' If the prior on the cutoff location was discrete each box corresponds to a distinct cutoff point.
#' If the prior was continuous each box correspond to an interval of cutoff values
#' (the number of intervals can be changed through nbins).
#' @param LoTTA_posterior - output of one of the LoTTA functions (LoTTA_fuzzy_CONT, LoTTA_fuzzy_BIN)
#' with all parameters sampled (the default option in those functions)
#' @param nbins - number of bins to aggregate the input data
#' @param probs - list of two quantiles that limit the range of cutoff values displayed on the plots
#' @param x_lab - label of the x-axis
#' @param y_lab1 - label of the y-axis of the bottom plot
#' @param y_lab2 - label of the y-axis of the top plot
#' @param title - title of the plot
#' @param axis.text - can be any value that is accepted in the argument _axis.text_ in the _theme_ function
#' of ggplot2 package,refer to ggplot2 manual for the possible values;
#' by default is changes font to a serif one axis.text=element_text(family = "sans",size = 10)
#' @param text - can be any value that is accepted in the argument _text_ in the _theme_ function
#' of ggplot2 package,refer to ggplot2 manual for the possible values;
#' by default is changes font to a serif one text=element_text(family='serif')
#' @param plot.theme - ggplot2 plot theme (see https://ggplot2.tidyverse.org/reference/ggtheme.html)
#' possibly with additional arguments, it takes the default value plot.theme=theme_classic(base_size = 14),
#' @param text - can be any value that is accepted in the argument _text_ in the _theme_ function
#' of ggplot2 package,refer to ggplot2 manual for the possible values;
#' by default is changes font to a serif one text=element_text(family='serif')
#' @param plot.title - can be any value that is accepted in the argument _plot.title_ in the _theme_ function
#' of ggplot2 package,refer to ggplot2 manual for the possible values;
#' by default it centers the plot title plot.title = element_text(hjust = 0.5)
#' @param ... -  other arguments of the _theme_ function, refer to ggplot2 manual
#' @return ggplot2 object
#' @import ggplot2
#' @import ggpubr

LoTTA_plot_effect_CONT<-function(LoTTA_posterior,nbins=10,probs = c(0.025,0.975),x_lab="Cutoff location",y_lab1="Treatment effect",y_lab2="Density of cutoff location",title="Cutoff location vs. Treatment effect",
                                 axis.text=element_text(family = "sans",size = 10),text=element_text(family='serif'),plot.theme=theme_classic(base_size = 14),plot.title = element_text(hjust = 0.5),...){
  c = as.numeric(LoTTA_posterior$Samples$c)
  eff = as.numeric(LoTTA_posterior$Samples$eff)
  dat = data.frame(c=c,eff=eff)
  bounds_c=quantile(c,probs = probs)
  group_c=seq(bounds_c[[1]],bounds_c[[2]],(bounds_c[[2]]-bounds_c[[1]])/nbins)
  lab=signif((group_c[2:length(group_c)]+group_c[1:length(group_c)-1])*0.5,2)
  dat$group_c=as.numeric(cut(c, group_c,right=TRUE,labels = lab))
  dat=na.omit(dat)
  dat$group_c=lab[dat$group_c]
  pl1=ggplot(data = dat, mapping = aes(x = group_c))+
    geom_boxplot(data=dat,aes(y = eff, group=group_c))+labs( y=y_lab1, title=title,x="",fill='posterior prob. of c') +
    plot.theme+theme(text = text,axis.text=axis.text,plot.title = plot.title,axis.ticks.x = element_blank(),...)
  pl2=ggplot(data = dat, mapping = aes(x = c)) +
    geom_density()+labs( y=y_lab2, x=x_lab) +
    plot.theme+theme(text = text,axis.text=axis.text,plot.title = plot.title,axis.ticks.x = element_blank(),...)
  ggarrange(pl1,pl2,nrow = 2,align = 'hv')
}
#' Function that visualizes the impact of the cutoff location on the treatment effect estimate.
#' It plots too figures. The bottom figure depicts the posterior density of the cutoff location.
#' The top figure depicts the box plot of the treatment effect given the cutoff point.
#' If the prior on the cutoff location was discrete each box corresponds to a distinct cutoff point.
#' If the prior was continuous each box correspond to an interval of cutoff values
#' (the number of intervals can be changed through nbins).
#' @param LoTTA_posterior - output of one of the LoTTA functions (LoTTA_fuzzy_CONT, LoTTA_fuzzy_BIN)
#' with all parameters sampled (the default option in those functions)
#' @param probs - list of two quantiles that limit the range of cutoff values displayed on the plots
#' @param x_lab - label of the x-axis
#' @param y_lab1 - label of the y-axis of the bottom plot
#' @param y_lab2 - label of the y-axis of the top plot
#' @param title - title of the plot
#' @param axis.text - can be any value that is accepted in the argument _axis.text_ in the _theme_ function
#' of ggplot2 package,refer to ggplot2 manual for the possible values;
#' by default is changes font to a serif one axis.text=element_text(family = "sans",size = 10)
#' @param text - can be any value that is accepted in the argument _text_ in the _theme_ function
#' of ggplot2 package,refer to ggplot2 manual for the possible values;
#' by default is changes font to a serif one text=element_text(family='serif')
#' @param plot.theme - ggplot2 plot theme (see https://ggplot2.tidyverse.org/reference/ggtheme.html)
#' possibly with additional arguments, it takes the default value plot.theme=theme_classic(base_size = 14),
#' @param text - can be any value that is accepted in the argument _text_ in the _theme_ function
#' of ggplot2 package,refer to ggplot2 manual for the possible values;
#' by default is changes font to a serif one text=element_text(family='serif')
#' @param plot.title - can be any value that is accepted in the argument _plot.title_ in the _theme_ function
#' of ggplot2 package,refer to ggplot2 manual for the possible values;
#' by default it centers the plot title plot.title = element_text(hjust = 0.5)
#' @param ... -  other arguments of the _theme_ function, refer to ggplot2 manual
#' @return ggplot2 object
#' @import ggplot2
#' @import ggpubr
LoTTA_plot_effect_DIS<-function(LoTTA_posterior,probs = c(0.025,0.975),x_lab="Cutoff location",y_lab1="Treatment effect",y_lab2="Density of cutoff location",title="Cutoff location vs. Treatment effect",
                                axis.text=element_text(family = "sans",size = 10),text=element_text(family='serif'),plot.theme=theme_classic(base_size = 14),plot.title = element_text(hjust = 0.5),...){
  c = as.numeric(LoTTA_posterior$Samples$c)
  eff = as.numeric(LoTTA_posterior$Samples$eff)
  bounds_c=quantile(c,probs = probs)
  dat = data.frame(c=c,eff=eff)
  dat=subset(dat,c<=bounds_c[[2]]&c>=bounds_c[[1]])
  pl1=ggplot(data = dat, mapping = aes(x = c))+
    geom_boxplot(data=dat,aes(y = eff, group=c))+labs( y=y_lab1, title=title,x="",fill='posterior prob. of c') +
    plot.theme+theme(text = text,axis.text=axis.text,plot.title = plot.title,axis.ticks.x = element_blank(),...)
  pl2=ggplot(data = dat, mapping = aes(x = c)) +
    geom_density()+labs( y=y_lab2, x=x_lab) +
    plot.theme+theme(text = text,axis.text=axis.text,plot.title = plot.title,axis.ticks.x = element_blank(),...)
  ggarrange(pl1,pl2,nrow = 2,align = 'hv')
}

#' LoTTA_plot_effect
#'
#' Function that visualizes the impact of the cutoff location on the treatment effect estimate.
#' It plots too figures. The bottom figure depicts the posterior density of the cutoff location.
#' The top figure depicts the box plot of the treatment effect given the cutoff point.
#' If the prior on the cutoff location was discrete each box corresponds to a distinct cutoff point.
#' If the prior was continuous each box correspond to an interval of cutoff values
#' (the number of intervals can be changed through nbins).
#' @param LoTTA_posterior - output of one of the LoTTA functions (LoTTA_fuzzy_CONT, LoTTA_fuzzy_BIN)
#' with all parameters sampled (the default option in those functions)
#' @param nbins - number of bins to aggregate the input data
#' @param probs - list of two quantiles that limit the range of cutoff values displayed on the plots
#' @param x_lab - label of the x-axis
#' @param y_lab1 - label of the y-axis of the bottom plot
#' @param y_lab2 - label of the y-axis of the top plot
#' @param title - title of the plot
#' @param axis.text - can be any value that is accepted in the argument _axis.text_ in the _theme_ function
#' of ggplot2 package,refer to ggplot2 manual for the possible values;
#' by default is changes font to a serif one axis.text=element_text(family = "sans",size = 10)
#' @param text - can be any value that is accepted in the argument _text_ in the _theme_ function
#' of ggplot2 package,refer to ggplot2 manual for the possible values;
#' by default is changes font to a serif one text=element_text(family='serif')
#' @param plot.theme - ggplot2 plot theme (see https://ggplot2.tidyverse.org/reference/ggtheme.html)
#' possibly with additional arguments, it takes the default value plot.theme=theme_classic(base_size = 14),
#' @param text - can be any value that is accepted in the argument _text_ in the _theme_ function
#' of ggplot2 package,refer to ggplot2 manual for the possible values;
#' by default is changes font to a serif one text=element_text(family='serif')
#' @param plot.title - can be any value that is accepted in the argument _plot.title_ in the _theme_ function
#' of ggplot2 package,refer to ggplot2 manual for the possible values;
#' by default it centers the plot title plot.title = element_text(hjust = 0.5)
#' @param ... -  other arguments of the _theme_ function, refer to ggplot2 manual
#' @return ggplot2 object
#' @import ggplot2
#' @import ggpubr
#' @export
#' @examples
#' \donttest{ # code to generate the data
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
#' funC <- function(X) {
#'   Y = X
#'   Y = 1 / (1 + exp(-2 * X)) - 0.5 + 0.4
#'   return(Y)
#' }
#'
#' funC_sample <- function(X) {
#'   Y = funC(X)+ rnorm(length(X), 0, 0.025)
#'   return(Y)
#' }
#'
#' # data generation
#' N=150 # try different dataset size
#' set.seed(1234)
#' X = sort(runif(N, -1, 1))
#' T = sample_prob55(X)
#' Y = funC_sample(X)
#' c_prior=list(clb=-0.25,cub=0.25)
#'
#' # running LoTTA function on sharp RDD with continuous outcomes;
#' # cutoff = 0, treatment effect = 0
#' # remember to check convergence and adjust burnin, sample and adapt if needed
#' out = LoTTA_fuzzy_CONT(X,T,Y,c_prior,burnin = 6000,sample = 2000,adapt=1000,seed=1234)
#' LoTTA_plot_effect(out,nbins = 15)}

LoTTA_plot_effect<-function(LoTTA_posterior,nbins=10,probs=c(0.025,0.975),x_lab="Cutoff location",y_lab1="Treatment effect",y_lab2="Density of cutoff location",title="Cutoff location vs. Treatment effect",
                            axis.text=element_text(family = "sans",size = 10),text=element_text(family='serif'),plot.theme=theme_classic(base_size = 14),plot.title = element_text(hjust = 0.5),...){
  if(is.null(LoTTA_posterior$Priors$`Cutoff prior`$cstart)==FALSE){
    pl1=LoTTA_plot_effect_DIS(LoTTA_posterior,probs,x_lab,y_lab1,y_lab2,title,
                              axis.text,text,plot.theme,plot.title,...)
    print(pl1)
  }

  if(is.null(LoTTA_posterior$Priors$`Cutoff prior`$clb)==FALSE){
    pl1=LoTTA_plot_effect_CONT(LoTTA_posterior,nbins,probs,x_lab,y_lab1,y_lab2,title,
                               axis.text,text,plot.theme,plot.title,...)
    print(pl1)
  }
  return(pl1)

}
