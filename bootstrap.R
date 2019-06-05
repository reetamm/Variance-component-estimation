#Put this file and the data generator in the same folder and set it as working directory
#setwd("D:/Google Drive/UMBC/2019 02 Spring/612/project")
#setwd("C:/Users/reetam/Google Drive/2019 02 Spring/612/project")
#install.packages("gridExtra")

library(lme4)
#library(ggplot2)
#library(gridExtra)
source("data_generator.R") #imports the contents of the data generator

###########################################################
######## Bootstrap Functions ##############################
###########################################################
####### 1. extractparam() is a paramter extraction function
####### 2. npbootsample() generates single nonparametric bootstrap sample
####### 3. pbootsample() generates sinlge parametric bootstrap sample
####### 4. blupbootsample() does 2 runs; run 1 extracts parameters as list. run 2 uses that list to generate blup sample.
####### 5. bootrem() generates M bootstrap estimates of mu, sigma2, sigmat2
####### 6. bootconf() gets the confidence intervals based on bootrem()
####### 7. bootperf() is a wrapper for bootconf and gives coverage and performance

# Function 1. extractparam()
# Extract the paramters for the model. Uses Seed = 1
# Returns a vector of length 3 with mean, sigma2 and sigmat2
extractparam = function(data,type='nonparametric')
{
  #insert error checks here to verify first variable is numeric, 2nd is group. drop any other variables
  names(data)=c("y","group")            # Change mdata2() nomenclature from grp to group
  size=as.vector(table(data$group))     # Check if each group has at least 5 obs           
  if(min(size)<5)
    print("Warning: at least one group has less than 5 observations")
  
  
  unique_groups = as.vector(unique(data$group)) #extract the group names
  n_i=NULL
  for(i in unique_groups)               # Run loop for each unique group. Vectorize within
  {
    n_i[i]=sum(data$group==i)
  }
  n_i = as.vector(n_i)
  
  if(type=='nonparametric' | type== 'np' | type == 'n')
    
  {
    param_np = list(unique_groups = unique_groups,n_i = n_i)
    return(param_np)
  }
  
  if(type == 'parametric' | type == 'p')
  {
  fit=lmer(y~ (1 | group), data=data)   
  mu_hat = as.data.frame(fixef(fit))[1,1] # Common mean
  v = as.data.frame(VarCorr(fit))       
  sigmat2_hat = v[1,4]                    # Group Variance - change 4 to 5 to get SD instead
  sigma2_hat = v[2,4]                     # Noise Variance - change 4 to 5 to get SD instead
  n_i_sq=sum(n_i^2)
  param_p = c(mu_hat,sigma2_hat,sigmat2_hat,n_i_sq)
  return(param_p)
  }
  
  if(type == 'blup')
  {
    
    fitted=lmer(y~ (1 | group), data=data)
    mu_hat=as.data.frame(fixef(fitted))[1,1]
    tau_hat = as.vector(ranef(fitted)$group[,1])
    resid = as.vector(resid(fitted))
    
    param_blup = list(unique_groups=unique_groups,n_i=n_i,mu_hat=mu_hat,tau_hat=tau_hat,resid=resid)
    return(param_blup)          
  }
  
  if(type == 'blue')
  {
    fitted = lm(y~group -1, data = data)
    mu_hat = mean(data$y)
    tau_hat = as.vector(fitted$coefficients - mu_hat)
    resid = as.vector(fitted$residuals)
    
    param_blue = list(unique_groups=unique_groups,n_i=n_i,mu_hat=mu_hat,tau_hat=tau_hat,resid=resid)
    return(param_blue)  
  }
}


# Function 2. npbootsample()
# non parametric bootstrap sample generation
# Returns a resampled dataframe with variables y and group
npbootsample = function(data)           # input dataframe with 1st variable numeric, 2nd variable factor.
{
  params = extractparam(data, type = 'nonparametric')
  unique_groups = params$unique_groups #extract the group names
  n_i = params$n_i
  ynew=NA
  groupnew=NA
  for(i in 1:length(n_i))               # Run loop for each unique group. Vectorize within
  {
    temp = sample(data$y[data$grp==unique_groups[i]],n_i[i],replace = T) # Resample from i-th group with replacement
    ynew = c(ynew,temp)                   # Append group i's resample 
    groupnew = c(groupnew,rep(i,n_i[i]))       # n_i observations of group i
  }
  newdata = data.frame(y=ynew,group=groupnew)     # Set data frame
  newdata = newdata[-1,]                # Drop the first NA row
  return(newdata)
}


# Function 3. pbootsample()
# parametric bootstrap sample generation
# Returns a resampled dataframe with variables y and group
pbootsample = function(data,k,n,b,random=F) # k, n, b inherited from the top
{
  params = extractparam(data, type = 'parametric')
  mu_hat=params[1];sigma2_hat=params[2];sigmat2_hat=params[3];n_i_sq=params[4]
  df1=n*k-k
  sigma2_star = sigma2_hat*rchisq(1, df1)/(df1)
  sigmat2_star=sigmat2_hat
  if(random==T)
  {
    df2=k-1
    constant = (n*k-n_i_sq/(n*k))/(k-1)
    sigmat2_star = (sigma2_hat+constant*sigmat2_hat)*rchisq(1,df2)/(df2)
  }
  
  data_star = mdata2(k,n,b,sig2=sigma2_star, sigt2=sigmat2_star,mu=mu_hat,mut=numeric(k))  # m inherited from top
  names(data_star) = c("y","group")     # Rename grp to group
  return(data_star)  
}

# Function 4. blupbootsample()
# BLUP Bootstrap sample generation. Can be run in 2 modes
# Mode 1: extract=T, data needed - will extract the parameters as a list
# Mode 2: extract=F, data not needed, list from mode 1 as input. outputs a single resample.
blupbootsample = function(data, effect.type = 'random')
{
  if(effect.type == 'random' | effect.type == 'r')
    params = extractparam(data, type = 'blup')
  if(effect.type == 'fixed' | effect.type == 'f')
    params = extractparam(data, type = 'blue')
  ystar=NULL
  groupnew=NULL
  for(i in 1:length(params$unique_groups))
  {
    tau_star=sample(params$tau_hat,1) # Sample a random tau from the vector of estimated tau
    resid_star = sample(params$resid,params$n_i[i]) # random innovations. How many? - equal to group size
    temp = params$mu_hat+tau_star+resid_star  # get y_ij for i-th group
    ystar = c(ystar,temp)                         # append y_ij
    groupnew = c(groupnew, rep(params$unique_groups[i],params$n_i[i]))
  }
  data_star = data.frame(y=ystar,group=groupnew)
  return(data_star)
}

# Function 5. bootrem()
# for getting bootstrap estimates. type = parametric (or p), nonparametric (or n), blup (or b)
bootrem = function(data = data, M=1000,type='n',verbose=F,effect.type='random')       # M is number of bootstrap samples, removed seed and added data as function input
{
  mu=numeric(M);sigma2=numeric(M);sigmat2=numeric(M)            # Initialize with blanks
  
  for(i in 1:M)
  {
    if(verbose==T)
      print(i);				    # iteration counter
    if(type=='nonparametric' | type=='n' | type == 'np')                 # Non parametric bootstrap
    {
      sample=npbootsample(data)         # Get the i-th nonparametric sample
    }                                   # OR
    if(type=='parametric'|type=='p')    # Parametric bootstrap
    {
      sample=pbootsample(data,k,n,b)
    }                                   # OR
    if(type=='blup'|type=='b')          # BLUP bootstrap
    {
      sample = blupbootsample(data,effect.type)  #data not needed. extract=F. Needs list of params from initial run
    }
    
    fit=lmer(y~ (1 | group),data=sample) # Run the ANOVA
    mu_1 = as.data.frame(fixef(fit))[1,1] 
    v = as.data.frame(VarCorr(fit))
    sigmat2_1 = v[1,4]
    sigma2_1 = v[2,4]
    mu[i] = mu_1
    sigmat2[i] = sigmat2_1
    sigma2[i] = sigma2_1
  }
  estimates = data.frame(mu=mu,sigma2=sigma2,sigmat2=sigmat2) 
  return(estimates)                     # Dataframe with 3 columns of M estimates
}


# Function 6. bootconf()
# Input the dataframe from bootrem(), vector of original estimates (potentially
# from extractparam()) and alpha (default 5%)
# two.sided gives 2 sided intervals for variance components
# sd=T gives square root of variance
bootconf = function(estimates, alpha=.05, two.sided=F, sd = F,plots=T) 
{
  mu_L = quantile(estimates$mu,alpha/2)
  mu_U = quantile(estimates$mu, 1 - alpha/2)
  if(two.sided==T)
  {
    sigma2_L = quantile(estimates$sigma2, alpha/2)
    sigma2_U = quantile(estimates$sigma2, 1 - alpha/2)
    sigmat2_L = quantile(estimates$sigmat2, alpha/2)
    sigmat2_U = quantile(estimates$sigmat2, 1 - alpha/2)
  }
  if(two.sided==F)
  {
    sigma2_L = quantile(estimates$sigma2, alpha)
    sigma2_U = Inf
    sigmat2_L = quantile(estimates$sigmat2, alpha)
    sigmat2_U = Inf
  }
  if(sd==F)
  {
    CI = matrix(c(mu_L,mu_U,sigma2_L,sigma2_U,sigmat2_L,sigmat2_U),ncol=2,byrow = T)
    colnames(CI) = c(as.character(alpha/2),as.character(1-alpha/2))
    rownames(CI) = c("(Intercept)",".sigma2",".sigmat2")
  }
  if(sd==T)
  {
    CI = matrix(c(mu_L,mu_U,sqrt(sigma2_L),sqrt(sigma2_U),sqrt(sigmat2_L),sqrt(sigmat2_U)),ncol=2,byrow = T)
    colnames(CI) = c(as.character(alpha/2),as.character(1-alpha/2))
    rownames(CI) = c("(Intercept)",".sigma",".sigmat")
  }
  if(plots==T)
  {
    p1 = ggplot(estimates,aes(x=sigma2))+geom_density()
    p2 = ggplot(estimates,aes(x=sigmat2))+geom_density()
    grid.arrange(p1,p2,nrow=2)
  }
  return(CI)
}

# Function 7: bootperf()
# Gives coverage and average length for sigma and sigma_tau based on 'max' bootstrap samples, each of size 'M'
# Can set whether we want SD or Var with 'sd' and 1 or 2 sided with 'two.sided'. Plots set to F, alpha at .05

bootperf = function(k,b,var,M,max,two.sided=F,sd=F,alpha=0.05,run.parametric=T,run.nonparametric=T,run.blup=T,verbose=T, effect.type = 'random')
{
  para=numeric(max);paracheck=numeric(max)
  nonpara=numeric(max);nonparacheck=numeric(max)
  blup=numeric(max);blupcheck=numeric(max)
  paratau=numeric(max);parataucheck=numeric(max)
  nonparatau=numeric(max);nonparataucheck=numeric(max)
  bluptau=numeric(max);bluptaucheck=numeric(max)
  if(sd==T)
    pow = 0.5
  else
    pow = 1
  i=1
  while(i <= max)
  {
    print(i)
    data=mdata2(k,10,b,sig2=0.5, sigt2=var,mu=10)
    fit=lmer(y~ (1 | grp), data=data)
    if(isSingular(fit)==F)
    {
    
    if(run.nonparametric==T)
    {
      estimates1=bootrem(data,M=M,type = 'n',verbose = F)
      CIa = bootconf(estimates1,alpha,two.sided,sd, plots = F)
      
      if(CIa[2,1]<=0.5^pow)
      {
        nonparacheck[i] = T
        nonpara[i] = CIa[2,1]
      }
      if(CIa[3,1]<=var^pow)
      {
        nonparataucheck[i]=T
        nonparatau[i]=CIa[3,1]
      }
      if(verbose == T)
        print(i)
    }
    
    if(run.parametric==T)
    {
      estimates2=bootrem(data=data,M=M,type = 'p',verbose = F)
      CIb = bootconf(estimates2,alpha,two.sided,sd, plots = F)
      
      if(CIb[2,1]<=0.5^pow)
      {
        paracheck[i] = T
        para[i] = CIb[2,1]
      }
      if(CIb[3,1]<=var^pow)
      {
        parataucheck[i]=T
        paratau[i]=CIb[3,1]
      }
      if(verbose == T)
      print(i)
    }
    
    if(run.blup==T)
    {
      estimates3=bootrem(data=data,M=M,type = 'b',verbose = F)
      CIc = bootconf(estimates3,alpha,two.sided,sd, plots = F)
      
      if(CIc[2,1]<=0.5^pow)
      {
        blupcheck[i] = T
        blup[i] = CIc[2,1]
      }
      if(CIc[3,1]<=var^pow)
      {
        bluptaucheck[i]=T
        bluptau[i]=CIc[3,1]
      }
      if(verbose == T)
        print(i)
    }
      i = i+1
    }
  }
  if(run.nonparametric==T)
  {
    nonparacoverage = sum(nonparacheck)
    nonparalength = sum(nonpara)/nonparacoverage
    nonparataucoverage = sum(nonparataucheck)
    nonparataulength = sum(nonparatau)/nonparataucoverage
  }
  else
  {
    nonparacoverage = NA
    nonparalength = NA
    nonparataucoverage = NA
    nonparataulength = NA
  }
  
  if(run.parametric==T)
  {
    paracoverage = sum(paracheck)
    paralength = sum(para)/paracoverage
    parataucoverage = sum(parataucheck)
    parataulength = sum(paratau)/parataucoverage
  }
  else
  {
    paracoverage = NA
    paralength = NA
    parataucoverage = NA
    parataulength = NA
  }
  
  if(run.blup==T)
  {
    blupcoverage = sum(blupcheck)
    bluplength = sum(blup)/blupcoverage
    bluptaucoverage = sum(bluptaucheck)
    bluptaulength = sum(bluptau)/bluptaucoverage
  }
  else
  {
    blupcoverage = NA
    bluplength = NA
    bluptaucoverage = NA
    bluptaulength = NA
  }
  out = matrix(c(nonparacoverage,nonparalength,nonparataucoverage,nonparataulength,paracoverage,paralength,parataucoverage,parataulength,blupcoverage,bluplength,bluptaucoverage,bluptaulength),ncol=4,byrow = T)
  out=as.data.frame(out)
  
  if(effect.type == 'random' | effect.type == 'r')
    rownames(out) = c("nonparam","param","BLUP")
  if(effect.type == 'fixed' | effect.type == 'f')
    rownames(out) = c("nonparam","param","BLUP(BLUE)")
  
  if(sd==T)
    names(out) = c("sigma coverage","sigma length","sigmatau coverage","sigmatau length")
  else
    names(out) = c("sigmasq coverage","sigmasq length","sigmatausq coverage","sigmatausq length")
  return(out)
}



#############################################################################
#############################################################################
###################### MODEL RUNS BEGIN HERE ###############################

# Generate data
k = 10
n = 10
sigmat2 = .1
b_i = c(1,2,3)
#alldata = vector('list',3)
#set.seed(0)
#alldata[[1]]=mdata2(k,n,b=1,sig2=0.5, sigt2=sigmat2,mu=10)
#set.seed(1)
#alldata[[2]]=mdata2(k,n,b=2,sig2=0.5, sigt2=sigmat2,mu=10)
#set.seed(2)
#alldata[[3]]=mdata2(k,n,b=3,sig2=0.5, sigt2=sigmat2,mu=10)
#write.csv(alldata[[1]],'alldata1.csv',row.names = F)
#write.csv(alldata[[2]],'alldata2.csv',row.names = F)
#write.csv(alldata[[3]],'alldata3.csv',row.names = F)

#alldata[[1]] = read.csv('alldata1.csv')
#alldata[[2]] = read.csv('alldata2.csv')
#alldata[[3]] = read.csv('alldata3.csv')

#extractparam(data)
#fit=lmer(y~ (1 | grp), data=data)
#summary(fit)

for(b in b_i)
 # for(k in k_i)
  #  for(sigt2 in sigt2_i)
    {

#confint(fit)

out1=bootperf(k,b,sigmat2,M=1000, max=200, two.sided = F, sd=F, alpha = 0.05, run.parametric = T, run.nonparametric = T, run.blup = T,verbose = F)
out2=bootperf(k,b,sigmat2,M=1000, max=200, two.sided = F, sd=F, alpha = 0.05, run.parametric = F, run.nonparametric = F, run.blup = T,verbose = F, effect.type = 'fixed')
outname = paste0("k",k,"b",b,'tausq',sigmat2,'.csv')
write.csv(rbind(out1,out2[3,]),outname)
}
################################################################################
################################################################################


# bootrem() generates an Mx3 dataframe of mu, sigma2 and sigmat2.
# type = n, p, b for nonparametric, parametric and blup respectively
# verbose will show counter for every step
#estimates1=bootrem(data=data,M=1000,type = 'n',verbose = T)
#estimates2=bootrem(data=data,M=1000,type = 'p',verbose = T)
#estimates3=bootrem(data=data,M=1000,type = 'b',verbose = T)

# even if we set two.sided=F, the table header does not change - working on it.
# sd = F will output variance. 
#bootconf(estimates1,alpha = .05,two.sided = F,sd=T, plots = T)
#bootconf(estimates2,alpha = .05,two.sided = F,sd=T, plots = T)
#bootconf(estimates3,alpha = .05,two.sided = F,sd=T, plots = T)

