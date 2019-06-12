source("data_generator.R") #imports the contents of the data generator


genci <- function(data, s=10000, alpha = 0.05){
  genres <- aov(y ~ grp, data = data) #do 1-way anova on the data
  mse <- summary(genres)[[1]][[3]][2] #extract mse
  mst <- summary(genres)[[1]][[3]][1] #extract mst

  N <- nrow(data) #total sample size
  params <- extractparam(data, type = 'np')

  a <- length(params$unique_groups) #number of groups
  nv <- params$n_i #extract sample size for each group
  n0 <- (1/(a - 1))*(sum(nv) - sum(nv^2)/sum(nv)) #compute n0 which is needed for unbalanced case

  chi1 <- rchisq(s, df = a - 1)
  chi2 <- rchisq(s, df = N - a)
  T <- (mst*(a - 1)/n0)/chi1 - ((N - a)*mse/n0)/chi2 #samples for sigsqtau, derivation on paper
  t1 <- as.numeric(quantile(T, probs = alpha)) #1-tailed ci
  t2 <- as.numeric(quantile(T, probs = c(alpha/2, 1 - alpha/2))) #2-tailed ci
  t3 <- mean(T)
  return(list(t1,t2,t3))
}

test = genci(mdata3())
test
