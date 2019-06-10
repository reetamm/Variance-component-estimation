#library(gplots)
#Settings:
#k=10,20
#N=10*k 50*k 100*k
#Imbalance= balanced, slight unbalanced, highly unbalanced



mdata <- function(k=10, n = 10, sig2 = 0.5, sigt2 = 1, mu = 10, mut = numeric(k)){
    y <- numeric(k*n);
    grp <- numeric(k*n);
    out <- data.frame(y,grp);

    a <- 0

    for (i in 1:k) {
        t <- rnorm(1, mut[i], sqrt(sigt2))
        for (j in 1:n) {
            out[j + a, 1] <- mu + t + rnorm(1, 0, sqrt(sig2));
            out[j + a, 2] <- i;
        }
        a <- a + n;
    }
    return(out);
}


mdata2 <- function(k = 10, n = 10, b = 1, sig2 = 0.5, m = 5, sigt2 = 1, mu = 10, mut = numeric(k)) {
    if (b == 1) {
        data <- mdata(k = k, n = n, sig2 = sig2, sigt2 = sigt2, mu = mu, mut = mut) #changed the name of the variable from out to data
    }

    if (b == 2) {
        n2 <- round((k*n - 2*m)/k)
        out1 <- mdata(k = k, n = n2, sig2 = sig2, sigt2 = sigt2, mu = mu, mut = mut)
        out2 <- mdata(k = 2, n = m, sig2 = sig2, sigt2 = sigt2, mu = mu, mut = mut)
        data <- rbind(out1, out2) #changed the name of the variable from out to data
    }

    if (b == 3) {
        n2 <- round((n - m)*k/2)
        out1 <- mdata(k = 2, n = n2, sig2 = sig2, sigt2 = sigt2, mu = mu, mut = mut)
        out2 <- mdata(k = k, n = m, sig2 = sig2, sigt2 = sigt2, mu = mu, mut = mut)
        data <- rbind(out1, out2) #changed the name of the variable from out to data
    }
    data$grp <- as.factor(data$grp)
    #plot=plotmeans(y ~ grp, data = data);
    #output=list(data=data,plot=plot)

    return(data)
}
#for mdata3, n_i is a vector containing the number of samples per group, k is number of groups
#by defaut, mdata3 generates balanced data with 10 groups and 10 samples per group
mdata3 <- function(n_i = NULL, k = 10, n = 10, b = 1, m = 5, sig2 = 0.5, sigt2 = 1, mu = 10, mut = numeric(k)) {
    if (sum(is.na(n_i)) == 0 & sum(n_i > 0) == length(n_i) & length(n_i) > 0) #as long as n_i is a valid vector of group sizes
    {
        k <- length(n_i) #extracts number of groups
        muv <- rnorm(k, mean = 0, sd = sqrt(sigt2)) #generates vector of means to be used for each group
        y <- 0
        grp <- 0
        out <- cbind(y, grp) #container of data
        for (i in 1:k) {
            yd <- mu + muv[i] + rnorm(n_i[i], mean = 0, sd = sqrt(sig2)) #data for group i
            out <- rbind(out, cbind(yd, numeric(n_i[i]) + i)) #index i
        }
        out <- out[-1, ] #removes first row
        out <- as.data.frame(out) #turns matrix to dataframe
        out$grp <- as.factor(as.character(out$grp)) #coerces grp to be treated as factor
    }
    else
    {
        out <- mdata2(k = k, n = n, b = b, sig2 = sig2, sigt2 = sigt2, m = 5, mu = 10, mut = numeric(k))
    }
    #plotmeans(y ~ grp, data = out); #plot
    return(out)
}

# Function 1. extractparam()
# Extract the paramters for the model. Uses Seed = 1
# Returns a vector of length 3 with mean, sigma2 and sigmat2
extractparam <- function(data, type = 'nonparametric')
  {
  #insert error checks here to verify first variable is numeric, 2nd is group. drop any other variables
  names(data) <- c("y", "group")            # Change mdata2() nomenclature from grp to group
  size <- as.vector(table(data$group))     # Check if each group has at least 5 obs
  if (min(size) < 5)
    print("Warning: at least one group has less than 5 observations")


  unique_groups <- as.vector(unique(data$group)) #extract the group names
  n_i <- NULL
  for (i in unique_groups)               # Run loop for each unique group. Vectorize within
  {
    n_i[i] <- sum(data$group == i)
  }
  n_i <- as.vector(n_i)

  if (type == 'nonparametric' | type == 'np' | type == 'n')

  {
    param_np <- list(unique_groups = unique_groups, n_i = n_i)
    return(param_np)
  }

  if (type == 'parametric' | type == 'p')
  {
  fit <- lmer(y ~ (1 | group), data = data)
  mu_hat <- as.data.frame(fixef(fit))[1,1] # Common mean
  v <- as.data.frame(VarCorr(fit))
  sigmat2_hat <- v[1,4]                    # Group Variance - change 4 to 5 to get SD instead
  sigma2_hat <- v[2,4]                     # Noise Variance - change 4 to 5 to get SD instead
  n_i_sq <- sum(n_i^2)
  param_p <- c(mu_hat,sigma2_hat,sigmat2_hat,n_i_sq)
  return(param_p)
  }

  if (type == 'blup')
  {

    fitted <- lmer(y~ (1 | group), data = data)
    mu_hat <- as.data.frame(fixef(fitted))[1,1]
    tau_hat <- as.vector(ranef(fitted)$group[ ,1])
    resid <- as.vector(resid(fitted))

    param_blup <- list(unique_groups = unique_groups, n_i = n_i, mu_hat = mu_hat, tau_hat = tau_hat, resid = resid)
    return(param_blup)
  }

  if (type == 'blue')
  {
    fitted <- lm(y ~ group - 1, data = data)
    mu_hat <- mean(data$y)
    tau_hat <- as.vector(fitted$coefficients - mu_hat)
    resid <- as.vector(fitted$residuals)

    param_blue <- list(unique_groups = unique_groups, n_i = n_i, mu_hat = mu_hat, tau_hat = tau_hat, resid = resid)
    return(param_blue)
  }
}
