#Put this file and the data generator in the same folder and set it as working directory
setwd("D:/Google Drive/UMBC/2019 02 Spring/612/project")
#setwd("C:/Users/reetam/Google Drive/2019 02 Spring/612/project")
#install.packages("gridExtra")
library(psych)
library(lme4)
#library(ggplot2)
#library(gridExtra)
source("data_generator.R") #imports the contents of the data generator


approxCI <- function(dataset, type='th', alpha = 0.05, two.sided = T)
  {
  params <- extractparam(dataset)
  groups <- params[[1]]
  n.i <- params[[2]]
  n <- sum(n.i)
  n.min <- min(n.i)
  n.max <- max(n.i)
  n.h <- harmonic.mean(params[[2]])
  k <- length(groups)

  fit <- lm(y ~ grp, data = dataset)
  ms2 <- summary(fit)$sigma

  ss3 <- 0
  y.bar <- numeric(k)
  for (i in 1:k)
    {
    y.bar[i] <- mean(dataset$y[dataset$grp == groups[i]])
    }
  ms3 <- var(y.bar)

  if (type == 'th')
    {
    if (two.sided == T)
      {
      th.lower <- ((k - 1)/qchisq((1 - alpha)/2, k - 1))*(ms3 - (ms2*qf((1 - alpha)/2, k - 1, n - k))/n.h)
      th.upper <- ((k - 1)/qchisq(alpha/2, k - 1))*(ms3 - (ms2*qf(alpha/2, k - 1, n - k))/n.h)
      }
    else
      {
      th.lower <- max(0, ((k - 1)/qchisq(1 - alpha, k - 1))*(ms3 - (ms2*qf(1 - alpha, k - 1, n - k))/n.h))
      th.upper <- Inf
      }
    ci <- c(th.lower, th.upper)
    }

  if (type == 'be')
    {
    if (two.sided == T)
      {
      L <- max(0, (ms3/(ms2*qf((1 - alpha)/2, k - 1,n - k))) - (1/n.min))
      U <- max(0, (ms3/(ms2*qf(alpha/2, k - 1, n - k))) - (1/n.max))

      be.lower <- (n.h*L/(1 + n.h*L))*((k - 1)*ms3)/qchisq((1 - alpha)/2, k - 1)
      be.upper <- (n.h*U/(1 + n.h*U))*((k - 1)*ms3)/qchisq(alpha, k - 1)
      }
    else
      {
      L <- max(0, (ms3/(ms2*qf((1 - alpha), k - 1, n - k))) - (1/n.min))

      be.lower <- (n.h*L/(1 + n.h*L))*((k - 1)*ms3)/qchisq((1 - alpha), k - 1)
      be.upper = Inf
      }
    ci = c(be.lower, be.upper)
    }
  return(ci)
}

coverage.th <- numeric(18)
coverage.be <- numeric(18)
length.th <- numeric(18)
length.be <- numeric(18)
setting = 1

for (k in c(10, 20))
  for (sigt2 in c(0.1, 1, 2))
    for (b in 1:3)
    {
      for (i in 1:1000)
      {
        print(paste(setting, k, sigt2, b, i))
        dataset <- mdata3(n_i = NULL, k = k, b = b, sig2 = 0.5, sigt2 = sigt2, mu = 10)
        ci.th <- approxCI(dataset, type = 'th', two.sided = F)
        ci.be <- approxCI(dataset, type = 'be', two.sided = F)
        if (ci.th[1] < sigt2)
        {
          coverage.th[setting] <- coverage.th[setting] + 1
          length.th[setting] <- length.th[setting] + ci.th[1]
        }
        if (ci.be[1] < sigt2)
        {
          coverage.be[setting] <- coverage.be[setting] + 1
          length.be[setting] <- length.be[setting] + ci.be[1]
        }
      }
      setting <- setting + 1
    }
coverage.be
coverage.th
length.be
length.th

out <- data.frame(coverage.th = coverage.th, length.th = length.th/coverage.th, coverage.be = coverage.be, length.be = length.be/coverage.be)
write.csv(out,'thbe.csv')
