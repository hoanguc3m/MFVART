library(fredr)
setwd("/home/hoanguc3m/Dropbox/WP16/Code/MFVART/vignette/")
Monthly <- read.csv("FredMD2024-10.csv")
Quarterly <- read.csv("FeadQD2024-10.csv")
# Remove Transform:,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,2,2,2,5,5,2,2,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,1,2,1,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,2,6,6,5,6,6,7,6,6,6,2,5,2,5,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,2,6,6,6,1

lag <- 1
t_max <- nrow(Monthly)
Time <- as.Date(Monthly$sasdate, format = "%m/%d/%Y")[(lag+1):t_max]

Var1 <- diff(log(Monthly$INDPRO), lag = lag) * 100
Var2 <- diff(log(Monthly$CPIAUCSL), lag = lag) * 100
#plot(Var2, type ="l")
Var3 <- (Monthly$UNRATE)[(lag+1):t_max]
plot(Var3, type ="l")

Var4 <- log(Monthly$VIXCLSx)[(lag+1):t_max]
plot(Var3, Var4)
Monthly <- data.frame(Time = Time, 
                      IP = Var1, 
                      Inflation = Var2, 
                      UNRATE = Var3,
                      VIX = Var4)

GDP <- Quarterly$GDPC1
GDP_growth <- 100*diff(log(GDP),lag = 1) # Using 1 lag instead of 4 lags.
Time_q <- as.Date(Quarterly$sasdate, format = "%m/%d/%Y")[2:nrow(Quarterly)]

Quarterly <- data.frame(Time = Time_q, GDP_growth = GDP_growth)
Data <- merge(x = Monthly, y = Quarterly, by = "Time", all.x = T)

y <- as.matrix(Data[44:nrow(Data),2:6])
y0 <- as.matrix(Data[43,2:6])
y0[,5] <- Data[41,6]
p = 1
K = 5
t_max <- nrow(y)
vars::VARselect( na.omit(y))

# load("/home/hoanguc3m/MEGA/WP16/RData/03SMN.RData")
# head(y)
library(MFVART)

prior <- get_prior(y, p = p, dist="Gaussian", SV = F, aggregation = "triangular", idq = c(5))
inits <- get_init(prior)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain1 <- BMFVAR.novol(y, K = K, p = p, dist = "Gaussian", y0 = y0, prior = prior, inits = inits)
B <- matrix(apply(Chain1$mcmc$param, MARGIN = 2, FUN = median)[1:30], nrow = 5)
A = diag(1, K, K)
A[upper.tri(A, diag=FALSE)] <- apply(Chain1$mcmc$param, MARGIN = 2, FUN = median)[31:40]
A <- t(A)
A_inv <- solve(A)
A_inv %*% t(A_inv)

prior <- get_prior(y, p = p, dist="Student", SV = F, aggregation = "triangular", idq = c(5))
inits <- get_init(prior)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain2 <- BMFVAR.novol(y, K = K, p = p, dist = "Student", y0 = y0, prior = prior, inits = inits)

plot(apply(Chain1$mcmc$param, MARGIN = 2, FUN = median)[1:45],
     apply(Chain2$mcmc$param, MARGIN = 2, FUN = median)[1:45])
abline(a = 0, b = 1)
apply(Chain2$mcmc$param, MARGIN = 2, FUN = median)[46]

prior <- get_prior(y, p = p, dist="OT", SV = F, aggregation = "triangular", idq = c(5))
inits <- get_init(prior)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain3 <- BMFVAR.novol(y, K = K, p = p, dist = "OT", y0 = y0, prior = prior, inits = inits)
plot(apply(Chain3$mcmc$param, MARGIN = 2, FUN = median)[1:45],
     apply(Chain2$mcmc$param, MARGIN = 2, FUN = median)[1:45])
abline(a = 0, b = 1)
apply(Chain3$mcmc$param, MARGIN = 2, FUN = median)[46:50]

############################################
prior <- get_prior(y, p = p, dist="Gaussian", SV = T, aggregation = "triangular", idq = c(5))
inits <- get_init(prior)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain5 <- BMFVAR.SV(y, K = K, p = p, dist = "Gaussian", y0 = y0, prior = prior, inits = inits)

prior <- get_prior(y, p = p, dist="Student", SV = T, aggregation = "triangular", idq = c(5))
inits <- get_init(prior)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain6 <- BMFVAR.SV(y, K = K, p = p, dist = "Student", y0 = y0, prior = prior, inits = inits)
median(Chain6$mcmc$param[,41])

apply(Chain6$mcmc$param, MARGIN = 2, FUN = median)[1:46]
plot(apply(Chain5$mcmc$param, MARGIN = 2, FUN = median)[1:45],
     apply(Chain6$mcmc$param, MARGIN = 2, FUN = median)[c(1:40,42:46)])
abline(a = 0, b = 1)

prior <- get_prior(y, p = p, dist="OT", SV = T, aggregation = "triangular", idq = c(5))
inits <- get_init(prior)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain7 <- BMFVAR.SV(y, K = K, p = p, dist = "OT", y0 = y0, prior = prior, inits = inits)
apply(Chain7$mcmc$param, MARGIN = 2, FUN = median)[41:45]
# nu1       nu2       nu3       nu4       nu5 
# 14.367142 19.200400  8.738173  5.907335 19.278837 

B <- matrix(apply(Chain7$mcmc$param, MARGIN = 2, FUN = median)[1:30], nrow = 5)
A = diag(1, K, K)
A[upper.tri(A, diag=FALSE)] <- apply(Chain7$mcmc$param, MARGIN = 2, FUN = median)[31:40]
A <- t(A)
A_inv <- solve(A)
A_inv %*% t(A_inv)

plot(apply(Chain6$mcmc$param, MARGIN = 2, FUN = median)[c(1:40,42:46)],
     apply(Chain7$mcmc$param, MARGIN = 2, FUN = median)[c(1:40,46:50)])
abline(a = 0, b = 1)

# save.image("/home/hoanguc3m/MEGA/WP16/RData/03KMN.RData")

prior <- get_prior(y, p = p, dist="Gaussian", SV = T, 
                   aggregation = "triangular", idq = c(5), r = 2)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain8 <- BMFVAR.FSV(y, K = K, p = p, dist = "Gaussian", y0 = y0, prior = prior, inits = inits)

plot(apply(get_post(Chain5$mcmc$param, element = "B"), 2, mean),
     apply(get_post(Chain8$mcmc$param, element = "B"), 2, mean))
abline(a = 0, b = 1)
plot(apply((Chain5$mcmc$y_miss), 2, mean),
     apply(Chain8$mcmc$y_miss, 2, mean))

prior <- get_prior(y, p = p, dist="Student", SV = T, 
                   aggregation = "triangular", idq = c(5), r = 2)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain9 <- BMFVAR.FSV(y, K = K, p = p, dist = "Student", y0 = y0, prior = prior, inits = inits)
plot(apply(get_post(Chain6$mcmc$param, element = "B"), 2, mean),
     apply(get_post(Chain9$mcmc$param, element = "B"), 2, mean))
abline(a = 0, b = 1)
plot(apply((Chain6$mcmc$y_miss), 2, mean),
     apply(Chain9$mcmc$y_miss, 2, mean))

prior <- get_prior(y, p = p, dist="OT", SV = T, 
                   aggregation = "triangular", idq = c(5), r = 2)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain10 <- BMFVAR.FSV(y, K = K, p = p, dist = "OT", y0 = y0, prior = prior, inits = inits)

plot(apply(get_post(Chain7$mcmc$param, element = "B"), 2, mean),
     apply(get_post(Chain10$mcmc$param, element = "B"), 2, mean))
abline(a = 0, b = 1)
plot(apply((Chain7$mcmc$y_miss), 2, mean),
     apply(Chain10$mcmc$y_miss, 2, mean))

save.image("/home/hoanguc3m/MEGA/WP16/RData/03KMN.RData")
