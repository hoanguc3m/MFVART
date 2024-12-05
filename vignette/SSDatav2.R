library(fredr)
setwd("/home/hoanguc3m/Dropbox/WP16/Code/MFVART/Data/")
Monthly <- read.csv("FredMD2024-10.csv")
Quarterly <- read.csv("FredQD2024-10.csv")
# Remove Transform:,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,2,2,2,5,5,2,2,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,1,2,1,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,2,6,6,5,6,6,7,6,6,6,2,5,2,5,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,2,6,6,6,1

lag <- 12
t_max <- nrow(Monthly)
Time <- as.Date(Monthly$sasdate, format = "%m/%d/%Y")[(lag+1):t_max]

UNR <- Monthly$UNRATE[(lag+1):t_max]
HRS <- Monthly$AWHMAN[(lag+1):t_max]
CPI <- diff(log(Monthly$CPIAUCSL), lag = lag) * 100
IP <- diff(log(Monthly$INDPRO), lag = lag) * 100
PCE <- diff(log(Monthly$PCEPI), lag = lag) * 100
FF <- Monthly$FEDFUNDS[(lag+1):t_max]
TB <- Monthly$GS10[(lag+1):t_max]
SP500 <- tail(diff(log(Monthly$S.P.500), lag = 1) * 100, t_max-lag)

Monthly <- data.frame(Time = Time, 
                      UNR = UNR, HRS = HRS, CPI = CPI,
                      IP = IP, PCE = PCE, FF = FF, 
                      TB = TB, SP500 = SP500)

GDP <- Quarterly$GDPC1
GDP_growth <- 100*diff(log(GDP),lag = 1) # Using 1 lag instead of 4 lags.
Time_q <- as.Date(Quarterly$sasdate, format = "%m/%d/%Y")[2:nrow(Quarterly)]

FPIx <- Quarterly$FPIx
INVFIX <- 100*diff(log(FPIx),lag = 1) # Using 1 lag instead of 4 lags.

GCEC <- Quarterly$GCEC1
GOV <- 100*diff(log(GCEC),lag = 1) # Using 1 lag instead of 4 lags.

Quarterly <- data.frame(Time = Time_q, GDP = GDP_growth, 
                        INVFIX = INVFIX, GOV = GOV)
Data <- merge(x = Monthly, y = Quarterly, by = "Time", all.x = T)


y <- as.matrix(Data[9:nrow(Data),2:12])
y0 <- as.matrix(Data[7:8,2:12])
y0[1,c(9:11)] <- as.numeric(Data[3,10:12])
y0[2,c(9:11)] <- as.numeric(Data[6,10:12])
p = 2
K = 11
t_max <- nrow(y)
vars::VARselect( na.omit(y))

# load("/home/hoanguc3m/MEGA/WP16/RData/03SSv2.RData")
library(MFVART)

prior <- get_prior(y, p = p, dist="Gaussian", SV = F, aggregation = "triangular", idq = c(9,10,11))
inits <- get_init(prior)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain1 <- BMFVAR.novol(y, K = K, p = p, dist = "Gaussian", y0 = y0, prior = prior, inits = inits)
2*11*11+11
B <- matrix(apply(Chain1$mcmc$param, MARGIN = 2, FUN = median)[1:253], nrow = 11)
A = diag(1, K, K)
2*11*11+11+10*11/2
A[upper.tri(A, diag=FALSE)] <- apply(Chain1$mcmc$param, MARGIN = 2, FUN = median)[254:308]
A <- t(A)

prior <- get_prior(y, p = p, dist="Student", SV = F, aggregation = "triangular", idq = c(9,10,11))
inits <- get_init(prior)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain2 <- BMFVAR.novol(y, K = K, p = p, dist = "Student", y0 = y0, prior = prior, inits = inits)
2*11*11+11+10*11/2+11
apply(Chain2$mcmc$param, MARGIN = 2, FUN = median)[320]
plot(apply(Chain1$mcmc$param, MARGIN = 2, FUN = median)[1:319],
     apply(Chain2$mcmc$param, MARGIN = 2, FUN = median)[1:319])
abline(a = 0, b = 1)

prior <- get_prior(y, p = p, dist="OT", SV = F, aggregation = "triangular", idq = c(9,10,11))
inits <- get_init(prior)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain3 <- BMFVAR.novol(y, K = K, p = p, dist = "OT", y0 = y0, prior = prior, inits = inits)
apply(Chain3$mcmc$param, MARGIN = 2, FUN = median)[320:330]
plot(apply(Chain3$mcmc$param, MARGIN = 2, FUN = median)[1:319],
     apply(Chain2$mcmc$param, MARGIN = 2, FUN = median)[1:319])
abline(a = 0, b = 1)
# cbind(apply(Chain3$mcmc$param, MARGIN = 2, FUN = median)[1:319],
#       apply(Chain2$mcmc$param, MARGIN = 2, FUN = median)[1:319])

############################################
prior <- get_prior(y, p = p, dist="Gaussian", SV = T, aggregation = "triangular", idq = c(9,10,11))
inits <- get_init(prior)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain5 <- BMFVAR.SV(y, K = K, p = p, dist = "Gaussian", y0 = y0, prior = prior, inits = inits)

prior <- get_prior(y, p = p, dist="Student", SV = T, aggregation = "triangular", idq = c(9,10,11))
inits <- get_init(prior)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain6 <- BMFVAR.SV(y, K = K, p = p, dist = "Student", y0 = y0, prior = prior, inits = inits)
median(Chain6$mcmc$param[,320])
plot(apply(Chain5$mcmc$param, MARGIN = 2, FUN = median)[1:319],
     apply(Chain6$mcmc$param, MARGIN = 2, FUN = median)[c(1:26,28:31)])
abline(a = 0, b = 1)

apply(Chain6$mcmc$param, MARGIN = 2, FUN = median)[1:31]
plot(apply(Chain5$mcmc$param, MARGIN = 2, FUN = median)[1:30],
     apply(Chain6$mcmc$param, MARGIN = 2, FUN = median)[c(1:26,28:31)])
abline(a = 0, b = 1)

prior <- get_prior(y, p = p, dist="OT", SV = T, aggregation = "triangular", idq = c(9,10,11))
inits <- get_init(prior)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain7 <- BMFVAR.SV(y, K = K, p = p, dist = "OT", y0 = y0, prior = prior, inits = inits)
apply(Chain7$mcmc$param, MARGIN = 2, FUN = median)[309:319]
# nu1       nu2       nu3       nu4       nu5       nu6       nu7       nu8       nu9 
# 8.508969  4.346246 18.090924 20.060630 27.558793 17.616148 28.698023  5.691893 17.160504 
# nu10      nu11 
# 22.793705 21.217010 

B <- matrix(apply(Chain7$mcmc$param, MARGIN = 2, FUN = median)[1:253], nrow = 11)
A = diag(1, K, K)
A[upper.tri(A, diag=FALSE)] <- apply(Chain7$mcmc$param, MARGIN = 2, FUN = median)[21:26]
A <- t(A)

plot(apply(Chain6$mcmc$param, MARGIN = 2, FUN = median)[c(1:308,310:320)],
     apply(Chain7$mcmc$param, MARGIN = 2, FUN = median)[c(1:308,320:330)])
abline(a = 0, b = 1)

####################################



prior <- get_prior(y, p = p, dist="Gaussian", SV = T, 
                   aggregation = "triangular", idq = c(9,10,11), r = 4)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain8 <- BMFVAR.FSV(y, K = K, p = p, dist = "Gaussian", y0 = y0, prior = prior, inits = inits)

plot(apply(get_post(Chain5$mcmc$param, element = "B"), 2, mean),
     apply(get_post(Chain8$mcmc$param, element = "B"), 2, mean))
abline(a = 0, b = 1)
plot(apply((Chain5$mcmc$y_miss), 2, mean),
     apply(Chain8$mcmc$y_miss, 2, mean))

prior <- get_prior(y, p = p, dist="Student", SV = T, 
                   aggregation = "triangular", idq = c(9,10,11), r = 4)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain9 <- BMFVAR.FSV(y, K = K, p = p, dist = "Student", y0 = y0, prior = prior, inits = inits)
plot(apply(get_post(Chain6$mcmc$param, element = "B"), 2, mean),
     apply(get_post(Chain9$mcmc$param, element = "B"), 2, mean))
abline(a = 0, b = 1)
plot(apply((Chain6$mcmc$y_miss), 2, mean),
     apply(Chain9$mcmc$y_miss, 2, mean))

prior <- get_prior(y, p = p, dist="OT", SV = T, 
                   aggregation = "triangular", idq = c(9,10,11), r = 4)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain10 <- BMFVAR.FSV(y, K = K, p = p, dist = "OT", y0 = y0, prior = prior, inits = inits)
plot(apply(get_post(Chain7$mcmc$param, element = "B"), 2, mean),
     apply(get_post(Chain10$mcmc$param, element = "B"), 2, mean))
abline(a = 0, b = 1)
plot(apply((Chain7$mcmc$y_miss), 2, mean),
     apply(Chain10$mcmc$y_miss, 2, mean))

plot(Chain1, element = "y_9")
plot(Chain5, element = "y_9")
plot(Chain8, element = "y_9")

plot(Chain2, element = "y_9")
plot(Chain6, element = "y_9")
plot(Chain9, element = "y_9")

plot(Chain3, element = "y_9")
plot(Chain7, element = "y_9")
plot(Chain10, element = "y_9")

c(Chain5$esttime, Chain6$esttime, Chain7$esttime)
c(Chain8$esttime, Chain9$esttime, Chain10$esttime)

save.image("/home/hoanguc3m/MEGA/WP16/RData/03SSv2.RData")

plot(Chain1, element = "y_10")
plot(Chain5, element = "y_10")
plot(Chain8, element = "y_10")

plot(Chain2, element = "y_10")
plot(Chain6, element = "y_10")
plot(Chain9, element = "y_10")

plot(Chain3, element = "y_10")
plot(Chain7, element = "y_10")
plot(Chain10, element = "y_10")

plot(Chain1, element = "y_11")
plot(Chain5, element = "y_11")
plot(Chain8, element = "y_11")

plot(Chain2, element = "y_11")
plot(Chain6, element = "y_11")
plot(Chain9, element = "y_11")

plot(Chain3, element = "y_11")
plot(Chain7, element = "y_11")
plot(Chain10, element = "y_11")
