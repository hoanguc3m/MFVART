library(fredr)
setwd("/home/hoanguc3m/Dropbox/WP16/Code/MFVART/Data/")
DatM_name <- "FredQD2024-10.csv"
DatQ_name <- "FredMD2024-10.csv"
headers <- read.csv(DatM_name, header = F, nrows = 1, as.is = T)
Monthly <- read.csv(DatM_name, skip = 2, header = F)
colnames(Monthly) <- headers

headers <- read.csv(DatQ_name, header = F, nrows = 1, as.is = T)
Quarterly <- read.csv(DatQ_name, skip = 3, header = F)
colnames(Quarterly) <- headers
# Remove Transform:,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,2,2,2,5,5,2,2,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,1,2,1,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,2,6,6,5,6,6,7,6,6,6,2,5,2,5,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,2,6,6,6,1

lag <- 1
t_max <- nrow(Monthly)
Time <- as.Date(Monthly$sasdate, format = "%m/%d/%Y")[(lag+1):t_max]

UNR <- diff(Monthly$UNRATE, lag = lag)
HRS <- Monthly$AWHMAN[(lag+1):t_max]
CPI <- diff(log(Monthly$CPIAUCSL), lag = lag) * 100
IP <- diff(log(Monthly$INDPRO), lag = lag) * 100
PCE <- diff(log(Monthly$PCEPI), lag = lag) * 100
FF <- diff(Monthly$FEDFUNDS, lag = lag)
TB <- diff(Monthly$GS10, lag = lag)
SP500 <- diff(log(Monthly$S.P.500), lag = lag) * 100

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


y <- as.matrix(Data[11:nrow(Data),2:12])
y0 <- as.matrix(Data[9:10,2:12])
y0[1,c(9:11)] <- as.numeric(Data[5,10:12])
y0[2,c(9:11)] <- as.numeric(Data[8,10:12])
p = 2
K = 11
t_max <- nrow(y)
vars::VARselect( na.omit(y))

# load("/home/hoanguc3m/MEGA/WP16/RData/03SS.RData")
head(y)
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
cbind(apply(Chain3$mcmc$param, MARGIN = 2, FUN = median)[1:319],
     apply(Chain2$mcmc$param, MARGIN = 2, FUN = median)[1:319])

############################################
prior <- get_prior(y, p = p, dist="Gaussian", SV = T, aggregation = "triangular", idq = c(9,10,11))
inits <- get_init(prior)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain5 <- BMFVAR.SV(y, K = K, p = p, dist = "Gaussian", y0 = y0, prior = prior, inits = inits)
# Chain5 can not be estimated due to chol matrix error

prior <- get_prior(y, p = p, dist="Student", SV = T, aggregation = "triangular", idq = c(9,10,11))
inits <- get_init(prior)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain6 <- BMFVAR.SV(y, K = K, p = p, dist = "Student", y0 = y0, prior = prior, inits = inits)
2*11*11+11+10*11/2+11

median(Chain6$mcmc$param[,309])

# plot(apply(Chain5$mcmc$param, MARGIN = 2, FUN = median)[1:319],
#      apply(Chain6$mcmc$param, MARGIN = 2, FUN = median)[c(1:26,28:31)])
# abline(a = 0, b = 1)

prior <- get_prior(y, p = p, dist="OT", SV = T, aggregation = "triangular", idq = c(9,10,11))
inits <- get_init(prior)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain7 <- BMFVAR.SV(y, K = K, p = p, dist = "OT", y0 = y0, prior = prior, inits = inits)
apply(Chain7$mcmc$param, MARGIN = 2, FUN = median)[309:319]
# nu1       nu2       nu3       nu4       nu5       nu6       nu7       nu8       nu9 
# 8.358775  4.386976 24.117296  5.437233 20.645196 14.272021 28.764083  6.105962 24.065667 
# nu10      nu11 
# 5.515587 26.645618 

B <- matrix(apply(Chain7$mcmc$param, MARGIN = 2, FUN = median)[1:253], nrow = 11)
A = diag(1, K, K)
A[upper.tri(A, diag=FALSE)] <- apply(Chain7$mcmc$param, MARGIN = 2, FUN = median)[21:26]
A <- t(A)


plot(apply(Chain6$mcmc$param, MARGIN = 2, FUN = median)[c(1:308,310:320)],
     apply(Chain7$mcmc$param, MARGIN = 2, FUN = median)[c(1:308,320:330)])
abline(a = 0, b = 1)

save.image("/home/hoanguc3m/MEGA/WP16/RData/03SS.RData")