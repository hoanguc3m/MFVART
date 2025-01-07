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

Var1 <- diff(log(Monthly$INDPRO), lag = lag) * 100
Var2 <- diff(log(Monthly$CPIAUCSL), lag = lag) * 100
#plot(Var2, type ="l")
Var3 <- (Monthly$UNRATE)[(lag+1):t_max]
plot(Var3, type ="l")

Var4 <- log(Monthly$VIXCLSx)[(lag+1):t_max]
Monthly <- data.frame(Time = Time, 
                      IP = Var1, 
                      Inflation = Var2, 
                      UNRATE = Var3,
                      VIX = Var4)

# GDP <- Quarterly$GDPC1
# GDP_growth <- 100*diff(log(GDP),lag = 1) # Using 1 lag instead of 4 lags.
# Time_q <- as.Date(Quarterly$sasdate, format = "%m/%d/%Y")[2:nrow(Quarterly)]
# 
# Quarterly <- data.frame(Time = Time_q, GDP_growth = GDP_growth)
# Data <- merge(x = Monthly, y = Quarterly, by = "Time", all.x = T)
Data <- Monthly
y <- as.matrix(Data[44:nrow(Data),2:5])
y0 <- as.matrix(Data[42:43,2:5])
p = 2
K = 4
t_max <- nrow(y)
vars::VARselect( na.omit(y))

library(MFVART)
prior <- fatBVARS::get_prior(y, p = p, dist="Gaussian", SV = F)
prior$r <- 2; prior$Vlhyper <- 1; prior$mulhyper <- 0; 
inits <- fatBVARS::get_init(prior)
inits <- fatBVARS::get_init(prior, samples = 10000, burnin = 5000, thin = 1)
Chain1 <- JCVARF.novol(y, K = K, p = p, dist = "Gaussian", y0 = y0, prior = prior, inits = inits)
Chain2 <- fatBVARS::BVAR.novol(y, K = K, p = p, dist = "Gaussian", y0 = y0, prior = prior, inits = inits)
C1_Bmat <- get_post(Chain1$mcmc, element = "B")
C2_Bmat <- get_post(Chain2$mcmc, element = "B")
plot(colMeans(C1_Bmat), colMeans(C2_Bmat))
abline(a=0,b = 1)
C1_Amat <- get_post(Chain1$mcmc, element = "L")
C2_Amat <- get_post(Chain2$mcmc, element = "a")
colMeans(C1_Amat)
colMeans(C2_Amat)
C1_Smat <- get_post(Chain1$mcmc, element = "sigma")
C2_Smat <- get_post(Chain2$mcmc, element = "sigma")
colMeans(C1_Smat)
colMeans(C2_Smat)


prior <- fatBVARS::get_prior(y, p = p, dist="Student", SV = F)
prior$r <- 2; prior$Vlhyper <- 1; prior$mulhyper <- 0; 
inits <- fatBVARS::get_init(prior)
inits <- fatBVARS::get_init(prior, samples = 10000, burnin = 5000, thin = 1)
Chain3 <- JCVARF.novol(y, K = K, p = p, dist = "Student", y0 = y0, prior = prior, inits = inits)
Chain4 <- fatBVARS::BVAR.novol(y, K = K, p = p, dist = "Student", y0 = y0, prior = prior, inits = inits)
C1_Bmat <- get_post(Chain3$mcmc, element = "B")
C2_Bmat <- get_post(Chain4$mcmc, element = "B")
plot(colMeans(C1_Bmat), colMeans(C2_Bmat))
abline(a=0,b = 1)
cbind(colMeans(C1_Bmat), colMeans(C2_Bmat), abs(colMeans(C1_Bmat)- colMeans(C2_Bmat)))
C1_Amat <- get_post(Chain3$mcmc, element = "L")
C2_Amat <- get_post(Chain4$mcmc, element = "a")
colMeans(C1_Amat)
colMeans(C2_Amat)
mean(get_post(Chain3$mcmc, element = "nu"))
mean(get_post(Chain4$mcmc, element = "nu"))
C1_Smat <- get_post(Chain3$mcmc, element = "sigma")
C2_Smat <- get_post(Chain4$mcmc, element = "sigma")
colMeans(C1_Smat)
colMeans(C2_Smat)

prior <- fatBVARS::get_prior(y, p = p, dist="OT", SV = F)
prior$r <- 2; prior$Vlhyper <- 1; prior$mulhyper <- 0; 
inits <- fatBVARS::get_init(prior)
inits <- fatBVARS::get_init(prior, samples = 10000, burnin = 5000, thin = 1)
Chain5 <- JCVARF.novol(y, K = K, p = p, dist = "OT", y0 = y0, prior = prior, inits = inits)
Chain6 <- fatBVARS::BVAR.novol(y, K = K, p = p, dist = "OT", y0 = y0, prior = prior, inits = inits)
C1_Bmat <- get_post(Chain5$mcmc, element = "B")
C2_Bmat <- get_post(Chain6$mcmc, element = "B")
plot(colMeans(C1_Bmat), colMeans(C2_Bmat))
abline(a=0,b = 1)
C1_Amat <- get_post(Chain5$mcmc, element = "L")
C2_Amat <- get_post(Chain6$mcmc, element = "a")
colMeans(C1_Amat)
colMeans(C2_Amat)
colMeans(get_post(Chain5$mcmc, element = "nu"))
colMeans(get_post(Chain6$mcmc, element = "nu"))



library(MFVART)
prior <- fatBVARS::get_prior(y, p = p, dist="Gaussian", SV = T)
prior$r <- 2; prior$Vlhyper <- 1; prior$mulhyper <- 0; 
inits <- fatBVARS::get_init(prior)
inits <- fatBVARS::get_init(prior, samples = 10000, burnin = 5000, thin = 1)
Chain1 <- JCVARF.SV(y, K = K, p = p, dist = "Gaussian", y0 = y0, prior = prior, inits = inits)
Chain2 <- fatBVARS::BVAR.SV(y, K = K, p = p, dist = "Gaussian", y0 = y0, prior = prior, inits = inits)
C1_Bmat <- get_post(Chain1$mcmc, element = "B")
C2_Bmat <- get_post(Chain2$mcmc, element = "B")
plot(colMeans(C1_Bmat), colMeans(C2_Bmat))
abline(a=0,b = 1)
C1_Amat <- get_post(Chain1$mcmc, element = "L")
C2_Amat <- get_post(Chain2$mcmc, element = "a")
colMeans(C1_Amat)
colMeans(C2_Amat)
C1_Smat <- get_post(Chain1$mcmc, element = "sigma")
C2_Smat <- get_post(Chain2$mcmc, element = "sigma")
colMeans(C1_Smat)
colMeans(C2_Smat)


prior <- fatBVARS::get_prior(y, p = p, dist="Student", SV = T)
prior$r <- 2; prior$Vlhyper <- 1; prior$mulhyper <- 0; 
inits <- fatBVARS::get_init(prior)
inits <- fatBVARS::get_init(prior, samples = 10000, burnin = 5000, thin = 1)
Chain3 <- JCVARF.SV(y, K = K, p = p, dist = "Student", y0 = y0, prior = prior, inits = inits)
Chain4 <- fatBVARS::BVAR.SV(y, K = K, p = p, dist = "Student", y0 = y0, prior = prior, inits = inits)
C1_Bmat <- get_post(Chain3$mcmc, element = "B")
C2_Bmat <- get_post(Chain4$mcmc, element = "B")
plot(colMeans(C1_Bmat), colMeans(C2_Bmat))
abline(a=0,b = 1)
C1_Amat <- get_post(Chain3$mcmc, element = "L")
C2_Amat <- get_post(Chain4$mcmc, element = "a")
colMeans(C1_Amat)
colMeans(C2_Amat)
C1_Smat <- get_post(Chain3$mcmc, element = "sigma")
C2_Smat <- get_post(Chain4$mcmc, element = "sigma")
colMeans(C1_Smat)
colMeans(C2_Smat)
mean(get_post(Chain3$mcmc, element = "nu"))
mean(get_post(Chain4$mcmc, element = "nu"))

prior <- fatBVARS::get_prior(y, p = p, dist="OT", SV = T)
prior$r <- 2; prior$Vlhyper <- 1; prior$mulhyper <- 0; 
inits <- fatBVARS::get_init(prior)
inits <- fatBVARS::get_init(prior, samples = 10000, burnin = 5000, thin = 1)
Chain5 <- JCVARF.SV(y, K = K, p = p, dist = "OT", y0 = y0, prior = prior, inits = inits)
Chain6 <- fatBVARS::BVAR.SV(y, K = K, p = p, dist = "OT", y0 = y0, prior = prior, inits = inits)
C1_Bmat <- get_post(Chain5$mcmc, element = "B")
C2_Bmat <- get_post(Chain6$mcmc, element = "B")
plot(colMeans(C1_Bmat), colMeans(C2_Bmat))
abline(a=0,b = 1)
C1_Amat <- get_post(Chain5$mcmc, element = "L")
C2_Amat <- get_post(Chain6$mcmc, element = "a")
colMeans(C1_Amat)
colMeans(C2_Amat)
C1_Smat <- get_post(Chain5$mcmc, element = "sigma")
C2_Smat <- get_post(Chain6$mcmc, element = "sigma")
colMeans(C1_Smat)
colMeans(C2_Smat)
colMeans(get_post(Chain5$mcmc, element = "nu"))
colMeans(get_post(Chain6$mcmc, element = "nu"))
