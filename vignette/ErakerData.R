library(fredr)
fredr_set_key("49698c425172d43d2058f3baa69f8fe0")

# 3-Month Treasury Bill Secondary Market Rate, Discount Basis
# TB3MS
# 
# Consumer Price Index for All Urban Consumers: All Items in U.S. City Average
# CPIAUCSL
# 
# Market Yield on U.S. Treasury Securities at 20-Year Constant Maturity, Quoted on an Investment Basis
# DGS20
# 
# Moody's Seasoned Aaa Corporate Bond Yield
# AAA
# 
# Moody's Seasoned Baa Corporate Bond Yield
# DBAA
# 
# Gross Domestic Product (GDP)
# GDP


TB3MS <- fredr(series_id = "TB3MS", # Industrial Production: Total Index
                observation_start = as.Date("1920-01-01"),
                observation_end = as.Date("2024-10-01"),
                frequency = "m",
                aggregation_method = "avg")

CPIAUCSL <- fredr(series_id = "CPIAUCSL", # Consumer Price Index for All Urban Consumers: All Items in U.S. City Average
                  observation_start = as.Date("1920-01-01"),
                  observation_end = as.Date("2024-10-01"),
                  frequency = "m",
                  aggregation_method = "avg")
logInfIndex <- log(CPIAUCSL$value)
Inflation <- 100*(logInfIndex[13:length(CPIAUCSL$value)] - 
               logInfIndex[1:(length(CPIAUCSL$value)-12)])


GS10 <- fredr(series_id = "GS10", # Treasury Securities at 10-Year Constant Maturity
                observation_start = as.Date("1920-01-01"),
                observation_end = as.Date("2024-10-01"),
                frequency = "m",
                aggregation_method = "avg")

AAA <- fredr(series_id = "AAA", # Moody's Seasoned AAA Corporate Bond Yield
               observation_start = as.Date("1920-01-01"),
               observation_end = as.Date("2024-10-01"),
               frequency = "m",
               aggregation_method = "avg")

BAA <- fredr(series_id = "BAA", # Moody's Seasoned BAA Corporate Bond Yield
             observation_start = as.Date("1920-01-01"),
             observation_end = as.Date("2024-10-01"),
             frequency = "m",
             aggregation_method = "avg")

GDP <- fredr(series_id = "GDPC1", # Gross Domestic Product (GDP)
             observation_start = as.Date("1920-01-01"),
             observation_end = as.Date("2024-10-01"),
             frequency = "q",
             aggregation_method = "avg")

t_max <- 859
Time <- tail(TB3MS$date, t_max)
RealRate <- tail(TB3MS$value, t_max) - tail(Inflation, t_max)
Slope <- tail(TB3MS$value, t_max) - tail(GS10$value, t_max)
Default_spread <- tail(BAA$value - AAA$value, t_max)
GDP_growth <- 100*diff(log(GDP$value),lag = 1) # Using 1 lag instead of 4 lags.
# GDP_growth <- GDP_growth[5:length(GDP_growth)]

Time_q <- tail(GDP$date, length(GDP_growth))

plot(tail(CPIAUCSL$date,t_max), RealRate, type = "l")
plot(tail(CPIAUCSL$date,t_max), Slope, type = "l")
plot(tail(CPIAUCSL$date,t_max), Default_spread, type = "l")

# Monthly <- as.data.frame(cbind(Time, RealRate, Slope, Default_spread))
                      
                      
Monthly <- data.frame(Time = Time, 
                         RealRate = RealRate, 
                         Slope = Slope, 
                         Default_spread = Default_spread)
Quarterly <- data.frame(Time = Time_q, GDP_growth = GDP_growth)

Data <- merge(x = Monthly, y = Quarterly, by = "Time", all.x = T)

y <- as.matrix(Data[4:nrow(Data),2:5])
y0 <- as.matrix(Data[3,2:5])
y0[,4] <- Data[1,5]
p = 1
K = 4
t_max <- nrow(y)
vars::VARselect( na.omit(y))

# load("/home/hoanguc3m/MEGA/WP16/RData/03Eraker.RData")
# head(y)

library(MFVART)

prior <- get_prior(y, p = p, dist="Gaussian", SV = F, aggregation = "triangular", idq = c(4))
inits <- get_init(prior)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain1 <- BMFVAR.novol(y, K = K, p = p, dist = "Gaussian", y0 = y0, prior = prior, inits = inits)
B <- matrix(apply(Chain1$mcmc$param, MARGIN = 2, FUN = median)[1:20], nrow = 4)
A = diag(1, K, K)
A[upper.tri(A, diag=FALSE)] <- apply(Chain7$mcmc$param, MARGIN = 2, FUN = median)[21:26]
A <- t(A)
A_inv <- solve(A)
A_inv %*% t(A_inv)

prior <- get_prior(y, p = p, dist="Student", SV = F, aggregation = "triangular", idq = c(4))
inits <- get_init(prior)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain2 <- BMFVAR.novol(y, K = K, p = p, dist = "Student", y0 = y0, prior = prior, inits = inits)
apply(Chain2$mcmc$param, MARGIN = 2, FUN = median)[1:31]
plot(apply(Chain1$mcmc$param, MARGIN = 2, FUN = median)[1:30],
     apply(Chain2$mcmc$param, MARGIN = 2, FUN = median)[1:30])
abline(a = 0, b = 1)
cbind(apply(Chain1$mcmc$param, MARGIN = 2, FUN = median)[1:30],
     apply(Chain2$mcmc$param, MARGIN = 2, FUN = median)[1:30])


prior <- get_prior(y, p = p, dist="OT", SV = F, aggregation = "triangular", idq = c(4))
inits <- get_init(prior)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain3 <- BMFVAR.novol(y, K = K, p = p, dist = "OT", y0 = y0, prior = prior, inits = inits)
apply(Chain3$mcmc$param, MARGIN = 2, FUN = median)[1:34]
plot(apply(Chain3$mcmc$param, MARGIN = 2, FUN = median)[1:30],
     apply(Chain2$mcmc$param, MARGIN = 2, FUN = median)[1:30])
abline(a = 0, b = 1)

############################################
prior <- get_prior(y, p = p, dist="Gaussian", SV = T, aggregation = "triangular", idq = c(4))
inits <- get_init(prior)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain5 <- BMFVAR.SV(y, K = K, p = p, dist = "Gaussian", y0 = y0, prior = prior, inits = inits)

prior <- get_prior(y, p = p, dist="Student", SV = T, aggregation = "triangular", idq = c(4))
inits <- get_init(prior)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain6 <- BMFVAR.SV(y, K = K, p = p, dist = "Student", y0 = y0, prior = prior, inits = inits)
median(Chain6$mcmc$param[,27])

apply(Chain6$mcmc$param, MARGIN = 2, FUN = median)[1:31]
plot(apply(Chain5$mcmc$param, MARGIN = 2, FUN = median)[1:30],
     apply(Chain6$mcmc$param, MARGIN = 2, FUN = median)[c(1:26,28:31)])
abline(a = 0, b = 1)

prior <- get_prior(y, p = p, dist="OT", SV = T, aggregation = "triangular", idq = c(4))
inits <- get_init(prior)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain7 <- BMFVAR.SV(y, K = K, p = p, dist = "OT", y0 = y0, prior = prior, inits = inits)
apply(Chain7$mcmc$param, MARGIN = 2, FUN = median)[27:30]
# nu1      nu2      nu3      nu4 
# 21.04975 27.83922 36.65009 14.40485 

B <- matrix(apply(Chain7$mcmc$param, MARGIN = 2, FUN = median)[1:20], nrow = 4)
A = diag(1, K, K)
A[upper.tri(A, diag=FALSE)] <- apply(Chain7$mcmc$param, MARGIN = 2, FUN = median)[21:26]
A <- t(A)


plot(apply(Chain6$mcmc$param, MARGIN = 2, FUN = median)[c(1:26,28:31)],
     apply(Chain7$mcmc$param, MARGIN = 2, FUN = median)[c(1:26,31:34)])
abline(a = 0, b = 1)

save.image("/home/hoanguc3m/MEGA/WP16/RData/03Eraker.RData")