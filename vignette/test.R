library(MFVART)
datagen <- sim.MFVAR.novol(dist="Gaussian", idq = c(4), p = 2, K = 4, t_max = 500)
# datagen <- sim.MFVAR.SV(dist="OT", idq = c(5))
# 
# datagen <- sim.MFVAR.SV(dist="Student", idq = c(5))

head(datagen$y)
head(datagen$y0)
tail(datagen$y)
y = datagen$y
y0 = datagen$y0
p = datagen$p
K = datagen$K

prior <- get_prior(y, p = p, dist="Gaussian", SV = F)
inits <- get_init(prior)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain1 <- BMFVAR.novol(y, K = K, p = p, dist = "Gaussian", y0 = y0, prior = prior, inits = inits)
hist(Chain1$mcmc$y_miss[,1])
datagen$y_true[is.na(datagen$y[,5]),5]

plot( datagen$y_true[is.na(datagen$y[,5]),5], apply(Chain1$mcmc$y_miss, MARGIN = 2, FUN = mean))
cor( datagen$y_true[is.na(datagen$y[,5]),5], apply(Chain1$mcmc$y_miss, MARGIN = 2, FUN = mean))
# cor( datagen$y_true[is.na(datagen$y[,5]),5], apply(My_miss, MARGIN = 1, FUN = mean))
# plot( datagen$y_true[is.na(datagen$y[,5]),5], apply(My_miss, MARGIN = 1, FUN = mean))
abline(a = 0, b = 1)

inits$A0 <-  datagen$A0
inits$B0 <-  datagen$B0
inits$sigma <- as.numeric(exp(0.5*datagen$h))
b_sample <- vec(inits$B0);
a_sample <- t(inits$A0)[upper.tri(inits$A0)]
sigma <- as.numeric(exp(0.5*datagen$h))
