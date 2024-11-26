library(MFVART)
#############################################
set.seed(0)
datagen <- sim.MFVAR.novol(dist="Gaussian", idq = c(4), p = 2, K = 4, t_max = 500, aggregation = "triangular")

head(datagen$y)
head(datagen$y0)
tail(datagen$y)
y = datagen$y
y0 = datagen$y0
p = datagen$p
K = datagen$K

prior <- get_prior(y, p = p, dist="Gaussian", SV = F, aggregation = "triangular", idq = c(4))
inits <- get_init(prior)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain1 <- BMFVAR.novol(y, K = K, p = p, dist = "Gaussian", 
                       y0 = y0, prior = prior, inits = inits)

prior <- get_prior(y, p = p, dist="Gaussian", SV = F, 
                   aggregation = "triangular", idq = c(4), r = 2)
inits <- get_init(prior)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain2 <- BMFVAR.FHV(y, K = K, p = p, dist = "Gaussian", 
                     y0 = y0, prior = prior, inits = inits)
plot(apply(get_post(Chain1$mcmc$param, element = "B"), 2, mean),
      apply(get_post(Chain2$mcmc$param, element = "B"), 2, mean))
abline(a = 0, b = 1)
plot(vec(datagen$B0),
     apply(get_post(Chain1$mcmc$param, element = "B"), 2, mean))
plot(vec(datagen$B0),
     apply(get_post(Chain2$mcmc$param, element = "B"), 2, mean))

plot(apply((Chain1$mcmc$y_miss), 2, mean),
     apply(Chain2$mcmc$y_miss, 2, mean))
abline(a = 0, b = 1)

apply(get_post(Chain1$mcmc$param, element = "a"), 2, mean)
apply(get_post(Chain2$mcmc$param, element = "L"), 2, mean)

###############

library(MFVART)

set.seed(0)

datagen <- sim.MFVAR.novol(dist="Student", idq = c(4), p = 2, K = 4, t_max = 500, aggregation = "triangular")
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

prior <- get_prior(y, p = p, dist="Student", SV = F, aggregation = "triangular", idq = c(4))
inits <- get_init(prior)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain3 <- BMFVAR.novol(y, K = K, p = p, dist = "Student", y0 = y0, prior = prior, inits = inits)

prior <- get_prior(y, p = p, dist="Student", SV = F, 
                   aggregation = "triangular", idq = c(4), r = 2)
inits <- get_init(prior)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain4 <- BMFVAR.FHV(y, K = K, p = p, dist = "Student", y0 = y0, prior = prior, inits = inits)

cbind(apply(get_post(Chain3$mcmc$param, element = "B"), 2, mean),
     apply(get_post(Chain4$mcmc$param, element = "B"), 2, mean))

plot(apply((Chain3$mcmc$y_miss), 2, mean),
      apply(Chain4$mcmc$y_miss, 2, mean))

cbind(apply(get_post(Chain3$mcmc$param, element = "nu"), 2, mean),
      apply(get_post(Chain4$mcmc$param, element = "nu"), 2, mean))
cbind(apply(get_post(Chain3$mcmc$param, element = "nu"), 2, quantile, prob = c(0.95)),
      apply(get_post(Chain4$mcmc$param, element = "nu"), 2, quantile, prob = c(0.95)))

apply(get_post(Chain3$mcmc$param, element = "a"), 2, mean)
apply(get_post(Chain4$mcmc$param, element = "L"), 2, mean)

################
library(MFVART)

set.seed(0)

datagen <- sim.MFVAR.novol(dist="OT", idq = c(4), p = 2, K = 4, t_max = 500, aggregation = "triangular")
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

prior <- get_prior(y, p = p, dist="OT", SV = F, aggregation = "triangular", idq = c(4))
inits <- get_init(prior)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain5 <- BMFVAR.novol(y, K = K, p = p, dist = "OT", y0 = y0, prior = prior, inits = inits)

prior <- get_prior(y, p = p, dist="OT", SV = F, 
                   aggregation = "triangular", idq = c(4), r = 2)
inits <- get_init(prior)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain6 <- BMFVAR.FHV(y, K = K, p = p, dist = "OT", y0 = y0, prior = prior, inits = inits)

plot(apply(get_post(Chain5$mcmc$param, element = "B"), 2, mean),
      apply(get_post(Chain6$mcmc$param, element = "B"), 2, mean))

plot(apply((Chain5$mcmc$y_miss), 2, mean),
     apply(Chain6$mcmc$y_miss, 2, mean))

cbind(apply(get_post(Chain5$mcmc$param, element = "nu"), 2, mean),
      apply(get_post(Chain6$mcmc$param, element = "nu"), 2, mean))
cbind(apply(get_post(Chain5$mcmc$param, element = "nu"), 2, quantile, prob = c(0.95)),
      apply(get_post(Chain6$mcmc$param, element = "nu"), 2, quantile, prob = c(0.95)))

apply(get_post(Chain5$mcmc$param, element = "a"), 2, mean)
apply(get_post(Chain6$mcmc$param, element = "L"), 2, mean)



#############################################
library(MFVART)
set.seed(0)

datagen <- sim.MFVAR.SV(dist="Gaussian", idq = c(4), p = 2, K = 4, t_max = 500, aggregation = "triangular")

head(datagen$y)
head(datagen$y0)
tail(datagen$y)
y = datagen$y
y0 = datagen$y0
p = datagen$p
K = datagen$K
t_max <- nrow(y)

prior <- get_prior(y, p = p, dist="Gaussian", SV = T, aggregation = "triangular", idq = c(4))
inits <- get_init(prior)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain7 <- BMFVAR.SV(y, K = K, p = p, dist = "Gaussian", y0 = y0, prior = prior, inits = inits)

prior <- get_prior(y, p = p, dist="Gaussian", SV = T, 
                   aggregation = "triangular", idq = c(4), r = 2)
inits <- get_init(prior)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain8 <- BMFVAR.FSV(y, K = K, p = p, dist = "Gaussian", y0 = y0, prior = prior, inits = inits)

plot(apply(get_post(Chain7$mcmc$param, element = "B"), 2, mean),
      apply(get_post(Chain8$mcmc$param, element = "B"), 2, mean))

plot(apply((Chain7$mcmc$y_miss), 2, mean),
     apply(Chain8$mcmc$y_miss, 2, mean))


#####################
library(MFVART)
set.seed(0)
datagen <- sim.MFVAR.SV(dist="Student", idq = c(4), p = 2, K = 4, t_max = 500, aggregation = "triangular")

head(datagen$y)
head(datagen$y0)
tail(datagen$y)
y = datagen$y
y0 = datagen$y0
p = datagen$p
K = datagen$K
t_max <- nrow(datagen$y)

prior <- get_prior(y, p = p, dist="Student", SV = T, aggregation = "triangular", idq = c(4))
inits <- get_init(prior)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain9 <- BMFVAR.SV(y, K = K, p = p, dist = "Student", y0 = y0, prior = prior, inits = inits)

prior <- get_prior(y, p = p, dist="Student", SV = T, 
                   aggregation = "triangular", idq = c(4), r = 2)
inits <- get_init(prior)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain10 <- BMFVAR.FSV(y, K = K, p = p, dist = "Student", y0 = y0, prior = prior, inits = inits)

plot(apply(get_post(Chain9$mcmc$param, element = "B"), 2, mean),
     apply(get_post(Chain10$mcmc$param, element = "B"), 2, mean))

plot(apply((Chain9$mcmc$y_miss), 2, mean),
     apply(Chain10$mcmc$y_miss, 2, mean))
apply(get_post(Chain9$mcmc$param, element = "nu"), 2, mean)
apply(Chain10$mcmc$param, 2, mean)[1:50] # 5.4874

#####################
library(MFVART)
set.seed(0)
datagen <- sim.MFVAR.SV(dist="OT", idq = c(4), p = 2, K = 4, t_max = 500, aggregation = "triangular")

head(datagen$y)
head(datagen$y0)
tail(datagen$y)
y = datagen$y
y0 = datagen$y0
p = datagen$p
K = datagen$K
t_max <- nrow(datagen$y)

prior <- get_prior(y, p = p, dist="OT", SV = T, aggregation = "triangular", idq = c(4))
inits <- get_init(prior)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain11 <- BMFVAR.SV(y, K = K, p = p, dist = "OT", y0 = y0, prior = prior, inits = inits)

prior <- get_prior(y, p = p, dist="OT", SV = T, 
                   aggregation = "triangular", idq = c(4), r = 2)
inits <- get_init(prior)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain12 <- BMFVAR.FSV(y, K = K, p = p, dist = "OT", y0 = y0, prior = prior, inits = inits)

plot(apply(get_post(Chain11$mcmc$param, element = "B"), 2, mean),
     apply(get_post(Chain12$mcmc$param, element = "B"), 2, mean))

plot(apply((Chain11$mcmc$y_miss), 2, mean),
     apply(Chain12$mcmc$y_miss, 2, mean))
apply(get_post(Chain11$mcmc$param, element = "nu"), 2, mean)
apply(get_post(Chain12$mcmc$param, element = "nu"), 2, mean)

y[,4] <- NA
yt <- t(datagen$y_true)
plot( yt[is.na(t(y))], 
      apply(Chain12$mcmc$y_miss, MARGIN = 2, FUN = mean))
abline(a = 0, b = 1)
