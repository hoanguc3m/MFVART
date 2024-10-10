datagen <- sim.MFVAR.SV(dist="Gaussian", idq = c(5), p = 3)
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

Make_SectMat <- function(y){
  n = ncol (y)
  t_max = nrow(y)
  
  yt <- t(y)
  yt_vec <- vec(yt)

  n_obs <- sum(!is.na(y))
  n_miss <- sum(is.na(y))
  yobs_vec <- na.exclude(yt_vec)
  ymiss_vec <- matrix(NA, nrow = n_miss, ncol = 1)
  
  nT <- n_obs+n_miss
  i_obs <- j_obs <- rep(0,n_obs)
  i_miss <- j_miss <- rep(0,n_miss)
  c_obs <- c_miss <- 0
  
  for (t in c(1:t_max))
    for (i in c(1:n)){
      if (is.na(yt[i,t])){
        c_miss <- c_miss + 1 
        i_miss[c_miss] <- i + n*(t-1)
        j_miss[c_miss] <- c_miss
      } else{
        c_obs <- c_obs + 1
        i_obs[c_obs] <- i + n*(t-1)
        j_obs[c_obs] <- c_obs
      }
    }
  S_obs <- sparseMatrix(i = i_obs, j = j_obs, x = 1, dims =c(nT,n_obs))
  S_miss <- sparseMatrix(i_miss, j_miss, x = 1, dims =c(nT,n_miss))
  # cbind(yt_vec, S_obs%*% yobs_vec + S_miss%*% ymiss_vec)
}
Make_HB <- function(B){
  b_intercept = B[,1]
  m =  B[,2:ncol(B)]
  B_t = t(matrix(t(m), nrow=nrow(m))[, c(matrix(1:ncol(m), nrow(m), byrow=T)) ])
  IB_t = rbind( diag(ncol = ncol(B_t), nrow = ncol(B_t)), B_t )
}