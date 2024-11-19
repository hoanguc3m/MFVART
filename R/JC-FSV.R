#' Bayesian inference of VAR model with factor RW-SV
#'
#' Bayesian inference of VAR model with factor RW-SV
#' \deqn{y_t = B x_t + L f_t + SQRT(w_t) H_t eps_t}
#' @param y The input data as matrix T x K where T is the number of observations and K is the number of variables
#' @param K The number of variables in JCVARF model.
#' @param p The number of lags in JCVARF model.
#' @param dist The variable specifies the JCVARF error distribution. It should be one of
#' c("Gaussian","Student","OT").
#' @param y0 The number of observations
#' @param prior The prior specification of JCVARF.
#' @param inits The initial values of JCVARF.
#' @return A coda object of the posterior samples.
#' @export
#' @examples
#' \dontrun{
#' datagen <- sim.VAR.SV(dist="Gaussian")
#' y <- datagen$y
#' prior <- get_prior(y, p = 2, dist="Gaussian", SV = T)
#' inits <- get_init(prior)
#' Chain1 <- JCVARF.SV(y, K = 5, p = 2, dist = "Gaussian", y0 = NULL, prior = prior, inits = inits)
#' plot(Chain1)
#' }
JCVARF.SV <- function(y, K, p, dist, y0 = NULL, prior = NULL, inits = NULL){
  if (!(dist %in% c("Gaussian","Student","OT") ))
    stop("dist is not implemented.")
  if (prior$SV == TRUE){
    Start = Sys.time()
    if (dist == "Gaussian") Chain <- JCVARF.Gaussian.SV(y, K, p, y0, prior, inits)
    if (dist == "Student") Chain <- JCVARF.Student.SV(y, K, p, y0, prior, inits)
    if (dist == "OT") Chain <- JCVARF.OT.SV(y, K, p, y0, prior, inits)
    elapsedTime = Sys.time() - Start
    print(elapsedTime)
    out <- list(mcmc = Chain,
                y = y,
                y0 = y0,
                K = K,
                p = p,
                dist = dist,
                prior = prior,
                inits = inits,
                esttime = elapsedTime)
    class(out) <- c("JCVARF")
    return(out)
  } else {
    warning("prior$SV is TRUE")
  }
}
#' @export
JCVARF.Gaussian.SV <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL){
  # Init regressors in the right hand side
  t_max <- nrow(y)
  yt = t(y)
  xt <- makeRegressor(y, y0, t_max, K, p)
  
  # Init prior and initial values
  m = K * p + 1
  if (is.null(prior)){
    prior <- get_prior(y, p, priorStyle = "Minnesota", dist = "Gaussian", SV = TRUE)
  }
  # prior B
  b_prior = matrix(prior$b_prior, nrow = K)
  V_b_prior = matrix(diag(prior$V_b_prior), nrow = K)
  # prior sigma
  sigma0_T0 <- prior$sigma_T0
  sigma0_S0 <- prior$sigma_S0
  # prior h0_mean
  h0_mean <- 2 * log(prior$sigma)
  
  # Initial values
  if (is.null(inits)){
    inits <- get_init(prior)
  }
  samples <- inits$samples
  B <- inits$B0
  h <- inits$h
  sigma_h <- inits$sigma_h
  
  # new precompute
  r <- prior$r
  Vlhyper <- prior$Vlhyper
  mulhyper <- prior$mulhyper
  
  sigma_hf <- diag(rep(1, r)) # Variance of factor
  n_load <- K*r-r*(r+1)/2
  L_idx <- c(1:(K*r))[lower.tri(matrix(1, nrow = K, ncol = r))]
  L <- inits$A0[,1:r]
  
  hf <- matrix(0, nrow = r, ncol = t_max)
  hf0 <- rep(0,r)
  # Output
  mcmc <- matrix(NA, nrow = m*K + length(L_idx) + K+r + K + K*t_max,
                 ncol = (samples - inits$burnin)%/% inits$thin)
  for (j in c(1:samples)){
    # Sample f
    e = matrix((yt - B %*%xt),ncol = 1)
    Xf <- kronecker(Diagonal(t_max), L)
    XfiSig <- Matrix::t(Xf) %*% Diagonal(t_max*K, as.vector(exp(-h)) )
    Kf <- Diagonal(t_max*r, as.vector(exp(-hf)) ) + XfiSig %*% Xf
    CKf <- chol(Kf)
    f <- backsolve( CKf,
                    backsolve( CKf, XfiSig %*% e,
                               upper.tri = T, transpose = T )
                    + rnorm(t_max*r) )
    # f_hat = backsolve(Kf, XfiSig %*% e)
    # f = f_hat +  backsolve(chol(Kf),  rnorm(t_max*r) )
    Fmat <- matrix(f, r, t_max)
    
    # sample L and B
    for (ii in 1:K) {
      if (ii <= r) {
        if (ii ==1){
          Zi <- xt
        } else {
          Zi <- rbind(xt, Fmat[1:(ii-1),, drop = FALSE])  
        }
        
        iVthetai <- Diagonal(x = 1/c(V_b_prior[ii,], 
                                     rep(1 / Vlhyper, ii-1)))
        thetai0 <- c(b_prior[ii,], rep(mulhyper, ii-1))
        ZiSigi <- Zi %*% Diagonal(t_max, exp(-h[ii,]))
        dthetai <- iVthetai %*% thetai0 + ZiSigi %*% (yt[ii,] - Fmat[ii,])
      } else {
        Zi <- rbind(xt, Fmat)
        iVthetai <- Diagonal(x = 1/c(V_b_prior[ii,], rep(Vlhyper,r)))
        thetai0 <- c(b_prior[ii,], rep(mulhyper, r))
        ZiSigi <- Zi %*% Diagonal(t_max, exp(-h[ii,]))
        dthetai <- iVthetai %*% thetai0 + ZiSigi %*% yt[ii,]
      }
      
      Kthetai <- iVthetai + ZiSigi %*% Matrix::t(Zi)
      CKthetai <- chol(Kthetai)
      
      thetai <- backsolve( CKthetai,
                           backsolve( CKthetai, dthetai,
                                      upper.tri = T, transpose = T )
                           + rnorm(m+min(ii-1, r)) )
      B[ii, ] <- thetai[1:m]
      if (ii > 1){
        L[ii, 1:min(ii-1, r)] <- thetai[(m+1):length(thetai)]  
      }
      
    }
    
    
    # Sample vol
    ytilde <- (yt - B %*% xt - L %*% Fmat)
    aux <- sample_h_ele(ytilde = ytilde, sigma_h = sigma_h, h0_mean = h0_mean, h = h, K = K, t_max = t_max, prior = prior)
    h <- aux$Sigtdraw
    h0 <- as.numeric(aux$h0)
    sqrtvol <- aux$sigt
    sigma_h <- aux$sigma_h
    
    aux <- sample_h_ele(ytilde = Fmat, sigma_h = sigma_hf, h0_mean = rep(0,r), 
                        h = hf, K = r, t_max = t_max, prior = prior)
    hf <- aux$Sigtdraw
    hf0 <- as.numeric(aux$h0)
    sqrtvolf <- aux$sigt
    sigma_hf <- aux$sigma_h
    
    
    if ((j > inits$burnin) & (j %% inits$thin == 0)){
      mcmc[, (j - inits$burnin) %/% inits$thin] <- c(vec(B), L[L_idx], diag(sigma_h), diag(sigma_hf), h0, as.numeric(h))
    }
    
    if (j %% 1000 == 0) { cat(" Iteration ", j, " \n")}
  }
  nameL <- matrix(paste("L", repcol(c(1:K),r), reprow(c(1:r),K), sep = "_"), ncol = r)
  nameL <- nameL[L_idx]
  row.names(mcmc) <- c( paste("B0",c(1:K), sep = ""),
                        sprintf("B%d_%d_%d",reprow(c(1:p),K*K), rep(repcol(c(1:K),K), p), rep(reprow(c(1:K),K), p)),
                        nameL,
                        paste("sigma_h",c(1:K), sep = ""),
                        paste("sigma_hf",c(1:r), sep = ""),
                        paste("lh0",c(1:K), sep = ""),
                        sprintf("h_%d_%d", repcol(c(1:K),t_max), reprow(c(1:t_max),K))
  )
  return(as.mcmc(t(mcmc)))
}

###########################################################################
#' @export
JCVARF.Student.SV <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL){
  # Init regressors in the right hand side
  t_max <- nrow(y)
  yt = t(y)
  xt <- makeRegressor(y, y0, t_max, K, p)
  
  # Init prior and initial values
  m = K * p + 1
  if (is.null(prior)){
    prior <- get_prior(y, p, priorStyle = "Minnesota", dist = "Student", SV = TRUE)
  }
  # prior B
  b_prior = matrix(prior$b_prior, nrow = K)
  V_b_prior = matrix(diag(prior$V_b_prior), nrow = K)
  # prior sigma
  sigma0_T0 <- prior$sigma_T0
  sigma0_S0 <- prior$sigma_S0
  # prior h0_mean
  h0_mean <- 2 * log(prior$sigma)
  # prior nu
  nu_gam_a = prior$nu_gam_a
  nu_gam_b = prior$nu_gam_b
  
  # Initial values
  if (is.null(inits)){
    inits <- get_init(prior)
  }
  samples <- inits$samples
  B <- inits$B0
  h <- inits$h
  sigma_h <- inits$sigma_h
  
  # Multi degrees of freedom
  nu <- inits$nu
  logsigma_nu <- 0
  acount_nu <- 0
  acount_w <- rep(0, t_max)
  # Init w as Gaussian
  w_sample <- rep(1, t_max)
  w <- reprow(w_sample, K)
  w_sqrt <- sqrt(w)
  
  # new precompute
  r <- prior$r
  Vlhyper <- prior$Vlhyper
  mulhyper <- prior$mulhyper
  
  sigma_hf <- diag(rep(1, r)) # Variance of factor
  n_load <- K*r-r*(r+1)/2
  L_idx <- c(1:(K*r))[lower.tri(matrix(1, nrow = K, ncol = r))]
  L <- inits$A0[,1:r]
  
  hf <- matrix(0, nrow = r, ncol = t_max)
  hf0 <- rep(0,r)
  
  # Output
  mcmc <-  matrix(NA, nrow = m*K + length(L_idx) + 1 + K+r + K + K*t_max + t_max,
                  ncol = (samples - inits$burnin)%/% inits$thin)
  for (j in c(1:samples)){
    # Sample f
    e = matrix((yt - B %*%xt),ncol = 1)
    Xf <- kronecker(Diagonal(t_max), L)
    XfiSig <- Matrix::t(Xf) %*% Diagonal(t_max*K, as.vector(exp(-h))/as.vector(w) )
    Kf <- Diagonal(t_max*r, as.vector(exp(-hf)) ) + XfiSig %*% Xf
    CKf <- chol(Kf)
    f <- backsolve( CKf,
                    backsolve( CKf, XfiSig %*% e,
                               upper.tri = T, transpose = T )
                    + rnorm(t_max*r) )
    # f_hat = backsolve(Kf, XfiSig %*% e)
    # f = f_hat +  backsolve(chol(Kf),  rnorm(t_max*r) )
    Fmat <- matrix(f, r, t_max)
    
    # sample L and B
    for (ii in 1:K) {
      if (ii <= r) {
        if (ii ==1){
          Zi <- xt
        } else {
          Zi <- rbind(xt, Fmat[1:(ii-1),, drop = FALSE])  
        }
        
        iVthetai <- Diagonal(x = 1/c(V_b_prior[ii,], 
                                     rep(1 / Vlhyper, ii-1)))
        thetai0 <- c(b_prior[ii,], rep(mulhyper, ii-1))
        ZiSigi <- Zi %*% Diagonal(t_max, exp(-h[ii,]) / w[ii,] )
        dthetai <- iVthetai %*% thetai0 + ZiSigi %*% (yt[ii,] - Fmat[ii,])
      } else {
        Zi <- rbind(xt, Fmat)
        iVthetai <- Diagonal(x = 1/c(V_b_prior[ii,], rep(Vlhyper,r)))
        thetai0 <- c(b_prior[ii,], rep(mulhyper, r))
        ZiSigi <- Zi %*% Diagonal(t_max, exp(-h[ii,]) / w[ii,])
        dthetai <- iVthetai %*% thetai0 + ZiSigi %*% yt[ii,]
      }
      
      Kthetai <- iVthetai + ZiSigi %*% Matrix::t(Zi)
      CKthetai <- chol(Kthetai)
      
      thetai <- backsolve( CKthetai,
                           backsolve( CKthetai, dthetai,
                                      upper.tri = T, transpose = T )
                           + rnorm(m+min(ii-1, r)) )
      B[ii, ] <- thetai[1:m]
      if (ii > 1){
        L[ii, 1:min(ii-1, r)] <- thetai[(m+1):length(thetai)]  
      }
      
    }
    
    
    
    # Sample vol
    ytilde <- (yt - B %*% xt - L %*% Fmat) / w_sqrt
    aux <- sample_h_ele(ytilde = ytilde, sigma_h = sigma_h, h0_mean = h0_mean, h = h, K = K, t_max = t_max, prior = prior)
    h <- aux$Sigtdraw
    h0 <- as.numeric(aux$h0)
    sqrtvol <- aux$sigt
    sigma_h <- aux$sigma_h
    
    aux <- sample_h_ele(ytilde = Fmat, sigma_h = sigma_hf, h0_mean = rep(0,r), 
                        h = hf, K = r, t_max = t_max, prior = prior)
    hf <- aux$Sigtdraw
    hf0 <- as.numeric(aux$h0)
    sqrtvolf <- aux$sigt
    sigma_hf <- aux$sigma_h  
    
    # Sample w
    u <-  (yt - B %*% xt - L %*% Fmat) 
    w_sample <- rinvgamma( n = t_max, shape = (nu+K)*0.5, rate = 0.5*( nu + colSums((u^2)/exp(h)) ) )
    w <- reprow(w_sample, K)
    w_sqrt <- sqrt(w)
    
    # Sample nu
    nu_temp = nu + exp(logsigma_nu)*rnorm(1)
    if (nu_temp > 4 && nu_temp < 100){
      num_mh = dgamma(nu_temp, shape = nu_gam_a, rate = nu_gam_b, log = T) +
        sum(dinvgamma(w_sample, shape = nu_temp*0.5, rate = nu_temp*0.5, log = T))
      denum_mh = dgamma(nu, shape = nu_gam_a, rate = nu_gam_b, log = T) +
        sum(dinvgamma(w_sample, shape = nu*0.5, rate = nu*0.5, log = T))
      alpha = num_mh - denum_mh;
      temp = log(runif(1));
      if (alpha > temp){
        nu = nu_temp
        acount_nu = acount_nu + 1
      }
      
    }
    
    if(j %% batchlength == 0 ){
      if (acount_nu > batchlength * TARGACCEPT){
        logsigma_nu = logsigma_nu + adaptamount(j %/% batchlength);
      }
      if (acount_nu < batchlength * TARGACCEPT){
        logsigma_nu = logsigma_nu - adaptamount(j %/% batchlength);
      }
      acount_nu = 0
    }
    
    if ((j > inits$burnin) & (j %% inits$thin == 0))
      mcmc[, (j - inits$burnin) %/% inits$thin] <- c(vec(B), L[L_idx], nu, diag(sigma_h), diag(sigma_hf), h0, as.numeric(h), as.numeric(w_sample))
    if (j %% 1000 == 0) {
      cat(" Iteration ", j, " ", logsigma_nu," ", round(nu,2)," \n")
      acount_w <- rep(0,t_max)
    }
  }
  nameL <- matrix(paste("L", repcol(c(1:K),r), reprow(c(1:r),K), sep = "_"), ncol = r)
  nameL <- nameL[L_idx]
  row.names(mcmc) <- c( paste("B0",c(1:K), sep = ""),
                        sprintf("B%d_%d_%d",reprow(c(1:p),K*K), rep(repcol(c(1:K),K), p), rep(reprow(c(1:K),K), p)),
                        nameL,
                        paste("nu"),
                        paste("sigma_h",c(1:K), sep = ""),
                        paste("sigma_hf",c(1:r), sep = ""),
                        paste("lh0",c(1:K), sep = ""),
                        sprintf("h_%d_%d", repcol(c(1:K),t_max), reprow(c(1:t_max),K)),
                        paste("w",c(1:t_max), sep = ""))
  
  return(as.mcmc(t(mcmc)))
}


###########################################################################
#' @export
JCVARF.OT.SV <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL){
  # Init regressors in the right hand side
  t_max <- nrow(y)
  yt = t(y)
  xt <- makeRegressor(y, y0, t_max, K, p)
  
  # Init prior and initial values
  m = K * p + 1
  if (is.null(prior)){
    prior <- get_prior(y, p, priorStyle = "Minnesota", dist = "OT", SV = TRUE)
  }
  # prior B
  b_prior = matrix(prior$b_prior, nrow = K)
  V_b_prior = matrix(diag(prior$V_b_prior), nrow = K)
  # prior sigma
  sigma0_T0 <- prior$sigma_T0
  sigma0_S0 <- prior$sigma_S0
  # prior h0_mean
  h0_mean <- 2 * log(prior$sigma)
  # prior nu
  nu_gam_a = prior$nu_gam_a
  nu_gam_b = prior$nu_gam_b
  
  # Initial values
  if (is.null(inits)){
    inits <- get_init(prior)
  }
  samples <- inits$samples
  B <- inits$B0
  h <- inits$h
  sigma_h <- inits$sigma_h
  
  # Multi degrees of freedom
  nu <- inits$nu
  logsigma_nu <- rep(0,K)
  acount_nu <- rep(0,K)
  acount_w <- rep(0, t_max)
  # Init w as Gaussian
  w_sample <- rep(1, t_max)
  w <- reprow(w_sample, K)
  w_sqrt <- sqrt(w)
  
  # new precompute
  r <- prior$r
  Vlhyper <- prior$Vlhyper
  mulhyper <- prior$mulhyper
  
  sigma_hf <- diag(rep(1, r)) # Variance of factor
  n_load <- K*r-r*(r+1)/2
  L_idx <- c(1:(K*r))[lower.tri(matrix(1, nrow = K, ncol = r))]
  L <- inits$A0[,1:r]
  
  hf <- matrix(0, nrow = r, ncol = t_max)
  hf0 <- rep(0,r)
  
  # Output
  mcmc <- matrix(NA, nrow = m*K + length(L_idx) + K + K+r + K + K*t_max + K*t_max,
                 ncol = (samples - inits$burnin)%/% inits$thin)
  for (j in c(1:samples)){
    # Sample f
    e = matrix((yt - B %*%xt),ncol = 1)
    Xf <- kronecker(Diagonal(t_max), L)
    XfiSig <- Matrix::t(Xf) %*% Diagonal(t_max*K, as.vector(exp(-h))/as.vector(w) )
    Kf <- Diagonal(t_max*r, as.vector(exp(-hf)) ) + XfiSig %*% Xf
    CKf <- chol(Kf)
    f <- backsolve( CKf,
                    backsolve( CKf, XfiSig %*% e,
                               upper.tri = T, transpose = T )
                    + rnorm(t_max*r) )
    # f_hat = backsolve(Kf, XfiSig %*% e)
    # f = f_hat +  backsolve(chol(Kf),  rnorm(t_max*r) )
    Fmat <- matrix(f, r, t_max)
    
    # sample L and B
    for (ii in 1:K) {
      if (ii <= r) {
        if (ii ==1){
          Zi <- xt
        } else {
          Zi <- rbind(xt, Fmat[1:(ii-1),, drop = FALSE])  
        }
        
        iVthetai <- Diagonal(x = 1/c(V_b_prior[ii,], 
                                     rep(1 / Vlhyper, ii-1)))
        thetai0 <- c(b_prior[ii,], rep(mulhyper, ii-1))
        ZiSigi <- Zi %*% Diagonal(t_max, exp(-h[ii,]) / w[ii,] )
        dthetai <- iVthetai %*% thetai0 + ZiSigi %*% (yt[ii,] - Fmat[ii,])
      } else {
        Zi <- rbind(xt, Fmat)
        iVthetai <- Diagonal(x = 1/c(V_b_prior[ii,], rep(Vlhyper,r)))
        thetai0 <- c(b_prior[ii,], rep(mulhyper, r))
        ZiSigi <- Zi %*% Diagonal(t_max, exp(-h[ii,]) / w[ii,])
        dthetai <- iVthetai %*% thetai0 + ZiSigi %*% yt[ii,]
      }
      
      Kthetai <- iVthetai + ZiSigi %*% Matrix::t(Zi)
      CKthetai <- chol(Kthetai)
      
      thetai <- backsolve( CKthetai,
                           backsolve( CKthetai, dthetai,
                                      upper.tri = T, transpose = T )
                           + rnorm(m+min(ii-1, r)) )
      B[ii, ] <- thetai[1:m]
      if (ii > 1){
        L[ii, 1:min(ii-1, r)] <- thetai[(m+1):length(thetai)]  
      }
      
    }
    
    
    
    # Sample vol
    ytilde <- (yt - B %*% xt - L %*% Fmat) / w_sqrt
    aux <- sample_h_ele(ytilde = ytilde, sigma_h = sigma_h, h0_mean = h0_mean, h = h, K = K, t_max = t_max, prior = prior)
    # aux <- sample_h_Chan(ytilde = ytilde, sigma_h = sigma_h, h0_mean = h0_mean, h = h, K = K, t_max = t_max, prior = prior)
    h <- aux$Sigtdraw
    h0 <- as.numeric(aux$h0)
    sqrtvol <- aux$sigt
    sigma_h <- aux$sigma_h
    
    aux <- sample_h_ele(ytilde = Fmat, sigma_h = sigma_hf, h0_mean = rep(0,r),
                        h = hf, K = r, t_max = t_max, prior = prior)
    # aux <- sample_h_Chan(ytilde = Fmat, sigma_h = sigma_hf, h0_mean = rep(0,r), 
    #                     h = hf, K = r, t_max = t_max, prior = prior)
    hf <- aux$Sigtdraw
    hf0 <- as.numeric(aux$h0)
    sqrtvolf <- aux$sigt
    sigma_hf <- aux$sigma_h  
    
    # Sample w
    u <-  (yt - B %*% xt - L %*% Fmat)
    w <- matrix( rinvgamma( n = K*t_max, shape = nu*0.5 + 1*0.5, rate = nu*0.5 + 0.5 * (u^2) / exp(h) ),
                 K, t_max )
    w_sqrt <- sqrt(w)
    
    
    # Sample nu
    nu_temp = nu + exp(logsigma_nu)*rnorm(K)
    for (k in c(1:K)){
      if (nu_temp[k] > 4 && nu_temp[k] < 100){
        num_mh = dgamma(nu_temp[k], shape = nu_gam_a, rate = nu_gam_b, log = T) +
          sum(dinvgamma(w[k,], shape = nu_temp[k]*0.5, rate = nu_temp[k]*0.5, log = T))
        denum_mh = dgamma(nu[k], shape = nu_gam_a, rate = nu_gam_b, log = T) +
          sum(dinvgamma(w[k,], shape = nu[k]*0.5, rate = nu[k]*0.5, log = T))
        alpha = num_mh - denum_mh;
        temp = log(runif(1));
        if (alpha > temp){
          nu[k] = nu_temp[k]
          acount_nu[k] = acount_nu[k] + 1
        }
        
      }
    }
    
    
    if(j %% batchlength == 0 ){
      for (jj in c(1:K)) {
        if (acount_nu[jj] > batchlength * TARGACCEPT){
          logsigma_nu[jj] = logsigma_nu[jj] + adaptamount(j %/% batchlength);
        }
        if (acount_nu[jj] < batchlength * TARGACCEPT){
          logsigma_nu[jj] = logsigma_nu[jj] - adaptamount(j %/% batchlength);
        }
        acount_nu[jj] = 0
      }
    }
    
    if ((j > inits$burnin) & (j %% inits$thin == 0))
      mcmc[, (j - inits$burnin) %/% inits$thin] <- c(vec(B), L[L_idx], nu, diag(sigma_h), diag(sigma_hf), h0, as.numeric(h), as.numeric(w))
    
    if (j %% 1000 == 0) {
      cat(" Iteration ", j, " ", logsigma_nu," ", round(nu,2)," \n")
      acount_w <- rep(0,t_max)
    }
  }
  nameL <- matrix(paste("L", repcol(c(1:K),r), reprow(c(1:r),K), sep = "_"), ncol = r)
  nameL <- nameL[L_idx]
  row.names(mcmc) <- c( paste("B0",c(1:K), sep = ""),
                        sprintf("B%d_%d_%d",reprow(c(1:p),K*K), rep(repcol(c(1:K),K), p), rep(reprow(c(1:K),K), p)),
                        nameL,
                        paste("nu",c(1:K), sep = ""),
                        paste("sigma_h",c(1:K), sep = ""),
                        paste("sigma_hf",c(1:r), sep = ""),
                        paste("lh0",c(1:K), sep = ""),
                        sprintf("h_%d_%d", repcol(c(1:K),t_max), reprow(c(1:t_max),K)),
                        sprintf("w_%d_%d", repcol(c(1:K),t_max), reprow(c(1:t_max),K)))
  
  return(as.mcmc(t(mcmc)))
}

#' @export
JCVARF.OT.SVinR <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL){
  # Init regressors in the right hand side
  t_max <- nrow(y)
  yt = t(y)
  xt <- makeRegressor(y, y0, t_max, K, p)
  
  # Init prior and initial values
  m = K * p + 1
  if (is.null(prior)){
    prior <- get_prior(y, p, priorStyle = "Minnesota", dist = "OT", SV = TRUE)
  }
  # prior B
  b_prior = matrix(prior$b_prior, nrow = K)
  V_b_prior = matrix(diag(prior$V_b_prior), nrow = K)
  # prior sigma
  sigma0_T0 <- prior$sigma_T0
  sigma0_S0 <- prior$sigma_S0
  # prior h0_mean
  h0_mean <- 2 * log(prior$sigma)
  # prior nu
  nu_gam_a = prior$nu_gam_a
  nu_gam_b = prior$nu_gam_b
  
  # Initial values
  if (is.null(inits)){
    inits <- get_init(prior)
  }
  samples <- inits$samples
  B <- inits$B0
  h <- inits$h
  sigma_h <- inits$sigma_h
  
  # Multi degrees of freedom
  nu <- inits$nu
  logsigma_nu <- rep(0,K)
  acount_nu <- rep(0,K)
  acount_w <- rep(0, t_max)
  # Init w as Gaussian
  w_sample <- rep(1, t_max)
  w <- reprow(w_sample, K)
  w_sqrt <- sqrt(w)
  
  # new precompute
  r <- prior$r
  Vlhyper <- prior$Vlhyper
  mulhyper <- prior$mulhyper
  
  sigma_hf <- diag(rep(1, r)) # Variance of factor
  n_load <- K*r-r*(r+1)/2
  L_idx <- c(1:(K*r))[lower.tri(matrix(1, nrow = K, ncol = r))]
  L <- inits$A0[,1:r]
  
  hf <- matrix(0, nrow = r, ncol = t_max)
  hf0 <- rep(0,r)
  
  # Output
  mcmc <- matrix(NA, nrow = m*K + length(L_idx) + K + K+r + K + K*t_max + K*t_max,
                 ncol = (samples - inits$burnin)%/% inits$thin)
  for (j in c(1:samples)){
    # Sample f
    e = matrix((yt - B %*%xt),ncol = 1)
    Xf <- kronecker(Diagonal(t_max), L)
    XfiSig <- Matrix::t(Xf) %*% Diagonal(t_max*K, as.vector(exp(-h))/as.vector(w) )
    Kf <- Diagonal(t_max*r, as.vector(exp(-hf)) ) + XfiSig %*% Xf
    CKf <- chol(Kf)
    f <- backsolve( CKf,
                    backsolve( CKf, XfiSig %*% e,
                               upper.tri = T, transpose = T )
                    + rnorm(t_max*r) )
    # f_hat = backsolve(Kf, XfiSig %*% e)
    # f = f_hat +  backsolve(chol(Kf),  rnorm(t_max*r) )
    Fmat <- matrix(f, r, t_max)
    
    # sample L and B
    for (ii in 1:K) {
      if (ii <= r) {
        if (ii ==1){
          Zi <- xt
        } else {
          Zi <- rbind(xt, Fmat[1:(ii-1),, drop = FALSE])  
        }
        
        iVthetai <- Diagonal(x = 1/c(V_b_prior[ii,], 
                                     rep(1 / Vlhyper, ii-1)))
        thetai0 <- c(b_prior[ii,], rep(mulhyper, ii-1))
        ZiSigi <- Zi %*% Diagonal(t_max, exp(-h[ii,]) / w[ii,] )
        dthetai <- iVthetai %*% thetai0 + ZiSigi %*% (yt[ii,] - Fmat[ii,])
      } else {
        Zi <- rbind(xt, Fmat)
        iVthetai <- Diagonal(x = 1/c(V_b_prior[ii,], rep(Vlhyper,r)))
        thetai0 <- c(b_prior[ii,], rep(mulhyper, r))
        ZiSigi <- Zi %*% Diagonal(t_max, exp(-h[ii,]) / w[ii,])
        dthetai <- iVthetai %*% thetai0 + ZiSigi %*% yt[ii,]
      }
      
      Kthetai <- iVthetai + ZiSigi %*% Matrix::t(Zi)
      CKthetai <- chol(Kthetai)
      
      thetai <- backsolve( CKthetai,
                           backsolve( CKthetai, dthetai,
                                      upper.tri = T, transpose = T )
                           + rnorm(m+min(ii-1, r)) )
      B[ii, ] <- thetai[1:m]
      if (ii > 1){
        L[ii, 1:min(ii-1, r)] <- thetai[(m+1):length(thetai)]  
      }
      
    }
    
    
    
    # Sample vol
    ytilde <- (yt - B %*% xt - L %*% Fmat) / w_sqrt
    # aux <- sample_h_ele(ytilde = ytilde, sigma_h = sigma_h, h0_mean = h0_mean, h = h, K = K, t_max = t_max, prior = prior)
    aux <- sample_h_Chan(ytilde = ytilde, sigma_h = sigma_h, h0_mean = h0_mean, h = h, K = K, t_max = t_max, prior = prior)
    h <- aux$Sigtdraw
    h0 <- as.numeric(aux$h0)
    sqrtvol <- aux$sigt
    sigma_h <- aux$sigma_h
    
    # aux <- sample_h_ele(ytilde = Fmat, sigma_h = sigma_hf, h0_mean = rep(0,r), 
    #                     h = hf, K = r, t_max = t_max, prior = prior)
    aux <- sample_h_Chan(ytilde = Fmat, sigma_h = sigma_hf, h0_mean = rep(0,r), 
                         h = hf, K = r, t_max = t_max, prior = prior)
    hf <- aux$Sigtdraw
    hf0 <- as.numeric(aux$h0)
    sqrtvolf <- aux$sigt
    sigma_hf <- aux$sigma_h  
    
    # Sample w
    u <-  (yt - B %*% xt - L %*% Fmat)
    w <- matrix( rinvgamma( n = K*t_max, shape = nu*0.5 + 1*0.5, rate = nu*0.5 + 0.5 * (u^2) / exp(h) ),
                 K, t_max )
    w_sqrt <- sqrt(w)
    
    
    # Sample nu
    nu_temp = nu + exp(logsigma_nu)*rnorm(K)
    for (k in c(1:K)){
      if (nu_temp[k] > 4 && nu_temp[k] < 100){
        num_mh = dgamma(nu_temp[k], shape = nu_gam_a, rate = nu_gam_b, log = T) +
          sum(dinvgamma(w[k,], shape = nu_temp[k]*0.5, rate = nu_temp[k]*0.5, log = T))
        denum_mh = dgamma(nu[k], shape = nu_gam_a, rate = nu_gam_b, log = T) +
          sum(dinvgamma(w[k,], shape = nu[k]*0.5, rate = nu[k]*0.5, log = T))
        alpha = num_mh - denum_mh;
        temp = log(runif(1));
        if (alpha > temp){
          nu[k] = nu_temp[k]
          acount_nu[k] = acount_nu[k] + 1
        }
        
      }
    }
    
    
    if(j %% batchlength == 0 ){
      for (jj in c(1:K)) {
        if (acount_nu[jj] > batchlength * TARGACCEPT){
          logsigma_nu[jj] = logsigma_nu[jj] + adaptamount(j %/% batchlength);
        }
        if (acount_nu[jj] < batchlength * TARGACCEPT){
          logsigma_nu[jj] = logsigma_nu[jj] - adaptamount(j %/% batchlength);
        }
        acount_nu[jj] = 0
      }
    }
    
    if ((j > inits$burnin) & (j %% inits$thin == 0))
      mcmc[, (j - inits$burnin) %/% inits$thin] <- c(vec(B), L[L_idx], nu, diag(sigma_h), diag(sigma_hf), h0, as.numeric(h), as.numeric(w))
    
    if (j %% 1000 == 0) {
      cat(" Iteration ", j, " ", logsigma_nu," ", round(nu,2)," \n")
      acount_w <- rep(0,t_max)
    }
  }
  nameL <- matrix(paste("L", repcol(c(1:K),r), reprow(c(1:r),K), sep = "_"), ncol = r)
  nameL <- nameL[L_idx]
  row.names(mcmc) <- c( paste("B0",c(1:K), sep = ""),
                        sprintf("B%d_%d_%d",reprow(c(1:p),K*K), rep(repcol(c(1:K),K), p), rep(reprow(c(1:K),K), p)),
                        nameL,
                        paste("nu",c(1:K), sep = ""),
                        paste("sigma_h",c(1:K), sep = ""),
                        paste("sigma_hf",c(1:r), sep = ""),
                        paste("lh0",c(1:K), sep = ""),
                        sprintf("h_%d_%d", repcol(c(1:K),t_max), reprow(c(1:t_max),K)),
                        sprintf("w_%d_%d", repcol(c(1:K),t_max), reprow(c(1:t_max),K)))
  
  return(as.mcmc(t(mcmc)))
}

