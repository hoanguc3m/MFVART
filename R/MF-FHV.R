#' Bayesian inference of MF VAR model with factor homogeneous vol
#'
#' Bayesian inference of MF VAR model with factor homogeneous vol
#' \deqn{y_t = B x_t + L f_t + SQRT(w_t)  Sigma eps_t}
#' @param y The input data as matrix T x K where T is the number of observations and K is the number of variables
#' @param K The number of variables in BVAR model.
#' @param p The number of lags in BVAR model.
#' @param r The number of factor in BVAR model used in prior arg.
#' @param dist The variable specifies the BVAR error distribution. It should be one of
#' c("Gaussian","Student","OT").
#' @param y0 The number of observations
#' @param prior The prior specification of BMFVAR.
#' @param inits The initial values of BMFVAR.
#' @return A coda object of the posterior samples.
#' @export
#' @examples
#' \dontrun{
#' datagen <- sim.MFVAR.FHV(dist="Gaussian", idq = c(4), p = 2, K = 4, t_max = 500)
#' y <- datagen$y
#' prior <- get_prior(y, p = 2, dist="Gaussian", SV = F)
#' inits <- get_init(prior)
#' Chain1 <- BMFVAR.FHV(y, K = 5, p = 2, dist = "Gaussian", y0 = NULL, prior = prior, inits = inits)
#' plot(Chain1)
#' }
BMFVAR.FHV <- function(y, K, p, dist, y0 = NULL, prior = NULL, inits = NULL){
  if (!(dist %in% c("Gaussian","Student",
                    "MT","OT") ))
    stop("dist is not implemented.")

  if (prior$SV == FALSE){
    Start = Sys.time()
    if (dist == "Gaussian") Chain <- BMFVAR.Gaussian.FHV(y, K, p, y0, prior, inits)
    if (dist == "Student") Chain <- BMFVAR.Student.FHV(y, K, p, y0, prior, inits)
    if (dist == "MT") Chain <- BMFVAR.MT.FHV(y, K, p, y0, prior, inits)
    if (dist == "OT") Chain <- BMFVAR.OT.FHV(y, K, p, y0, prior, inits)
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
    class(out) <- c("MFVART-Factor")

  } else {
    warning("prior$SV is TRUE")
  }
  return(out)
}

###########################################################################
#' @export
BMFVAR.Gaussian.FHV <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL){
  # Init regressors in the right hand side
  t_max <- nrow(y)
  y_raw <- y
  
  # Init prior and initial values
  m = K * p + 1
  if (is.null(prior)){
    prior <- get_prior(y, p, priorStyle = "Minnesota", dist = "Gaussian", 
                       SV = FALSE, aggregation = "identity", idq = K)
  }
  
  aggregation <- prior$aggregation
  idq <- prior$idq
  sdeviation <- prior$sdeviation
  ymiss_invprior <- 1/prior$ymiss_prior
  
  if (aggregation == "identity"){
    yt_vec <- vec(t(y))
    yobs_vec <- na.exclude(yt_vec)  
  }
  
  if (aggregation == "average" || aggregation == "triangular"){
    y[,idq] <- NA
    
    yt_vec <- vec(t(y))
    yobs_vec <- na.exclude(yt_vec)  
    z_obs_vec <- na.exclude(vec(t(y_raw[,idq])))  
    
    y_temp <- y
    y_temp[,idq] <- 0
    y_temp <- matrix(t(y_temp)[is.na(t(y))], ncol = 1)
    S_select <- Make_SectMat(y_temp)
    # Account for the missing in other columns rather than in low freq columns
    M_aggregate <- kronecker(Make_bandSparse(t_max, aggregation)[!is.na(y_raw[,idq[1] ]),],
                             Diagonal(length(idq))) %*% Matrix::t(S_select$S_obs)
  }
  
  
  # prior B
  b_prior = matrix(prior$b_prior, nrow = K)
  V_b_prior = matrix(diag(prior$V_b_prior), nrow = K)
  # prior sigma
  sigma0_T0 <- prior$sigma_T0
  sigma0_S0 <- prior$sigma_S0

  # Initial values
  if (is.null(inits)){
    inits <- get_init(prior)
  }
  samples <- inits$samples
  B <- inits$B0
  
  sigma <- inits$sigma
  sigma2 <- sigma^2
  
  # new precompute
  r <- prior$r
  Vlhyper <- prior$Vlhyper
  mulhyper <- prior$mulhyper
  
  sigmaf2 <- rep(1, r) # Variance of factor
  n_load <- K*r-r*(r+1)/2
  L_idx <- c(1:(K*r))[lower.tri(matrix(1, nrow = K, ncol = r))]
  L <- inits$A0[,1:r]
  Fmat <- matrix(0, nrow = r, ncol = t_max)
  
  S_select <- Make_SectMat(y)
  n_obs <- sum(!is.na(y))
  n_miss <- sum(is.na(y))
  n_cond <- length(y_raw[,idq]) - sum(is.na(y_raw[,idq]))
  S_obs <- S_select$S_obs
  S_miss <- S_select$S_miss
  K_prior_y_miss <- Matrix::t(M_aggregate) %*% Diagonal(n_cond, 1/sdeviation) %*% M_aggregate + Diagonal(n_miss, ymiss_invprior)
  b_prior_y_miss <- Matrix::t(M_aggregate) %*% Diagonal(n_cond, 1/sdeviation) %*% ( z_obs_vec)
  
  # Output
  mcmc <- matrix(NA, nrow = m*K + length(L_idx) + K+r,
                 ncol = (samples - inits$burnin)%/% inits$thin)
  My_miss <- matrix(NA, nrow = n_miss,
                    ncol = (samples - inits$burnin)%/% inits$thin)
                
  for (j in c(1:samples)){
    
    # # Sample y missing
    {
      H_B <- Make_HB(B, K, p ,t_max)
      c_B <- Make_cB(y0, B, K, p ,t_max, LF = as.numeric(L %*% Fmat) )
      # IA <- kronecker(Diagonal(t_max), Diagonal(K) ) # t_max * K diag
      c_const <- c_B
      IAH <- H_B
      G_obs <-  IAH %*% S_obs
      G_miss <-  IAH %*% S_miss
      K_y_miss <- Matrix::t(G_miss) %*% Diagonal(t_max*K, 1/sigma^2) %*% G_miss
      b_y_miss <- Matrix::t(G_miss) %*% Diagonal(t_max*K, 1/sigma^2) %*% ( c_const - G_obs %*% yobs_vec)
      
      if (aggregation == "identity"){
        # sample from N( K^-1 b, K^-1)
        U.chol <- Matrix::chol( K_y_miss )
        ymiss_vec <- backsolve(U.chol, backsolve(U.chol, b_y_miss,
                                                 upper.tri = T, transpose = T ) + rnorm(n_miss) ) 
      }
      if (aggregation == "average" || aggregation == "triangular"){
        K_bar_y_miss <- K_y_miss + K_prior_y_miss
        b_bar_y_miss <- b_y_miss + b_prior_y_miss
        U.chol <- Matrix::chol( K_bar_y_miss )
        ymiss_vec <- backsolve(U.chol, backsolve(U.chol, b_bar_y_miss,
                                                 upper.tri = T, transpose = T ) + rnorm(n_miss) ) 
      }
      
      yt_vec <- S_obs%*% yobs_vec + S_miss%*% ymiss_vec
    }


    
    # sample B and L
    {
      yt = matrix(yt_vec, nrow = K)
      y = t(yt)
      xt <- makeRegressor( y, y0, t_max, K, p)
      
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
          ZiSigi <- Zi %*% Diagonal(t_max, 1/sigma2[ii])
          dthetai <- iVthetai %*% thetai0 + ZiSigi %*% (yt[ii,] - Fmat[ii,])
        } else {
          Zi <- rbind(xt, Fmat)
          iVthetai <- Diagonal(x = 1/c(V_b_prior[ii,], rep(Vlhyper,r)))
          thetai0 <- c(b_prior[ii,], rep(mulhyper, r))
          ZiSigi <- Zi %*% Diagonal(t_max, 1/sigma2[ii])
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
    }
    
    {
      # Sample f
      e = matrix((yt - B %*%xt),ncol = 1)
      Xf <- kronecker(Diagonal(t_max), L)
      XfiSig <- Matrix::t(Xf) %*% Diagonal(t_max*K, 1/sigma2)
      Kf <- Diagonal(t_max*r, 1/sigmaf2) + XfiSig %*% Xf
      CKf <- chol(Kf)
      f <- backsolve( CKf,
                      backsolve( CKf, XfiSig %*% e,
                                 upper.tri = T, transpose = T )
                      + rnorm(t_max*r) )
      Fmat <- matrix(f, r, t_max)
    }
    
    # Sample sigma
    sigma2 <- rep(0,K)
    u_ort <- yt - B %*%xt - L %*% Fmat
    sigma_post_a <- sigma0_T0 + rep(t_max,K)
    sigma_post_b <- sigma0_S0 + rowSums(u_ort^2)
    for (i in c(1:K)){
      sigma2[i] <- rinvgamma(1, shape = sigma_post_a[i] * 0.5, rate = sigma_post_b[i] * 0.5)
    }
    sigma <- sqrt(sigma2)

    sigmaf2 <- rep(0,r)
    sigma_post_a <- sigma0_T0 + rep(t_max,r)
    sigma_post_b <- sigma0_S0 + rowSums(Fmat^2)
    for (i in c(1:r)){
      sigmaf2[i] <- rinvgamma(1, shape = sigma_post_a[i] * 0.5, rate = sigma_post_b[i] * 0.5)
    }
    sigmaf <- sqrt(sigmaf2)

    if ((j > inits$burnin) & (j %% inits$thin == 0)){
      mcmc[, (j - inits$burnin) %/% inits$thin] <- c(vec(B), L[L_idx], sigma, sigmaf)
      My_miss[, (j - inits$burnin) %/% inits$thin] <- ymiss_vec
    }
      
    if (j %% 1000 == 0) { cat(" Iteration ", j, " \n")}
  }
  nameL <- matrix(paste("L", repcol(c(1:K),r), reprow(c(1:r),K), sep = "_"), ncol = r)
  nameL <- nameL[L_idx]
  
  if (aggregation == "average" || aggregation == "triangular"){
    y_raw[,idq] <- NA
  }
  nameY <- matrix(paste("y", repcol(c(1:K),t_max), reprow(c(1:t_max),K), sep = "_"), nrow = K)
  row.names(My_miss) <- nameY[is.na(t(y_raw))]
  
  row.names(mcmc) <- c( paste("B0",c(1:K), sep = ""),
                        sprintf("B%d_%d_%d",reprow(c(1:p),K*K), rep(repcol(c(1:K),K), p), rep(reprow(c(1:K),K), p)),
                        nameL,
                        paste("sigma",c(1:K), sep = ""),
                        paste("sigmaf",c(1:r), sep = ""))
  return(list( param = as.mcmc(t(mcmc)),
               y_miss = as.mcmc(t(My_miss))
         ))
}

###########################################################################
#' @export
BMFVAR.Student.FHV <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL){
  # Init regressors in the right hand side
  t_max <- nrow(y)
  y_raw <- y
  yt_vec <- vec(t(y))
  yobs_vec <- na.exclude(yt_vec)
  
  # Init prior and initial values
  m = K * p + 1
  if (is.null(prior)){
    prior <- get_prior(y, p, priorStyle = "Minnesota", dist = "Student", 
                       SV = FALSE, aggregation = "identity", idq = K)
  }
  
  aggregation <- prior$aggregation
  idq <- prior$idq
  sdeviation <- prior$sdeviation
  ymiss_invprior <- 1/prior$ymiss_prior
  
  if (aggregation == "identity"){
    yt_vec <- vec(t(y))
    yobs_vec <- na.exclude(yt_vec)  
  }
  
  if (aggregation == "average" || aggregation == "triangular"){
    y[,idq] <- NA
    
    yt_vec <- vec(t(y))
    yobs_vec <- na.exclude(yt_vec)  
    z_obs_vec <- na.exclude(vec(t(y_raw[,idq])))  
    
    y_temp <- y
    y_temp[,idq] <- 0
    y_temp <- matrix(t(y_temp)[is.na(t(y))], ncol = 1)
    S_select <- Make_SectMat(y_temp)
    # Account for the missing in other columns rather than in low freq columns
    M_aggregate <- kronecker(Make_bandSparse(t_max, aggregation)[!is.na(y_raw[,idq[1] ]),],
                             Diagonal(length(idq))) %*% Matrix::t(S_select$S_obs)
  }
  
  # prior B
  b_prior = matrix(prior$b_prior, nrow = K)
  V_b_prior = matrix(diag(prior$V_b_prior), nrow = K)
  # prior sigma
  sigma0_T0 <- prior$sigma_T0
  sigma0_S0 <- prior$sigma_S0
  
  # prior nu
  nu_gam_a = prior$nu_gam_a
  nu_gam_b = prior$nu_gam_b
  
  # Initial values
  if (is.null(inits)){
    inits <- get_init(prior)
  }
  
  samples <- inits$samples
  B <- Vec_to_Mat(inits$B0, K, p)
  sigma <- inits$sigma
  sigma2 <- sigma^2
  
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
  
  sigmaf2 <- rep(1, r) # Variance of factor
  n_load <- K*r-r*(r+1)/2
  L_idx <- c(1:(K*r))[lower.tri(matrix(1, nrow = K, ncol = r))]
  L <- inits$A0[,1:r]
  Fmat <- matrix(0, nrow = r, ncol = t_max)
  
  S_select <- Make_SectMat(y)
  n_obs <- sum(!is.na(y))
  n_miss <- sum(is.na(y))
  n_cond <- length(y_raw[,idq]) - sum(is.na(y_raw[,idq]))
  S_obs <- S_select$S_obs
  S_miss <- S_select$S_miss
  K_prior_y_miss <- Matrix::t(M_aggregate) %*% Diagonal(n_cond, 1/sdeviation) %*% M_aggregate + Diagonal(n_miss, ymiss_invprior)
  b_prior_y_miss <- Matrix::t(M_aggregate) %*% Diagonal(n_cond, 1/sdeviation) %*% ( z_obs_vec)
  
  # Output
  mcmc <- matrix(NA, nrow = m*K + length(L_idx) + K+r + 1 + t_max,
                 ncol = (samples - inits$burnin)%/% inits$thin)
  My_miss <- matrix(NA, nrow = n_miss,
                    ncol = (samples - inits$burnin)%/% inits$thin)
  
  for (j in c(1:samples)){
    
    # # Sample y missing
    {
      H_B <- Make_HB(B, K, p ,t_max)
      c_B <- Make_cB(y0, B, K, p ,t_max, LF = as.numeric(L %*% Fmat) )
      # IA <- kronecker(Diagonal(t_max), Diagonal(K) ) # t_max * K diag
      c_const <- c_B
      IAH <- H_B
      G_obs <-  IAH %*% S_obs
      G_miss <-  IAH %*% S_miss
      Wh <- rep(1/w_sample, each = K) * rep(1/sigma^2, t_max)
      K_y_miss <- Matrix::t(G_miss) %*% Diagonal(t_max*K, Wh) %*% G_miss
      b_y_miss <- Matrix::t(G_miss) %*% Diagonal(t_max*K, Wh) %*% ( c_const - G_obs %*% yobs_vec)
      
      if (aggregation == "identity"){
        # sample from N( K^-1 b, K^-1)
        U.chol <- Matrix::chol( K_y_miss )
        ymiss_vec <- backsolve(U.chol, backsolve(U.chol, b_y_miss,
                                                 upper.tri = T, transpose = T ) + rnorm(n_miss) ) 
      }
      if (aggregation == "average" || aggregation == "triangular"){
        K_bar_y_miss <- K_y_miss + K_prior_y_miss
        b_bar_y_miss <- b_y_miss + b_prior_y_miss
        U.chol <- Matrix::chol( K_bar_y_miss )
        ymiss_vec <- backsolve(U.chol, backsolve(U.chol, b_bar_y_miss,
                                                 upper.tri = T, transpose = T ) + rnorm(n_miss) ) 
      }
      
      yt_vec <- S_obs%*% yobs_vec + S_miss%*% ymiss_vec
    }
    
    # sample B and L
    {
      yt = matrix(yt_vec, nrow = K)
      y = t(yt)
      xt <- makeRegressor( y, y0, t_max, K, p)
      
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
          ZiSigi <- Zi %*% Diagonal(t_max, 1/sigma2[ii]/ w[ii,] )
          dthetai <- iVthetai %*% thetai0 + ZiSigi %*% (yt[ii,] - Fmat[ii,])
        } else {
          Zi <- rbind(xt, Fmat)
          iVthetai <- Diagonal(x = 1/c(V_b_prior[ii,], rep(Vlhyper,r)))
          thetai0 <- c(b_prior[ii,], rep(mulhyper, r))
          ZiSigi <- Zi %*% Diagonal(t_max, 1/sigma2[ii]/ w[ii,])
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
    }
    
    {
      # Sample f
      e = matrix((yt - B %*%xt),ncol = 1)
      Xf <- kronecker(Diagonal(t_max), L)
      XfiSig <- Matrix::t(Xf) %*% Diagonal(t_max*K, 1/sigma2 /as.numeric(w) )
      Kf <- Diagonal(t_max*r, 1/sigmaf2) + XfiSig %*% Xf
      CKf <- chol(Kf)
      f <- backsolve( CKf,
                      backsolve( CKf, XfiSig %*% e,
                                 upper.tri = T, transpose = T )
                      + rnorm(t_max*r) )
      Fmat <- matrix(f, r, t_max)
    }
    # Sample sigma
    sigma2 <- rep(0,K)
    u_ort <- (yt - B %*%xt - L %*% Fmat)/w_sqrt
    sigma_post_a <- sigma0_T0 + rep(t_max,K)
    sigma_post_b <- sigma0_S0 + rowSums(u_ort^2)
    for (i in c(1:K)){
      sigma2[i] <- rinvgamma(1, shape = sigma_post_a[i] * 0.5, rate = sigma_post_b[i] * 0.5)
    }
    sigma <- sqrt(sigma2)
    
    sigmaf2 <- rep(0,r)
    sigma_post_a <- sigma0_T0 + rep(t_max,r)
    sigma_post_b <- sigma0_S0 + rowSums(Fmat^2)
    for (i in c(1:r)){
      sigmaf2[i] <- rinvgamma(1, shape = sigma_post_a[i] * 0.5, rate = sigma_post_b[i] * 0.5)
    }
    sigmaf <- sqrt(sigmaf2)
    
    # Sample w
    u <- (yt - B %*%xt - L %*% Fmat)
    shape_w <- nu*0.5 + K*0.5
    rate_w <- as.numeric(nu*0.5 + 0.5 * colSums((u / sigma)^2))
    w_sample <- rinvgamma( n=t_max, shape = shape_w, rate = rate_w )
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
    
    if ((j > inits$burnin) & (j %% inits$thin == 0)){
      mcmc[, (j - inits$burnin) %/% inits$thin] <- c(vec(B), L[L_idx], sigma, sigmaf, nu, as.numeric(w_sample))
      My_miss[, (j - inits$burnin) %/% inits$thin] <- ymiss_vec
    }
    
    if (j %% 1000 == 0) {
      cat(" Iteration ", j, " ", logsigma_nu," ", round(nu,2)," \n")
      acount_w <- rep(0,t_max)
    }
  }
  nameL <- matrix(paste("L", repcol(c(1:K),r), reprow(c(1:r),K), sep = "_"), ncol = r)
  nameL <- nameL[L_idx]
  
  if (aggregation == "average" || aggregation == "triangular"){
    y_raw[,idq] <- NA
  }
  nameY <- matrix(paste("y", repcol(c(1:K),t_max), reprow(c(1:t_max),K), sep = "_"), nrow = K)
  row.names(My_miss) <- nameY[is.na(t(y_raw))]
  
  row.names(mcmc) <- c( paste("B0",c(1:K), sep = ""),
                        sprintf("B%d_%d_%d",reprow(c(1:p),K*K), rep(repcol(c(1:K),K), p), rep(reprow(c(1:K),K), p)),
                        nameL,
                        paste("sigma",c(1:K), sep = ""),
                        paste("sigmaf",c(1:r), sep = ""),
                        paste("nu"),
                        paste("w",c(1:t_max), sep = ""))
  
  return(list( param = as.mcmc(t(mcmc)),
               y_miss = as.mcmc(t(My_miss))
  ))
}
###########################################################################
#' @export
BMFVAR.OT.FHV <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL){
  # Init regressors in the right hand side
  t_max <- nrow(y)
  y_raw <- y
  
  # Init prior and initial values
  m = K * p + 1
  if (is.null(prior)){
    prior <- get_prior(y, p, priorStyle = "Minnesota", dist = "OT", 
                       SV = FALSE, aggregation = "identity", idq = K)
  }
  
  aggregation <- prior$aggregation
  idq <- prior$idq
  sdeviation <- prior$sdeviation
  ymiss_invprior <- 1/prior$ymiss_prior
  
  if (aggregation == "identity"){
    yt_vec <- vec(t(y))
    yobs_vec <- na.exclude(yt_vec)  
  }
  
  if (aggregation == "average" || aggregation == "triangular"){
    y[,idq] <- NA
    
    yt_vec <- vec(t(y))
    yobs_vec <- na.exclude(yt_vec)  
    z_obs_vec <- na.exclude(vec(t(y_raw[,idq])))  
    
    y_temp <- y
    y_temp[,idq] <- 0
    y_temp <- matrix(t(y_temp)[is.na(t(y))], ncol = 1)
    S_select <- Make_SectMat(y_temp)
    # Account for the missing in other columns rather than in low freq columns
    M_aggregate <- kronecker(Make_bandSparse(t_max, aggregation)[!is.na(y_raw[,idq[1] ]),],
                             Diagonal(length(idq))) %*% Matrix::t(S_select$S_obs)
  }
  
  # prior B
  b_prior = matrix(prior$b_prior, nrow = K)
  V_b_prior = matrix(diag(prior$V_b_prior), nrow = K)
  # prior sigma
  sigma0_T0 <- prior$sigma_T0
  sigma0_S0 <- prior$sigma_S0
  # prior A
  a_prior = prior$a_prior
  V_a_prior = prior$V_a_prior
  # prior nu
  nu_gam_a = prior$nu_gam_a
  nu_gam_b = prior$nu_gam_b
  
  # Initial values
  if (is.null(inits)){
    inits <- get_init(prior)
  }
  samples <- inits$samples
  B <- Vec_to_Mat(inits$B0, K, p)
  sigma <- inits$sigma
  sigma2 <- sigma^2
  
  # Multi degrees of freedom
  nu <- inits$nu
  logsigma_nu <- rep(0,K)
  acount_nu <- rep(0,K)
  acount_w <- rep(0, t_max)
  # Init w as Gaussian
  w_sample <- rep(1, t_max)
  w <- reprow(w_sample, K)
  w_sqrt <- sqrt(w)
  w_sqrt_inv <- 1/w_sqrt
  
  # new precompute
  r <- prior$r
  Vlhyper <- prior$Vlhyper
  mulhyper <- prior$mulhyper
  
  sigmaf2 <- rep(1, r) # Variance of factor
  n_load <- K*r-r*(r+1)/2
  L_idx <- c(1:(K*r))[lower.tri(matrix(1, nrow = K, ncol = r))]
  L <- inits$A0[,1:r]
  Fmat <- matrix(0, nrow = r, ncol = t_max)
  
  S_select <- Make_SectMat(y)
  n_obs <- sum(!is.na(y))
  n_miss <- sum(is.na(y))
  n_cond <- length(y_raw[,idq]) - sum(is.na(y_raw[,idq]))
  S_obs <- S_select$S_obs
  S_miss <- S_select$S_miss
  K_prior_y_miss <- Matrix::t(M_aggregate) %*% Diagonal(n_cond, 1/sdeviation) %*% M_aggregate + Diagonal(n_miss, ymiss_invprior)
  b_prior_y_miss <- Matrix::t(M_aggregate) %*% Diagonal(n_cond, 1/sdeviation) %*% ( z_obs_vec)
  
  # Output
  mcmc <- matrix(NA, nrow = m*K + length(L_idx) + K+r + K + K*t_max,
                 ncol = (samples - inits$burnin)%/% inits$thin)
  My_miss <- matrix(NA, nrow = n_miss,
                    ncol = (samples - inits$burnin)%/% inits$thin)
  
  for (j in c(1:samples)){
    
    # # Sample y missing
    H_B <- Make_HB(B, K, p ,t_max)
    c_B <- Make_cB(y0, B, K, p ,t_max, LF = as.numeric(L %*% Fmat) )
    # IA <- kronecker(Diagonal(t_max), Diagonal(K) ) # t_max * K diag
    c_const <- c_B
    IAH <- H_B
    G_obs <-  IAH %*% S_obs
    G_miss <-  IAH %*% S_miss
    Wh <- 1/vec(w) * rep(1/sigma^2, t_max)
    K_y_miss <- Matrix::t(G_miss) %*% Diagonal(t_max*K, Wh) %*% G_miss
    b_y_miss <- Matrix::t(G_miss) %*% Diagonal(t_max*K, Wh) %*% ( c_const - G_obs %*% yobs_vec)
    
    if (aggregation == "identity"){
      # sample from N( K^-1 b, K^-1)
      U.chol <- Matrix::chol( K_y_miss )
      ymiss_vec <- backsolve(U.chol, backsolve(U.chol, b_y_miss,
                                               upper.tri = T, transpose = T ) + rnorm(n_miss) ) 
    }
    if (aggregation == "average" || aggregation == "triangular"){
      K_bar_y_miss <- K_y_miss + K_prior_y_miss
      b_bar_y_miss <- b_y_miss + b_prior_y_miss
      U.chol <- Matrix::chol( K_bar_y_miss )
      ymiss_vec <- backsolve(U.chol, backsolve(U.chol, b_bar_y_miss,
                                               upper.tri = T, transpose = T ) + rnorm(n_miss) ) 
    }
    
    yt_vec <- S_obs%*% yobs_vec + S_miss%*% ymiss_vec
    
    # sample B and L
    {
      yt = matrix(yt_vec, nrow = K)
      y = t(yt)
      xt <- makeRegressor( y, y0, t_max, K, p)
      
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
          ZiSigi <- Zi %*% Diagonal(t_max, 1/sigma2[ii]/ w[ii,] )
          dthetai <- iVthetai %*% thetai0 + ZiSigi %*% (yt[ii,] - Fmat[ii,])
        } else {
          Zi <- rbind(xt, Fmat)
          iVthetai <- Diagonal(x = 1/c(V_b_prior[ii,], rep(Vlhyper,r)))
          thetai0 <- c(b_prior[ii,], rep(mulhyper, r))
          ZiSigi <- Zi %*% Diagonal(t_max, 1/sigma2[ii]/ w[ii,])
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
    }
    
    {
      # Sample f
      e = matrix((yt - B %*%xt),ncol = 1)
      Xf <- kronecker(Diagonal(t_max), L)
      XfiSig <- Matrix::t(Xf) %*% Diagonal(t_max*K, 1/sigma2 /as.numeric(w) )
      Kf <- Diagonal(t_max*r, 1/sigmaf2) + XfiSig %*% Xf
      CKf <- chol(Kf)
      f <- backsolve( CKf,
                      backsolve( CKf, XfiSig %*% e,
                                 upper.tri = T, transpose = T )
                      + rnorm(t_max*r) )
      Fmat <- matrix(f, r, t_max)
    }
    
    # Sample sigma
    u_ort <- (yt - B %*%xt - L %*% Fmat)/w_sqrt
    sigma_post_a <- sigma0_T0 + rep(t_max,K)
    sigma_post_b <- sigma0_S0 + apply(u_ort^2, 1, sum)
    for (i in c(1:K)){
      sigma2[i] <- rinvgamma(1, shape = sigma_post_a[i] * 0.5, rate = sigma_post_b[i] * 0.5)
    }
    sigma <- sqrt(sigma2)
    
    sigmaf2 <- rep(0,r)
    sigma_post_a <- sigma0_T0 + rep(t_max,r)
    sigma_post_b <- sigma0_S0 + rowSums(Fmat^2)
    for (i in c(1:r)){
      sigmaf2[i] <- rinvgamma(1, shape = sigma_post_a[i] * 0.5, rate = sigma_post_b[i] * 0.5)
    }
    sigmaf <- sqrt(sigmaf2)
    
    # Sample w
    u <- (yt - B %*%xt - L %*% Fmat)
    w <- matrix( rinvgamma( n = K*t_max, shape = (nu+1)*0.5, rate = 0.5*( nu + (u^2)/sigma2 ) ), K, t_max )
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
    if ((j > inits$burnin) & (j %% inits$thin == 0)){
      mcmc[, (j - inits$burnin) %/% inits$thin] <- c(vec(B), L[L_idx], sigma, sigmaf, nu, as.numeric(w))
      My_miss[, (j - inits$burnin) %/% inits$thin] <- ymiss_vec
    }
    
    if (j %% 1000 == 0) {
      cat(" Iteration ", j, " ", logsigma_nu," ", min(acount_w)," ", max(acount_w)," ", mean(acount_w), " ", round(nu,2), " \n")
      acount_w <- rep(0,t_max)
    }
  }
  nameL <- matrix(paste("L", repcol(c(1:K),r), reprow(c(1:r),K), sep = "_"), ncol = r)
  nameL <- nameL[L_idx]
  
  if (aggregation == "average" || aggregation == "triangular"){
    y_raw[,idq] <- NA
  }
  nameY <- matrix(paste("y", repcol(c(1:K),t_max), reprow(c(1:t_max),K), sep = "_"), nrow = K)
  row.names(My_miss) <- nameY[is.na(t(y_raw))]
  
  row.names(mcmc) <- c( paste("B0",c(1:K), sep = ""),
                        sprintf("B%d_%d_%d",reprow(c(1:p),K*K), rep(repcol(c(1:K),K), p), rep(reprow(c(1:K),K), p)),
                        nameL,
                        paste("sigma",c(1:K), sep = ""),
                        paste("sigmaf",c(1:r), sep = ""),
                        paste("nu",c(1:K), sep = ""),
                        sprintf("w_%d_%d", repcol(c(1:K),t_max), reprow(c(1:t_max),K))
  )
  return(list( param = as.mcmc(t(mcmc)),
               y_miss = as.mcmc(t(My_miss))
  ))
}


