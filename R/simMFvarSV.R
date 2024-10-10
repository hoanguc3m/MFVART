#' Simulate data from MF VAR model with SV
#'
#' This function simulate data from a MF VAR model with SV.
#' \deqn{y_t = B x_t + SQRT(w_t) A^(-1) Sigma_t eps_t}
#' @param dist The variable specifies the BVAR error distribution. It should be one of
#' c("Gaussian","Student","Skew.Student","Skew.Student", "MT","Skew.MT","MST").
#' @param K The number of variables in BVAR model.
#' @param p The number of lags in BVAR model.
#' @param t_max The number of observations
#' @param b0 The coefficients of B matrix.
#' @param a0 The coefficients of A matrix.
#' @param h The initial log diag(Sigma_t) matrix. sigma_h is set at diag(0.03, K)
#' @param nu The degree of freedom.
#' @param gamma The skewness from \eqn{[-gamma, gamma]}.
#' @param idq The indicator positions of the quarter variables
#' @param freq The frequency of observed quarter variables.
#' @param seednum The default seed.
#' @param burn_in The discarded observations.
#' @return A list of simulated data with its specification.
#' @export
#' @examples
#' \dontrun{
#' datagen <- sim.MFVAR.SV(dist="Gaussian")
#' y <- datagen$y
#' }
sim.MFVAR.SV <- function(dist, K = 5, p = 2, t_max = 1000,
                          b0 = 0.5, a0 = 0.5, h = 0, sigma_h = NULL,
                          y0 = matrix(0, ncol = K, nrow = p),
                          nu = 6, gamma = 0.5, sigma_G = NULL,
                          idq = rep(0,K), freq = 3, seednum = 0, burn_in = 0){
  if (!(dist %in% c("Gaussian","Student",
                    "MT","OT") ))
    stop("dist is not implemented.")

  if (dist == "Gaussian") datagen <- sim.MFVAR.Gaussian.SV(K, p, t_max, b0, a0, h, sigma_h, y0, idq, freq, seednum, burn_in)
  if (dist == "Student") datagen <- sim.MFVAR.Student.SV(K, p, t_max, b0, a0, h, sigma_h, y0, nu, idq, freq, seednum, burn_in)
  if (dist == "MT") datagen <- sim.MFVAR.MT.SV(K, p, t_max, b0, a0, h, sigma_h, y0, nu, idq, freq, seednum, burn_in)
  if (dist == "OT") datagen <- sim.MFVAR.OT.SV(K, p, t_max, b0, a0, h, sigma_h, y0, nu, idq, freq, seednum, burn_in)

  return(datagen)
}

#' @export
sim.MFVAR.Gaussian.SV <- function(K = 5, p = 2, t_max = 1000,
                                b0 = 0.5, a0 = 0.5, h = 0, sigma_h = NULL,
                                y0 = matrix(0, ncol = K, nrow = p),
                                idq = rep(0,K), freq = 3, seednum = 0, burn_in = 0){
  t_max = t_max + burn_in
  set.seed(seednum)
  # Sample matrix coefficient B
  B0 <- cbind(rep(0,K))
  if (p > 0){
    if (length(b0) == 1) {
      for (i in c(1:p)){
        B0 <- cbind(B0, b0^i*diag(K))
      }
    } else {
      B0 <- matrix(b0, nrow = K)
    }
  } else {
    B0 <- matrix(b0, nrow = K)
  }


  # Sample matrix corr A0
  if (length(a0) == 1) {
    A0 <- matrix(a0, K, K)
    diag(A0) <- 1
    A0[upper.tri(A0)] <- 0
  } else {
    A0 <- matrix(0, nrow = K, ncol = K)
    A0[upper.tri(A0)] <- a0
    A0 <- t(A0)
    diag(A0) <- 1
  }
  # Sample matrix variance h
  if (length(h) == 1){
    h <- rep(h,K)
  }
  # No skew
  # No tail
  # Volatility volatility
  if (is.null(sigma_h)){
    Vh <- seq(3e-2, 3e-2, length.out = K)
  } else {
    Vh <-  sigma_h
  }

  ystar <- tail(y0, p)
  y_mean <- matrix(NA, nrow = t_max, ncol = K)
  y_var <- matrix(NA, nrow = t_max, ncol = 0.5*K*(K+1))
  volatility <- matrix(NA, nrow = t_max, ncol = K)

  eps <- matrix(rnorm(t_max*K), ncol = K)
  volatility <- matrix(NA, nrow = t_max, ncol = K)
  inv_A0 <- solve(A0)

  for (i in c(1:t_max)){
    h <- h +  Vh * rnorm(K)
    Sigma_t <-  inv_A0 %*% diag(exp(0.5*h), nrow = K)
    Sigma2_t <- Sigma_t %*% t(Sigma_t)
    y_var[i,] <- Sigma_t[lower.tri(Sigma_t, diag = T)]
    volatility[i,] <- diag(Sigma2_t)
    if (p > 0){
      xt <- rbind(1, vec( t(ystar[(p+i-1):i,])))
    } else {
      xt <- matrix(1,nrow = K)
    }
    y_mean[i,] <- B0 %*% xt
    ysim <-  B0 %*% xt + Sigma_t %*% eps[i,]
    ystar <- rbind(ystar, t(ysim))
  }

  t_max = t_max - burn_in
  y = y_true = as.matrix(ystar[(p+burn_in+1):(p+burn_in+t_max),], nrow = t_max)
  y[-seq(1, t_max, by = freq) ,idq] = NA
  
  list(y = y, y_true = y_true,
       y0 = y0, y_mean = y_mean, y_var = y_var, volatility = volatility,
       K = K, p = p, t_max = t_max,
       A0 = A0, B0 = matrix(B0, nrow = K),
       Vh = Vh,
       dist = "Gaussian", SV = TRUE)
}

#' @export
sim.MFVAR.Student.SV <- function(K = 5, p = 2, t_max = 1000,
                               b0 = 0.5, a0 = 0.5, h = 0, sigma_h = NULL,
                               y0 = matrix(0, ncol = K, nrow = p),
                               nu = 6, idq = rep(0,K), freq = 3, seednum = 0, burn_in = 0){
  t_max = t_max + burn_in
  set.seed(seednum)
  # Sample matrix coefficient B
  B0 <- cbind(rep(0,K))
  if (p > 0){
    if (length(b0) == 1) {
      for (i in c(1:p)){
        B0 <- cbind(B0, b0^i*diag(K))
      }
    } else {
      B0 <- matrix(b0, nrow = K)
    }
  } else {
    B0 <- matrix(b0, nrow = K)
  }


  # Sample matrix corr A0
  if (length(a0) == 1) {
    A0 <- matrix(a0, K, K)
    diag(A0) <- 1
    A0[upper.tri(A0)] <- 0
  } else {
    A0 <- matrix(0, nrow = K, ncol = K)
    A0[upper.tri(A0)] <- a0
    A0 <- t(A0)
    diag(A0) <- 1
  }
  # Sample matrix variance h
  if (length(h) == 1){
    h <- rep(h,K)
  }
  # No skew
  # Tail of student
  w_t <- rinvgamma(t_max, shape = nu/2, rate = nu/2)
  w_sqrt_t <- sqrt(w_t)

  # Volatility volatility
  if (is.null(sigma_h)){
    Vh <- seq(3e-2, 3e-2, length.out = K)
  } else {
    Vh <-  sigma_h
  }

  ystar <- tail(y0, p)
  y_mean <- matrix(NA, nrow = t_max, ncol = K)
  y_var <- matrix(NA, nrow = t_max, ncol = 0.5*K*(K+1))
  volatility <- matrix(NA, nrow = t_max, ncol = K)

  eps <- matrix(rnorm(t_max*K), ncol = K)
  volatility <- matrix(NA, nrow = t_max, ncol = K)
  inv_A0 <- solve(A0)

  for (i in c(1:t_max)){
    h <- h +  Vh * rnorm(K)
    Sigma_t <-  diag(w_sqrt_t[i], nrow = K) %*% inv_A0 %*% diag(exp(0.5*h), nrow = K)
    Sigma2_t <- Sigma_t %*% t(Sigma_t)
    y_var[i,] <- Sigma_t[lower.tri(Sigma_t, diag = T)]
    volatility[i,] <- diag(Sigma2_t)
    if (p > 0){
      xt <- rbind(1, vec( t(ystar[(p+i-1):i,])))
    } else {
      xt <- matrix(1,nrow = K)
    }
    y_mean[i,] <- B0 %*% xt
    ysim <-  B0 %*% xt + Sigma_t %*% eps[i,]
    ystar <- rbind(ystar, t(ysim))
  }

  t_max = t_max - burn_in
  y = y_true = as.matrix(ystar[(p+burn_in+1):(p+burn_in+t_max),], nrow = t_max)
  y[-seq(1, t_max, by = freq) ,idq] = NA
  
  list(y = y, y_true = y_true,
       y0 = y0, y_mean = y_mean, y_var = y_var, volatility = volatility,
       K = K, p = p, t_max = t_max,
       A0 = A0, B0 = matrix(B0, nrow = K),
       nu = nu, w = w_t[(burn_in+1):(burn_in+t_max)],
       Vh = Vh,
       dist = "Student", SV = TRUE)
}


#' @export
sim.MFVAR.MT.SV <- function(K = 5, p = 2, t_max = 1000,
                                    b0 = 0.5, a0 = 0.5, h = 0, sigma_h = NULL,
                                    y0 = matrix(0, ncol = K, nrow = p),
                                    nu = 6, idq = rep(0,K), freq = 3, seednum = 0, burn_in = 0){
  t_max = t_max + burn_in
  set.seed(seednum)
  # Sample matrix coefficient B
  B0 <- cbind(rep(0,K))
  if (p > 0){
    if (length(b0) == 1) {
      for (i in c(1:p)){
        B0 <- cbind(B0, b0^i*diag(K))
      }
    } else {
      B0 <- matrix(b0, nrow = K)
    }
  } else {
    B0 <- matrix(b0, nrow = K)
  }


  # Sample matrix corr A0
  if (length(a0) == 1) {
    A0 <- matrix(a0, K, K)
    diag(A0) <- 1
    A0[upper.tri(A0)] <- 0
  } else {
    A0 <- matrix(0, nrow = K, ncol = K)
    A0[upper.tri(A0)] <- a0
    A0 <- t(A0)
    diag(A0) <- 1
  }
  # Sample matrix variance h
  if (length(h) == 1){
    h <- rep(h,K)
  }

  # No skew
  # Tail of student
  if (length(nu) == 1) nu = rep(nu,K)
  w_t <- mapply(rinvgamma, n = t_max, shape = nu/2, rate = nu/2)
  w_sqrt_t <- sqrt(w_t)

  # Volatility volatility
  if (is.null(sigma_h)){
    Vh <- seq(3e-2, 3e-2, length.out = K)
  } else {
    Vh <-  sigma_h
  }

  ystar <- tail(y0, p)
  y_mean <- matrix(NA, nrow = t_max, ncol = K)
  y_var <- matrix(NA, nrow = t_max, ncol = 0.5*K*(K+1))
  volatility <- matrix(NA, nrow = t_max, ncol = K)

  eps <- matrix(rnorm(t_max*K), ncol = K)
  volatility <- matrix(NA, nrow = t_max, ncol = K)
  inv_A0 <- solve(A0)

  for (i in c(1:t_max)){
    h <- h +  Vh * rnorm(K)
    Sigma_t <-  diag(w_sqrt_t[i,]) %*% inv_A0 %*% diag(exp(0.5*h), nrow = K)
    Sigma2_t <- Sigma_t %*% t(Sigma_t)
    y_var[i,] <- Sigma_t[lower.tri(Sigma_t, diag = T)]
    volatility[i,] <- diag(Sigma2_t)
    if (p > 0){
      xt <- rbind(1, vec( t(ystar[(p+i-1):i,])))
    } else {
      xt <- matrix(1,nrow = K)
    }
    y_mean[i,] <- B0 %*% xt
    ysim <-  B0 %*% xt + Sigma_t %*% eps[i,]
    ystar <- rbind(ystar, t(ysim))
  }

  t_max = t_max - burn_in
  y = y_true = as.matrix(ystar[(p+burn_in+1):(p+burn_in+t_max),], nrow = t_max)
  y[-seq(1, t_max, by = freq) ,idq] = NA
  
  list(y = y, y_true = y_true,
       y0 = y0, y_mean = y_mean, y_var = y_var, volatility = volatility,
       K = K, p = p, t_max = t_max,
       A0 = A0, B0 = matrix(B0, nrow = K),
       nu = nu, w = w_t[(burn_in+1):(burn_in+t_max),],
       Vh = Vh,
       dist = "MT", SV = TRUE)
}


#' @export
sim.MFVAR.OT.SV <- function(K = 5, p = 2, t_max = 1000,
                                        b0 = 0.5, a0 = 0.5, h = 0, sigma_h = NULL,
                                        y0 = matrix(0, ncol = K, nrow = p),
                                        nu = 6, idq = rep(0,K), freq = 3, seednum = 0, burn_in = 0){
  t_max = t_max + burn_in
  set.seed(seednum)
  # Sample matrix coefficient B
  B0 <- cbind(rep(0,K))
  if (p > 0){
    if (length(b0) == 1) {
      for (i in c(1:p)){
        B0 <- cbind(B0, b0^i*diag(K))
      }
    } else {
      B0 <- matrix(b0, nrow = K)
    }
  } else {
    B0 <- matrix(b0, nrow = K)
  }


  # Sample matrix corr A0
  if (length(a0) == 1) {
    A0 <- matrix(a0, K, K)
    diag(A0) <- 1
    A0[upper.tri(A0)] <- 0
  } else {
    A0 <- matrix(0, nrow = K, ncol = K)
    A0[upper.tri(A0)] <- a0
    A0 <- t(A0)
    diag(A0) <- 1
  }
  # Sample matrix variance h
  if (length(h) == 1){
    h <- rep(h,K)
  }

  # No skew
  # Tail of student
  if (length(nu) == 1) nu = rep(nu,K)
  w_t <- mapply(rinvgamma, n = t_max, shape = nu/2, rate = nu/2)
  w_sqrt_t <- sqrt(w_t)

  # Volatility volatility
  if (is.null(sigma_h)){
    Vh <- seq(3e-2, 3e-2, length.out = K)
  } else {
    Vh <-  sigma_h
  }

  ystar <- tail(y0, p)
  y_mean <- matrix(NA, nrow = t_max, ncol = K)
  y_var <- matrix(NA, nrow = t_max, ncol = 0.5*K*(K+1))
  volatility <- matrix(NA, nrow = t_max, ncol = K)

  eps <- matrix(rnorm(t_max*K), ncol = K)
  volatility <- matrix(NA, nrow = t_max, ncol = K)
  inv_A0 <- solve(A0)

  for (i in c(1:t_max)){
    h <- h +  Vh * rnorm(K)
    Sigma_t <-  inv_A0 %*% diag(w_sqrt_t[i,] * exp(0.5*h))
    Sigma2_t <- Sigma_t %*% t(Sigma_t)
    y_var[i,] <- Sigma_t[lower.tri(Sigma_t, diag = T)]
    volatility[i,] <- diag(Sigma2_t)
    if (p > 0){
      xt <- rbind(1, vec( t(ystar[(p+i-1):i,])))
    } else {
      xt <- matrix(1,nrow = K)
    }
    y_mean[i,] <- B0 %*% xt
    ysim <-  B0 %*% xt + Sigma_t %*% eps[i,]
    ystar <- rbind(ystar, t(ysim))
  }

  t_max = t_max - burn_in
  y = y_true = as.matrix(ystar[(p+burn_in+1):(p+burn_in+t_max),], nrow = t_max)
  y[-seq(1, t_max, by = freq) ,idq] = NA
  
  list(y = y, y_true = y_true,
       y0 = y0, y_mean = y_mean, y_var = y_var, volatility = volatility,
       K = K, p = p, t_max = t_max,
       A0 = A0, B0 = matrix(B0, nrow = K),
       nu = nu, w = w_t[(burn_in+1):(burn_in+t_max),],
       Vh = Vh,
       dist = "OT", SV = TRUE)
}

