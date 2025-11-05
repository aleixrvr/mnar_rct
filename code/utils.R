library(MASS)

range_tolerance <- 0.05
logistic <- plogis

bound_to_margin <- function(prob, prob_margin){
  ifelse(prob < prob_margin, prob_margin, 
         ifelse(prob > (1 - prob_margin), 1 - prob_margin, 
                prob))
}


generate_data <- function(
    base_s, 
    impact_t_s, 
    impact_x_s, 
    impact_c_s, 
    base_o, 
    impact_t_o, 
    impact_x_o, 
    impact_c_o, 
    add_impact_s_o, 
    base_m, 
    impact_o_m, 
    impact_x_m, 
    impact_s_m, 
    p_c, 
    p_x, 
    p_t,
    n,
    prob_margin= 0.1
) {
  p_s_ntncnx <- base_s
  p_s_ntncx <- base_s + impact_x_s
  p_s_ntcnx <- base_s + impact_c_s
  p_s_ntcx <- base_s + impact_x_s + impact_c_s
  p_s_tncnx <- base_s + impact_t_s
  p_s_tncx <- base_s + impact_t_s + impact_x_s
  p_s_tcnx <- base_s + impact_t_s + impact_c_s
  p_s_tcx <- base_s + impact_t_s + impact_x_s + impact_c_s
  
  p_o_ntncnx <- base_o
  p_o_ntncx <- base_o + impact_x_o
  p_o_ntcnx <- base_o + impact_c_o
  p_o_ntcx <- base_o + impact_x_o + impact_c_o
  p_o_tncnx <- base_o + impact_t_o
  p_o_tncx <- base_o + impact_t_o + impact_x_o
  p_o_tcnx <- base_o + impact_t_o + impact_c_o
  p_o_tcx <- base_o + impact_t_o + impact_x_o + impact_c_o
  
  # p_a_nsnonx <- 1 - (base_m + impact_no_m)
  # p_a_nsnox <- 1 - (base_m + impact_no_m + impact_x_m)
  # p_a_nsonx <- 1 - base_m 
  # p_a_nsox <- 1 - (base_m + impact_x_m)
  # p_a_snonx <- 1 - (base_m + impact_s_m + impact_no_m)
  # p_a_snox <- 1 - (base_m + impact_s_m + impact_no_m + impact_x_m)
  # p_a_sonx <- 1 - (base_m + impact_s_m)
  # p_a_sox <- 1 - (base_m + impact_s_m + impact_x_m)
  
  t <- rbinom(n, 1, p_t)
  x <- rbinom(n, 1, p_x)
  c <- rbinom(n, 1, p_c)
  
  s_expect <- p_s_ntncnx*(1-t)*(1-c)*(1-x) + p_s_ntncx*(1-t)*(1-c)*x + 
    p_s_ntcnx*(1-t)*c*(1-x) + p_s_ntcx*(1-t)*c*x  
  s_expect <- s_expect + p_s_tncnx*t*(1-c)*(1-x) + p_s_tncx*t*(1-c)*x + 
    p_s_tcnx*t*c*(1-x) + p_s_tcx*t*c*x  
  s <- rbinom(n, 1, bound_to_margin(s_expect, prob_margin))
  
  o_expect <- base_o + impact_t_o*t + impact_x_o*x + impact_c_o*c + 
    (impact_t_o + impact_x_o + impact_c_o)*t*x*c
  o_expect <- o_expect + add_impact_s_o*s 
  o <- rbinom(n, 1, bound_to_margin(o_expect, prob_margin))
  
  m_expect <- base_m + impact_o_m*o + impact_x_m*x + impact_s_m*s + 
    (impact_o_m + impact_x_m + impact_s_m)*o*x*s
  
  a <- rbinom(n, 1, bound_to_margin(1 - m_expect, prob_margin))
  
  
  # p_o_t <- p_o_tcx*p_c*p_x*p_t + p_o_tcnx*p_c*(1-p_x)*p_t + p_o_tncx*(1-p_c)*p_x*p_t + p_o_tncnx*(1-p_c)*(1-p_x)*p_t
  # p_o_nt <- p_o_ntcx*p_c*p_x*(1-p_t) + p_o_ntcnx*p_c*(1-p_x)*(1-p_t) +  p_o_ntncx*(1-p_c)*p_x*(1-p_t) + p_o_ntncnx*(1-p_c)*(1-p_x)*(1-p_t)
  
  
  p_o_t <- (p_o_tcx + add_impact_s_o)*p_s_tcx*p_c*p_x + 
    (p_o_tcnx + add_impact_s_o)*p_s_tcnx*p_c*(1-p_x) + 
    (p_o_tncx + add_impact_s_o)*p_s_tncx*(1-p_c)*p_x + 
    (p_o_tncnx + add_impact_s_o)*p_s_tncnx*(1-p_c)*(1-p_x)
  p_o_t <- p_o_t + p_o_tcx*(1-p_s_tcx)*p_c*p_x + 
    p_o_tcnx*(1-p_s_tcnx)*p_c*(1-p_x) + 
    p_o_tncx*(1-p_s_tncx)*(1-p_c)*p_x + 
    p_o_tncnx*(1-p_s_tncnx)*(1-p_c)*(1-p_x)
  
  p_o_nt <- (p_o_ntcx + add_impact_s_o)*p_s_ntcx*p_c*p_x + 
    (p_o_ntcnx + add_impact_s_o)*p_s_ntcnx*p_c*(1-p_x) + 
    (p_o_ntncx + add_impact_s_o)*p_s_ntncx*(1-p_c)*p_x + 
    (p_o_ntncnx + add_impact_s_o)*p_s_ntncnx*(1-p_c)*(1-p_x)
  p_o_nt <- p_o_nt + p_o_ntcx*(1-p_s_ntcx)*p_c*p_x + 
    p_o_ntcnx*(1-p_s_ntcnx)*p_c*(1-p_x) + 
    p_o_ntncx*(1-p_s_ntncx)*(1-p_c)*p_x + 
    p_o_ntncnx*(1-p_s_ntncnx)*(1-p_c)*(1-p_x)
  
  y <- ifelse(a==1, o, NA)
  dt <- data.table(t, x, c, s, o, a, y)
  
  return(list(
    dt = dt,
    p_o_t = p_o_t,
    p_o_nt = p_o_nt
  ))
}

generate_data_ <- function(
    base_s, 
    impact_t_s, 
    impact_x_s, 
    impact_c_s, 
    base_o, 
    impact_t_o, 
    impact_x_o, 
    impact_c_o, 
    add_impact_s_o, 
    base_m, 
    impact_no_m, 
    impact_x_m, 
    impact_s_m, 
    p_c, 
    p_x, 
    p_t,
    n,
    prob_margin= 0.1
) {
  p_s_ntncnx <- base_s
  p_s_ntncx <- base_s + impact_x_s
  p_s_ntcnx <- base_s + impact_c_s
  p_s_ntcx <- base_s + impact_x_s + impact_c_s
  p_s_tncnx <- base_s + impact_t_s
  p_s_tncx <- base_s + impact_t_s + impact_x_s
  p_s_tcnx <- base_s + impact_t_s + impact_c_s
  p_s_tcx <- base_s + impact_t_s + impact_x_s + impact_c_s
  
  p_o_ntncnx <- base_o
  p_o_ntncx <- base_o + impact_x_o
  p_o_ntcnx <- base_o + impact_c_o
  p_o_ntcx <- base_o + impact_x_o + impact_c_o
  p_o_tncnx <- base_o + impact_t_o
  p_o_tncx <- base_o + impact_t_o + impact_x_o
  p_o_tcnx <- base_o + impact_t_o + impact_c_o
  p_o_tcx <- base_o + impact_t_o + impact_x_o + impact_c_o
  
  p_a_nsnonx <- 1 - (base_m + impact_no_m)
  p_a_nsnox <- 1 - (base_m + impact_no_m + impact_x_m)
  p_a_nsonx <- 1 - base_m 
  p_a_nsox <- 1 - (base_m + impact_x_m)
  p_a_snonx <- 1 - (base_m + impact_s_m + impact_no_m)
  p_a_snox <- 1 - (base_m + impact_s_m + impact_no_m + impact_x_m)
  p_a_sonx <- 1 - (base_m + impact_s_m)
  p_a_sox <- 1 - (base_m + impact_s_m + impact_x_m)
  
  t <- rbinom(n, 1, p_t)
  x <- rbinom(n, 1, p_x)
  c <- rbinom(n, 1, p_c)
  
  s_expect <- p_s_ntncnx*(1-t)*(1-c)*(1-x) + p_s_ntncx*(1-t)*(1-c)*x + 
    p_s_ntcnx*(1-t)*c*(1-x) + p_s_ntcx*(1-t)*c*x  
  s_expect <- s_expect + p_s_tncnx*t*(1-c)*(1-x) + p_s_tncx*t*(1-c)*x + 
    p_s_tcnx*t*c*(1-x) + p_s_tcx*t*c*x  
  s <- rbinom(n, 1, bound_to_margin(s_expect, prob_margin))
  
  o_expect <- p_o_ntncnx*(1-t)*(1-c)*(1-x) + p_o_ntncx*(1-t)*(1-c)*x +
    p_o_ntcnx*(1-t)*c*(1-x) + p_o_ntcx*(1-t)*c*x
  o_expect <- o_expect + p_o_tncnx*t*(1-c)*(1-x) + p_o_tncx*t*(1-c)*x +
    p_o_tcnx*t*c*(1-x) + p_o_tcx*t*c*x
  o <- rbinom(n, 1, bound_to_margin(o_expect, prob_margin))
  
  a_expect <- p_a_nsnonx*(1-s)*(1-o)*(1-x) + p_a_nsnox*(1-s)*(1-o)*x +
    p_a_nsonx*(1-s)*o*(1-x) + p_a_nsox*(1-s)*o*x
  a_expect <- a_expect + p_a_snonx*s*(1-o)*(1-x) + p_a_snox*s*(1-o)*x +
    p_a_sonx*s*o*(1-x) + p_a_sox*s*o*x
  a <- rbinom(n, 1, bound_to_margin(a_expect, prob_margin))
  
  
  # p_o_t <- p_o_tcx*p_c*p_x*p_t + p_o_tcnx*p_c*(1-p_x)*p_t + p_o_tncx*(1-p_c)*p_x*p_t + p_o_tncnx*(1-p_c)*(1-p_x)*p_t
  # p_o_nt <- p_o_ntcx*p_c*p_x*(1-p_t) + p_o_ntcnx*p_c*(1-p_x)*(1-p_t) +  p_o_ntncx*(1-p_c)*p_x*(1-p_t) + p_o_ntncnx*(1-p_c)*(1-p_x)*(1-p_t)
  
  
  p_o_t <- (p_o_tcx + add_impact_s_o)*p_s_tcx*p_c*p_x + 
    (p_o_tcnx + add_impact_s_o)*p_s_tcnx*p_c*(1-p_x) + 
    (p_o_tncx + add_impact_s_o)*p_s_tncx*(1-p_c)*p_x + 
    (p_o_tncnx + add_impact_s_o)*p_s_tncnx*(1-p_c)*(1-p_x)
  p_o_t <- p_o_t + p_o_tcx*(1-p_s_tcx)*p_c*p_x + 
    p_o_tcnx*(1-p_s_tcnx)*p_c*(1-p_x) + 
    p_o_tncx*(1-p_s_tncx)*(1-p_c)*p_x + 
    p_o_tncnx*(1-p_s_tncnx)*(1-p_c)*(1-p_x)
  
  p_o_nt <- (p_o_ntcx + add_impact_s_o)*p_s_ntcx*p_c*p_x + 
    (p_o_ntcnx + add_impact_s_o)*p_s_ntcnx*p_c*(1-p_x) + 
    (p_o_ntncx + add_impact_s_o)*p_s_ntncx*(1-p_c)*p_x + 
    (p_o_ntncnx + add_impact_s_o)*p_s_ntncnx*(1-p_c)*(1-p_x)
  p_o_nt <- p_o_nt + p_o_ntcx*(1-p_s_ntcx)*p_c*p_x + 
    p_o_ntcnx*(1-p_s_ntcnx)*p_c*(1-p_x) + 
    p_o_ntncx*(1-p_s_ntncx)*(1-p_c)*p_x + 
    p_o_ntncnx*(1-p_s_ntncnx)*(1-p_c)*(1-p_x)
  
  y <- ifelse(a==1, o, NA)
  dt <- data.table(t, x, c, s, o, a, y)
  
  return(list(
    dt = dt,
    p_o_t = p_o_t,
    p_o_nt = p_o_nt
  ))
}


generate_data_simple <- function(
  n=2000, p_t=0.5, p_x=0.5,
  t_o=0.2, o_a=-0.2, x_o=0.2, x_a=0.1
){
  t <- rbinom(n, 1, p_t)
  x <- rbinom(n, 1, p_x)

  p_o <- 0.5 + (t_o*t - t_o*(1-t) + x_o*x - x_o*(1-x))/2
  o <- rbinom(n, 1, p_o)
  p_a <- 0.5 + (o_a*o - o_a*(1-o) + x_a*x - x_a*(1-x))/2
  a <- rbinom(n, 1, p_a)
  y <- ifelse(a==1, o, NA)
  
  dt <- data.table(t, x, c, s=1, o, a, y)
  
  p_o_t <- (0.5 + (t_o + x_o)/2)*p_x + (0.5 + (t_o - x_o)/2)*(1-p_x)
  p_o_nt <- (0.5 + (-t_o + x_o)/2)*p_x + (0.5 + (-t_o - x_o)/2)*(1-p_x)
  
  return(list(
    dt = dt,
    p_o_t = p_o_t,
    p_o_nt = p_o_nt
  ))
}

generate_data_high_dimensional <- function(
    n, m, p_t, p_x, p_c, corr_c, corr_x, sigma, 
    intercept_o, base_o, impact_t_o, 
    intercept_s, base_s, impact_t_s, 
    intercept_ca, intercept_s_ca, intercept_o_ca, 
    covariate_importance, add_impact_s_o, n_bins = 3,
    s_only_0 = FALSE, a_margin = 0.1
) {
  t <- rbinom(n, 1, p_t)
  x_vars <- matrix(rbinom(n*m, 1, p_x), nrow = n)
  colnames(x_vars) <- paste("x", 1:m, sep="")
  c_vars <- matrix(rbinom(n*m, 1, p_c), nrow = n)
  colnames(c_vars) <- paste("c", 1:m, sep="")
  
  Sigma <- matrix(rep(corr_c, 4*4), 4,  4)
  diag(Sigma) <- rep(1, 4)
  Sigma <- Sigma*sigma
  coefs_c <- mvrnorm(m, rep(0, 4), Sigma)/m
  coefs_c <- matrix(coefs_c, ncol = 4)
  
  Sigma <- matrix(rep(corr_x, 4*4), 4,  4)
  diag(Sigma) <- rep(1, 4)
  Sigma <- Sigma*sigma
  coefs_x <- mvrnorm(m, rep(0, 4), Sigma)/m
  coefs_x <- matrix(coefs_x, ncol = 4)
  coefs_x <- coefs_x[order(-abs(coefs_x[, 1])), , drop=FALSE]
  
  if(s_only_0){
    s <- 0
  }else{
    s_linear <- intercept_s - x_vars%*%coefs_x[, 2] - c_vars%*%coefs_c[, 2] +
      impact_t_s*t
    s <- as.numeric(runif(n) < logistic(s_linear)) 
  }
  
  o_linear <- intercept_o +
    impact_t_o*(s + t)*
      covariate_importance*(x_vars%*%coefs_x[, 1] + c_vars%*%coefs_c[, 1])/m +
    impact_t_o*(s + t) +
    covariate_importance*(x_vars%*%coefs_x[, 1] + c_vars%*%coefs_c[, 1])/m +
    add_impact_s_o*s
  o <- as.numeric(runif(n) < logistic(o_linear))
  
  base_ca <- intercept_ca + covariate_importance*(x_vars%*%coefs_x[, 4] + c_vars%*%coefs_c[, 4])/m
  a_o1s0_prop <- rep(1, n)
  a_o1s0 <- rep(1, n)
  a_o1s1_prop <- logistic(base_ca + intercept_s_ca + intercept_o_ca)[, 1] %>% 
    bin_values(a_margin, n_bins) %>% round(2)
  a_o1s1 <- as.numeric(runif(n) < a_o1s1_prop)
  a_o0s0_prop <- logistic(base_ca)[, 1] %>% 
    bin_values(a_margin, n_bins) %>% round(2)
  a_o0s0 <- as.numeric(runif(n) < a_o0s0_prop)
  a_o0s1_prop <- logistic(base_ca  + intercept_s_ca)[, 1] %>% 
    bin_values(a_margin, n_bins) %>% round(2)
  a_o0s1 <- as.numeric(runif(n) < a_o0s1_prop)
  
  a <- a_o0s0*(1-s)*(1-o) + a_o1s0*(1-s)*o +  a_o0s1*s*(1-o) +  a_o1s1*s*o
  y <- ifelse(a==1, o, NA)
  
  return( data.table(
    t, x_vars, c_vars, s, o, a, y, 
    a_o0s0_prop, a_o0s1_prop, a_o1s0_prop, a_o1s1_prop,
    a_o0s0, a_o0s1, a_o1s0, a_o1s1
  ) )
}

bin_values <- function(values, a_margin, n_bins=3){
  if( sd(values) < range_tolerance ){
    mode_value <- names(sort(-table(values)))[1] %>% as.numeric
    return(rep(mode_value, length(values)))
  }
  
  if( a_margin < 1/n_bins ){
    a_margin <- 1/(2*n_bins)
  }
  
  r <- range(values)

  probs_bins <- seq(0, 1, length.out = n_bins + 1)
  breaks <- quantile(values, probs = probs_bins, na.rm = TRUE)
  mids <- (head(breaks, -1) + tail(breaks, -1)) / 2 
  res <- cut(values,
             breaks = breaks,
             include.lowest = TRUE,
             labels = mids %>% round(2)) %>% as.character() %>% as.numeric()
  return(res)
}

calc_ate <- function(dt_res_t, col_name){
  ate <- dt_res_t[t==1, get(col_name)] - dt_res_t[t==0, get(col_name)]
  new_row <- data.table(t = 'delta', value = ate)
  setnames(new_row, 'value', col_name)
  dt_res_t <- rbind(dt_res_t, new_row)
  dt_res_t
}
