library(emmeans)

MIN_PROB <- 1e-5

condition_on <- function(dt_, cols_w, val_w){
  dt_[Reduce("&", Map(function(col, val) get(col) == val, cols_w, val_w))]
}

rho <- function(
    val_o, val_t, col_t, dtp, p_w_t, delta=FALSE, p_t_global=NA, p_w
  ){
  p_t <- dtp[, mean(get(col_t)==val_t)]  
  dtp_a <- condition_on(dtp, "a", 1)
  pt_oa <- dtp_a[(y == val_o), mean(get(col_t)==val_t)]
  pt_noa <- dtp_a[(y != val_o), mean(get(col_t)==val_t)]
  
  if( delta ){
    if(is.na(pt_oa)){
      dt_tw <- condition_on(dtp, col_t, 1)
      lb_1 <- calc_p_o_t_lb(dt_tw, val_o)
      ub_1 <- lb_1 + calc_p_o_t_ub_2(dt_tw)
      
      dt_tw <- condition_on(dtp, col_t, 0)
      lb_0 <- calc_p_o_t_lb(dt_tw, val_o)
      ub_0 <- lb_0 + calc_p_o_t_ub_2(dt_tw)
      
      ub <- ub_1 - lb_0
      lb <- lb_1 - ub_0
      
      hat_p_o_t <- NA
      lb <- lb * p_w
      ub <- ub * p_w
    }else{
      hat_p_o_t_0 <- (pt_oa-p_t_global)/(p_t_global*(1-p_t_global))
      hat_p_o_t <- hat_p_o_t_0*(p_t - pt_noa)/(pt_oa - pt_noa)
      hat_p_o_t <- max(min(hat_p_o_t, 1), -1)
      
      lb <- hat_p_o_t_0*dtp_a[, mean(y)]*dtp[, mean(a)]
      lb <- max(-1, lb)
      ub <- lb + hat_p_o_t_0*(1-dtp[, mean(a)])
      ub <- min(1, ub)
      
      hat_p_o_t <- hat_p_o_t * p_w
      lb <- lb * p_w
      ub <- ub * p_w
    }
  }else{
    hat_p_o_t <- pt_oa*(p_t - pt_noa)/(p_t*(pt_oa - pt_noa))
    hat_p_o_t <- max(min(hat_p_o_t, 1), 0)
    
    dt_tw <- condition_on(dtp, col_t, val_t)
    
    lb <- calc_p_o_t_lb(dt_tw, val_o)
    ub <- lb + calc_p_o_t_ub_2(dt_tw)
    
    hat_p_o_t <- hat_p_o_t * p_w_t
    lb <- lb * p_w_t
    ub <- ub * p_w_t
  }
  
  if(!is.na(hat_p_o_t)){
    if( hat_p_o_t > ub ){
      hat_p_o_t <- ub
    }else if( hat_p_o_t < lb ){
      hat_p_o_t <- lb
    }
    ub <- hat_p_o_t
    lb <- hat_p_o_t
  }
  
  return(list(
    value=hat_p_o_t,
    ub=ub,
    lb=lb
  ))
  
}

calc_p_o_t_lb <- function(dtp, val_o){
  res <- condition_on(dtp, "a", 1)[, mean(y == val_o)] * dtp[, mean(a == 1)]
  if( is.na(res) ){
    res <- 0
  }
  return(res)
}

calc_p_o_t_ub_2 <- function(dtp){
  dtp[, mean(a == 0)]
}


calc_p_sign <- function(dt_w, val_o, n_sample=100){
  dt_wa <- condition_on(dt_w, "a", 1)
  dt_wao <- dt_wa[(y == val_o)]
  dt_wano <- dt_wa[(y != val_o)]
  
  delta_0 <- dt_wao[, mean(t)] - dt_wano[, mean(t)]
  
  deltas <- c()
  for(n_s in 1:n_sample){
    dt_wao_ <- dt_wao[sample(1:nrow(dt_wao), replace=TRUE)]
    dt_wano_ <- dt_wano[sample(1:nrow(dt_wano), replace=TRUE)]
    deltas <- c(deltas, dt_wao_[, mean(t)] - dt_wano_[, mean(t)])
  }
  
  q <- mean(na.omit(delta_0*deltas > 0))
  q <- max(1/2, q)
  return(q)
}

calc_phi_o_t <- function(
    dt_, val_o, val_t, col_t, cols_w, apply_bounds=TRUE, estimator=NULL, 
    n_sample=100){
  vals_t_w <- dt_[(get(col_t)==val_t), unique(.SD), .SDcols=cols_w]
  t_n <- condition_on(dt_, col_t, val_t)[, .N]
  
  phi_o_t <- 0
  phi_o_t_lb <- 0
  phi_o_t_ub <- 0
  prop_bounds <- 0
  for(row_n in 1:nrow(vals_t_w)){
    val_w <- vals_t_w[row_n, ]
    
    dt_w <- condition_on(dt_, cols_w, val_w)
    p_w_t <- condition_on(dt_w, col_t, val_t)[, .N]/t_n
    estimate <- rho(val_o, val_t, col_t, dt_w, p_w_t)
    hat_p_o_t_w <- estimate$value
    
    if(!is.null(estimator)){
      new_estimator <- estimator(dt_w, val_t)
      if(is.na(hat_p_o_t_w)){
        hat_p_o_t_w <- new_estimator
      }else{
        p_sign <- calc_p_sign(dt_w, val_o, n_sample=n_sample)
        hat_p_o_t_w <- hat_p_o_t_w*(2*p_sign - 1) + 2*(1-p_sign)*new_estimator
      }
    }
    
    if(apply_bounds){
      hat_p_o_t_w_lb <- estimate$lb
      hat_p_o_t_w_ub <- estimate$ub
      
      if(is.na(hat_p_o_t_w)){
        hat_p_o_t_w <- (hat_p_o_t_w_lb + hat_p_o_t_w_ub)/2
        prop_bounds <- prop_bounds + p_w_t
      }else{
        if( hat_p_o_t_w > hat_p_o_t_w_ub){
          hat_p_o_t_w <- hat_p_o_t_w_ub
          prop_bounds <- prop_bounds + p_w_t
        }else if(hat_p_o_t_w < hat_p_o_t_w_lb){
          hat_p_o_t_w <- hat_p_o_t_w_lb
          prop_bounds <- prop_bounds + p_w_t
        }
      }
      phi_o_t_lb <- phi_o_t_lb + hat_p_o_t_w_lb
      phi_o_t_ub <- phi_o_t_ub + hat_p_o_t_w_ub
    }
    
    phi_o_t <- phi_o_t + hat_p_o_t_w
  }
  
  if(!apply_bounds){
    prop_bounds <- NA
    phi_o_t_lb <- NA
    phi_o_t_ub <- NA
  }
  
  return(list(
    phi_o_t = phi_o_t,
    phi_o_t_lb = phi_o_t_lb,
    phi_o_t_ub = phi_o_t_ub,
    prop_bounds = prop_bounds
  ))
}

calc_phi_o_t_delta <- function(dt_, val_o, col_t, cols_w, apply_bounds=TRUE, estimator=NULL, n_sample=100){
  val_t <- 1
  vals_t_w <- dt_[(get(col_t)==val_t), unique(.SD), .SDcols=cols_w]
  t_n <- condition_on(dt_, col_t, val_t)[, .N]
  n <- dt_[, .N]
  p_t_global <- t_n/n
  
  phi_o_t <- 0
  phi_o_t_lb <- 0
  phi_o_t_ub <- 0
  prop_bounds <- 0
  for(row_n in 1:nrow(vals_t_w)){
    val_w <- vals_t_w[row_n, ]
    
    dt_w <- condition_on(dt_, cols_w, val_w)
    p_w <- dt_w[, .N]/n
    p_w_t <- condition_on(dt_w, col_t, val_t)[, .N]/t_n
    estimate <- rho(val_o, val_t, col_t, dt_w, p_w_t, delta=TRUE, p_t_global, p_w)
    hat_p_o_t_w <- estimate$value
    
    if(!is.null(estimator)){
      new_estimator <- estimator(dt_w, val_t=1) - estimator(dt_w, val_t=0)
      if(is.na(hat_p_o_t_w)){
        hat_p_o_t_w <- new_estimator
      }else{
        p_sign <- calc_p_sign(dt_w, val_o, n_sample=n_sample)
        hat_p_o_t_w <- hat_p_o_t_w*(2*p_sign - 1) + 2*(1-p_sign)*new_estimator
      }
    }

    if(apply_bounds){
      hat_p_o_t_w_lb <- estimate$lb
      hat_p_o_t_w_ub <- estimate$ub
      
      if(is.na(hat_p_o_t_w)){
        hat_p_o_t_w <- (hat_p_o_t_w_lb + hat_p_o_t_w_ub)/2
        prop_bounds <- prop_bounds + dt_w[, .N]/n
      }else{
        if( hat_p_o_t_w > hat_p_o_t_w_ub){
          hat_p_o_t_w <- hat_p_o_t_w_ub
          prop_bounds <- prop_bounds + dt_w[, .N]/n
        }else if(hat_p_o_t_w < hat_p_o_t_w_lb){
          hat_p_o_t_w <- hat_p_o_t_w_lb
          prop_bounds <- prop_bounds + dt_w[, .N]/n
        }
      }
      phi_o_t_lb <- phi_o_t_lb + hat_p_o_t_w_lb
      phi_o_t_ub <- phi_o_t_ub + hat_p_o_t_w_ub
    }
    
    phi_o_t <- phi_o_t + hat_p_o_t_w
  }
  
  if(!apply_bounds){
    prop_bounds <- NA
    phi_o_t_lb <- NA
    phi_o_t_ub <- NA
  }
  
  return(list(
    phi_o_t = phi_o_t,
    phi_o_t_lb = phi_o_t_lb,
    phi_o_t_ub = phi_o_t_ub,
    prop_bounds = prop_bounds
  ))
}

calc_logodds_linear_models <- function(dt_, col_t, cols_w){
  
  interaction <- paste(c(col_t, cols_w), collapse = "*")
  formula_ <- paste("y~", interaction, sep="")
  model <- glm(formula_, family=binomial, data = dt_) 
  
  interaction <- paste(cols_w, collapse = " : ")
  formula_ <- paste("~", col_t, "|", interaction, sep=" ")
  contrasts_ <- pairs(emmeans(model, formula(formula_)))
  
  logodds_results <- data.frame(contrasts_@grid, summary(contrasts_)$estimate)
  colnames(logodds_results)[ncol(logodds_results)] <- "logodds"
  logodds_results <- merge(dt_[, .(prob=.N/dt_[, .N]), cols_w], logodds_results, by=cols_w)
  logodds <- logodds_results[, prob*logodds] %>% sum
  
  return(logodds)
}

calc_odds <- function(dt_w, col_t){
  p_t_o <- dt_w[y==1, mean(get(col_t))]
  p_t_no <- dt_w[y==0, mean(get(col_t))]
  if( is.na(p_t_o) & is.na(p_t_no)){
    return(NA)
  }
  if( p_t_o == 1 ){
    p_t_o <- 1 - MIN_PROB
    p_t_no <- MIN_PROB
  }else if( p_t_o == 0 ){
    p_t_o <- MIN_PROB
    p_t_no <- 1 - MIN_PROB
  }
  odds <- (p_t_o/(1-p_t_o)) / (p_t_no/(1-p_t_no))
  return(odds)
}

calc_logodds <- function(dt_, col_t, cols_w){
  results <- dt_[, .(odds=calc_odds(.SD, col_t), prob_w=.N/dt_[, .N]), cols_w]
  logodds <- -results[, prob_w*log(odds)] %>% sum
  return(logodds)
}

mnar_estimator <- function(
    dt_, col_t, cols_w, 
    apply_bounds=TRUE, estimator=NULL, n_sample=100,
    delta=TRUE
  ){
  
  val_o <- 1
  vals_t <- dt_[, unique(get(col_t))]
  
  res_t <- c()
  res_phi <- c()
  res_phi_ub <- c()
  res_phi_lb <- c()
  res_prop_bounds <- c()
  for(val_t in vals_t){
    res_ <- calc_phi_o_t(
      dt_, val_o, val_t, col_t, cols_w, 
      apply_bounds=apply_bounds, estimator=estimator, n_sample=n_sample
    )
    res_t <- c(res_t, val_t)
    res_phi <- c(res_phi, res_$phi_o_t)
    res_phi_lb <- c(res_phi_lb, res_$phi_o_t_lb)
    res_phi_ub <- c(res_phi_ub, res_$phi_o_t_ub)
    res_prop_bounds <- c(res_prop_bounds, res_$prop_bounds)
  }
  
  if( delta ){
    res_ <- calc_phi_o_t_delta(
      dt_, val_o, col_t, cols_w, 
      apply_bounds=apply_bounds, estimator=estimator, n_sample=n_sample
    )
    res_t <- c(res_t, "delta")
    res_phi <- c(res_phi, res_$phi_o_t)
    res_phi_lb <- c(res_phi_lb, res_$phi_o_t_lb)
    res_phi_ub <- c(res_phi_ub, res_$phi_o_t_ub)
    res_prop_bounds <- c(res_prop_bounds, res_$prop_bounds)
  }
  
  results <- data.table(
    p_o_t_mnar=res_phi, p_o_t_lb=res_phi_lb, p_o_t_ub=res_phi_ub,
    prop_bounds=res_prop_bounds
  )
  results[, c(col_t):= res_t]
  
  if(delta){
    lo_table <- data.table(p_o_t_mnar = calc_logodds_linear_models(dt_, col_t, cols_w))
    lo_table[, c(col_t) := 'logodds']
    results <- rbind(results, lo_table, fill=TRUE)
  }
  
  if( "delta" %in% results[, get(col_t)]){
    if( results[t=='delta', is.na(p_o_t_mnar)]){
      results[t=='delta', p_o_t_mnar := results[t==1, p_o_t_mnar] - results[t==0, p_o_t_mnar]]
    }
  }
  return(results)
}


mar_estimator <- function(dt_, col_t, cols_w){
  p_y_tsxa <- dt_[a == 1, .(p_y=mean(y)), c(col_t, cols_w) ]
  p_sx_t0 <- dt_[get(col_t)==0, .(p_sx = .N/dt_[get(col_t)==0, .N], t=0), cols_w]
  p_sx_t1 <- dt_[get(col_t)==1, .(p_sx = .N/dt_[get(col_t)==1, .N], t=1), cols_w]
  p_sx_t <- rbind(p_sx_t0, p_sx_t1)
  if( col_t != "t"){
    p_sx_t[, c(col_t) := t][, t := NULL]
  }
  res <- merge(p_y_tsxa, p_sx_t, by = c(col_t, cols_w) )
  
  return(res[, .(p_o_m=p_y*p_sx, t=get(col_t))][, .(p_o_t_mar=sum(p_o_m)), t])
}

select_cov <- function(dt, covariates){
  df <- as.data.frame(cbind(
    y=dt[a==1, y], 
    dt[a==1, .SD, .SDcols = covariates]
  ))
  model <- rpart(y~., data=df,cp =0, minsplit=2, maxdepth=1)
  return(model$frame$var[1])
}

select_cov_ <- function(dt, X){
  inds <- dt[, a == 1]
  df <- as.data.frame(cbind(y=dt[a==1, y], X[inds]))
  model <- rpart(y~., data=df,cp =0, minsplit=2, maxdepth=1)
  return(model$frame$var[1])
}

mnar_estimator_covariate <- function(
    dt_, col_t, col_s, cols_ca, col_x, 
    apply_bounds=TRUE, estimator=NULL, n_sample=100, delta=TRUE
){
  vals_t <- dt_[, unique(get(col_t))]
  
  result <- list()
  for(val_t in vals_t){
    dt_t <- condition_on(dt_, 't', val_t)
    phi <- mnar_estimator(
      dt_t, col_t=col_x, cols_w=c('s', cols_ca), 
      apply_bounds=TRUE, delta=FALSE
    )
    
    p_x_t <- dt_t[, .(p=.N/dt_t[, .N]), col_x]
    res <- merge(phi, p_x_t, by=col_x) %>% 
      .[, .SD*p, .SDcols = c('p_o_t_mnar', 'p_o_t_lb', 'p_o_t_ub')] %>% 
      .[, lapply(.SD, sum)] %>% 
      .[, c(col_t) := val_t]
    result[[val_t %>% as.character]] <- res
  }
  result <- rbindlist(result)
  
  colnames(result) <- c(
    "p_o_t_mnar_cov", "p_o_t_lb_cov", "p_o_t_ub_cov", "t"
  )
  result[, t := as.character(t)]
  
  if( delta ){
    delta_result <- data.table(
      p_o_t_mnar_cov = result[t==1, p_o_t_mnar_cov] - result[t==0, p_o_t_mnar_cov],
      p_o_t_ub_cov = result[t==1, p_o_t_ub_cov] - result[t==0, p_o_t_lb_cov],
      p_o_t_lb_cov = result[t==1, p_o_t_lb_cov] - result[t==0, p_o_t_ub_cov],
      t = "delta"
    )
    result <- rbind(result, delta_result)
  }

  return(result)
}


