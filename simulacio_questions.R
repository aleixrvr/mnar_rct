suppressMessages({
  library(data.table)
  library(magrittr)
  library(ggplot2)
  library(snowfall)
  library(rpart)
  # source('code/utils.R')
  # source('code/estimator.R')
})

set.seed(1234)

sims_n <- 500
increment <- 0.2
max_impact_t_o <- 3
m <- 5
log_file <- "log_sim_mnar_questions.txt"
parallel <- F

cols_ca_0 <- c('a_o0s0_prop', 'a_o1s0_prop')

ns <- c(2000, 5000)

p_t <- 0.5
p_c <- 0.5
p_x <- 0.5

intercept_o <- 0
base_o <- 0.2
impact_t_o <- 2
intercept_s <- 0
base_s <- 0.2
impact_t_s <- 0.4
intercept_s_ca <- 8
intercept_o_ca <- 8
intercept_ca <- 0
corr_c <- 0.8
corr_x <- 0.8
sigma <- 100
covariate_importance <- 2.5

impacts_t_o <- seq(0, max_impact_t_o, increment)

capture.output({
  suppressMessages({
    sfInit(parallel = parallel, cpus = 10, slaveOutfile=log_file)
    sfExportAll()
    sfLibrary(data.table)
    sfLibrary(magrittr)
    sfLibrary(MASS)
    sfLibrary(rpart)
    invisible(sfLibrary(snowfall))
  })
}, file = nullfile())

results_list <- sfLapply(impacts_t_o, function(impact_t_o){
  results <- data.table()
  sfCat(paste("Impact ", impact_t_o), sep="\n")
  for(n in ns){
    for(sim_n in 1:sims_n){
      dt_full <- generate_data_high_dimensional(
        n, m, p_t, p_x, p_c, corr_c, corr_x, sigma, 
        intercept_o, base_o, impact_t_o, 
        intercept_s, base_s, impact_t_s, 
        intercept_ca, intercept_s_ca, intercept_o_ca,
        covariate_importance, s_only_0 = FALSE
      ) 
      
      dt <- dt_full[, .(
        t, a, y, s, a_o0s0_prop, a_o1s0_prop
      )]
      
      result_mnar <- mnar_estimator(
        dt, col_t='t', cols_w=c(cols_ca_0, 's'), 
        estimator=NULL, n_sample=100,
        delta=TRUE
      )
      
      X <- dt_full[, .SD, .SDcols = (names(dt_full) %like% "x")]
      dt <- cbind(
        dt,
        X  
      )
      
      tryCatch({
        col_x <- select_cov(dt, X)
      }, error = function(e) {
        browser()
      })
        
      result_mnar_covariate <- mnar_estimator_covariate(
        dt, col_t='t', col_s='s', 
        cols_ca = cols_ca_0, 
        col_x = col_x
      )
      
      true_o <- dt_full[, .(p_o_t=mean(o)), t]
      delta <- true_o[t==1, p_o_t] - true_o[t==0, p_o_t]
      result <- rbind(
        true_o,
        data.table(t="delta", p_o_t = delta)
      ) %>% merge(result_mnar) %>% 
        merge(result_mnar_covariate, all.x=TRUE)
      
      result$impact_t_o <- impact_t_o
      result$n <- n
      
      results <- rbind(results, result)
    }
  }
  return(results)
})
suppressMessages(sfStop())

do.call(rbind, results_list) %>%
  .[, .(
    mnar = abs(p_o_t_mnar - p_o_t),
    mnar_covariates = abs(p_o_t_mnar_cov - p_o_t),
    impact_t_o,
    t, n
  )] %>%
  melt(id.vars=c("impact_t_o", "t", "n"), variable.name = "method", value.name = "error") %>%
  .[, .(mean_error = mean(error), l_ci = quantile(error, 0.025), u_ci = quantile(error, 0.975)), .(method, impact_t_o, t, n)] %>%
  ggplot(aes(x=impact_t_o, y=mean_error, color=method, group=method)) +
  geom_line() +
  geom_ribbon(aes(ymin = l_ci, ymax = u_ci, alpha = 0.001, fill=method, linetype=method)) +
  facet_grid(n~t) +
  scale_alpha(guide = 'none') ->
  res_plot
print(res_plot)
