
library(data.table)
library(magrittr)
library(ggplot2)
library(snowfall)
source('utils.R')
source('estimators.R')

sims_n <- 100
increment <- 0.05
max_impact_t_o <- 0.3

log_file <- "log_sim_mnar_simple.txt"
parallel <- T

# CONFIGURATION

ns <- c(1000, 2000, 5000)

n=2000
p_t=0.5
p_x=0.5
o_a=-0.2 
x_o=0.2 
x_a=0.1
impacts_t_o <- seq(0, max_impact_t_o, increment)

sfInit(parallel = parallel, cpus = 10, slaveOutfile=log_file)
sfExportAll()
sfLibrary(data.table)
sfLibrary(magrittr)
invisible(sfLibrary(snowfall))

results_list <- sfLapply(impacts_t_o, function(impact_t_o){
  results <- data.table()
  sfCat(paste("Impact ", impact_t_o), sep="\n")
  for(n in ns){
    for(sim_n in 1:sims_n){
      data_info <- generate_data_simple(
        n=n, p_t=p_t, p_x=p_x,
        t_o=impact_t_o, o_a=o_a, x_o=x_o, x_a=x_a
      )
      
      dt_full <- data_info$dt
      p_o_t <- data_info$p_o_t
      p_o_nt <- data_info$p_o_nt
      dt <- copy(dt_full)
      dt[, o :=NULL]
      dt[, c :=NULL]
      
      result_mar <- mar_estimator(dt, col_t='t', cols_w=c("x"))
      result_mnar <- mnar_estimator(
        dt, col_t='t', cols_w=c('x', 's'), 
        apply_bounds=TRUE, estimator=NULL, n_sample=100,
        delta=TRUE
      )
      
      alt_estimator <- function(dt_w, val_t){
        dt_w[(a==1) & (t==val_t), mean(y)]
      }
      result_mnar_alt <- mnar_estimator(
        dt, col_t='t', cols_w=c('x', 's'), apply_bounds=TRUE, 
        estimator=alt_estimator, n_sample=1000,
        delta=TRUE
      )
      setnames(result_mnar_alt, 'p_o_t_mnar', 'p_o_t_mnar_alt')
      setnames(result_mnar_alt, 'prop_bounds', 'prop_bounds_alt')
      result_mnar_alt[, p_o_t_lb:=NULL][, p_o_t_ub:=NULL]
      
      missings_prop <- dt[, .(prop_miss=mean(a==1)), t]
      new_row <- data.table(t = 'delta', prop_miss = mean(dt[, (a==1)]))
      missings_prop <- rbind(missings_prop, new_row)
      
      dt_full_res <- calc_ate(
        dt_full[, .(sample_p_o=mean(o)), t], "sample_p_o"
      )
      
      dt_full_res_obs <- calc_ate(
        dt_full[, .(observed_p_o=mean(y, na.rm = TRUE)), t], "observed_p_o"
      )
      
      result_mar <- calc_ate(
        result_mar, "p_o_t_mar"
      )
      
      result <- data.table(
          t=c(1, 0, 'delta'), p_o_expected = c(p_o_t, p_o_nt, p_o_t - p_o_nt)
        ) %>%
        merge(missings_prop, by="t") %>% 
        merge(dt_full_res, by="t") %>%
        merge(dt_full_res_obs, by="t") %>%
        merge(result_mar, by="t") %>%
        merge(result_mnar, by="t") %>%
        merge(result_mnar_alt, by="t")
      result$impact_t_o <- impact_t_o
      result$n <- n
      
      results <- rbind(results, result)
    }
  }
  return(results)
})
sfStop()

do.call(rbind, results_list) %>%
  .[, .(
    mnar = abs(p_o_t_mnar - sample_p_o),
    mnar_smoothed = abs(p_o_t_mnar_alt - sample_p_o),
    mar = abs(p_o_t_mar - sample_p_o),
    observed = abs(observed_p_o - sample_p_o),
    impact_t_o,
    t, n
  )] %>%
  melt(id.vars=c("impact_t_o", "t", "n"), variable.name = "method", value.name = "error") %>%
  .[, .(
    mean_error = mean(error), 
    l_ci = quantile(error, probs=0.025),
    u_ci = quantile(error, probs=0.975)
  ), .(method, impact_t_o, t, n)] %>%
  ggplot(aes(x=impact_t_o, y=mean_error, color=method, group=method)) +
  geom_line() +
  geom_ribbon(
    aes(ymin = l_ci, ymax = u_ci, fill=method, linetype=method), alpha=0.3
  ) +
  facet_grid(n~t) +
  scale_alpha(guide = 'none') ->
  res_plot
print(res_plot)

do.call(rbind, results_list) %>%
  .[, .(
    impact_t_o, t, n, prop_miss, prop_bounds
  )] %>% 
  melt(id.vars=c("impact_t_o", "t", "n"), variable.name = "type", value.name = "proportion") %>%
  .[, .(
    proportion = mean(proportion), 
    l_ci = quantile(proportion, probs=0.025),
    u_ci = quantile(proportion, probs=0.975)
  ), .(type, impact_t_o, t, n)] %>%
  ggplot(aes(x=impact_t_o, y=proportion, color=type, group=type)) +
  geom_line() +
  geom_ribbon(
    aes(ymin = l_ci, ymax = u_ci, fill=type, linetype=type), alpha=0.3
  ) +
  facet_grid(n~t) +
  scale_alpha(guide = 'none') ->
  res_plot_miss
print(res_plot_miss)
