suppressMessages({
library(data.table)
library(magrittr)
library(ggplot2)
library(snowfall)
# source('code/utils.R')
# source('code/estimator.R')
})

set.seed(1234)

sims_n <- 100
increment <- 0.05
max_impact_t_o <- 0.4
log_file <- "log_sim_mnar.txt"
parallel <- T

# CONFIGURATION
# Example:
# x - age
# c - good genetics

ns <- c(2000, 5000)

# CONF 1
base_m <- 0
impact_x_m <- 0.5
impact_s_m <- 0.1
impact_no_m <- 0.4 # put 0 for MAR

p_t <- 0.5
p_c <- 0.5
p_x <- 0.5

base_s <- 0.1
impact_t_s <- 0.25
# impact_c_s <- -0.1
impact_c_s <- 0.15
impact_x_s <- 0.05

base_s <- 0.0
impact_t_s <- 0.0
impact_c_s <- 0.0
impact_x_s <- 0.0

base_o <- 0.2
impact_c_o <- 0.15
impact_x_o <- 0.25
# impact_x_o <- -0.15

# CONF 2
# base_m <- 0
# impact_x_m <- 0.05
# impact_s_m <- 0.0
# impact_no_m <- 0.4 # put 0 for MAR
#
# p_t <- 0.5
# p_c <- 0.5
# p_x <- 0.5
#
# base_s <- 0.0
# impact_t_s <- 0.0
# impact_c_s <- 0.0
# impact_x_s <- 0.0
#
# base_o <- 0.1
# impact_c_o <- 0.0
# impact_x_o <- 0.1

impacts_t_o <- seq(0, max_impact_t_o, increment)
# impacts_t_o <- seq(max_impact_t_o - increment, max_impact_t_o, increment)

capture.output({
  suppressMessages({
    sfInit(parallel = parallel, cpus = 10, slaveOutfile=log_file)
    sfExportAll()
    sfLibrary(data.table)
    sfLibrary(magrittr)
    invisible(sfLibrary(snowfall))
  })
}, file = nullfile())

results_list <- sfLapply(impacts_t_o, function(impact_t_o){
  results <- data.table()
  sfCat(paste("Impact ", impact_t_o), sep="\n")
  for(n in ns){
    for(sim_n in 1:sims_n){
      data_info <- generate_data(
        base_s=base_s,
        impact_t_s=impact_t_s,
        impact_x_s=impact_x_s,
        impact_c_s=impact_c_s,
        base_o=base_o,
        impact_t_o=impact_t_o,
        impact_x_o=impact_x_o,
        impact_c_o=impact_c_o,
        base_m=base_m,
        impact_no_m=impact_no_m,
        impact_x_m=impact_x_m,
        impact_s_m=impact_s_m,
        p_c=p_c,
        p_x=p_x,
        p_t=p_t,
        n=n
      )
      
      dt_full <- data_info$dt
      p_o_t <- data_info$p_o_t
      p_o_nt <- data_info$p_o_nt
      dt <- copy(dt_full)
      dt[, o :=NULL]
      dt[, c :=NULL]
      
      result_mar <- mar_estimator(dt, col_t='t', cols_w=c("x"))
      result_mnar <- mnar_estimator(dt, col_t='t', cols_w=c('x', 's'), apply_bounds=TRUE, estimator=NULL, n_sample=100)
      
      alt_estimator <- function(dt_w, val_t){
        dt_w[(a==1) & (t==val_t), mean(y)]
      }
      result_mnar_alt <- mnar_estimator(
        dt, col_t='t', cols_w=c('x', 's'), apply_bounds=TRUE, estimator=alt_estimator, n_sample=1000
      )
      setnames(result_mnar_alt, 'p_o_t_mnar', 'p_o_t_mnar_alt')
      result_mnar_alt[, p_o_t_lb:=NULL][, p_o_t_ub:=NULL]
      
      result <- data.table(t=c(1, 0), p_o_expected = c(p_o_t, p_o_nt)) %>%
        merge(dt_full[, .(sample_p_o=mean(o)), t], by="t") %>%
        merge(dt_full[, .(observed_p_o=mean(y, na.rm = TRUE)), t], by="t") %>%
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
suppressMessages(sfStop())

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
  .[, .(mean_error = mean(error), sd_error = sd(error)), .(method, impact_t_o, t, n)] %>%
  .[, .(
    mean_error,
    l_ci = mean_error - 1.96*sd_error,
    u_ci = mean_error + 1.96*sd_error,
    impact_t_o, method, t, n
  )] %>%
  ggplot(aes(x=impact_t_o, y=mean_error, color=method, group=method)) +
  geom_line() +
  geom_ribbon(aes(ymin = l_ci, ymax = u_ci, alpha = 0.001, fill=method, linetype=method)) +
  facet_grid(n~t) +
  scale_alpha(guide = 'none') ->
  res_plot
print(res_plot)
