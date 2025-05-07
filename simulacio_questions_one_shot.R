library(data.table)
library(magrittr)
library(ggplot2)
library(snowfall)
library(rpart)
library(glue)
# source('code/utils.R')
# source('code/estimator.R')

# set.seed(1234)

sims_n <- 100
n <- 2000
m <- 50
log_file <- "log_sim_questions.txt"

p_t <- 0.5
p_c <- 0.5
p_x <- 0.5

intercept_o <- 0
base_o <- 0.2
impact_t_o <- 2
# intercept_s <- -10000
intercept_s <- 0
base_s <- 0.2
impact_t_s <- 0.4
# intercept_ca <- -5
intercept_s_ca <- 8
intercept_o_ca <- 8
intercept_ca <- 0
# intercept_s_ca <- -0.5
# intercept_o_ca <- 1.5
corr_c <- 0.8
corr_x <- 0.8
sigma <- 100

cols_ca <- c('a_o0s0_prop', 'a_o0s1_prop', 'a_o1s0_prop', 'a_o1s1_prop')
cols_ca_ <- c('a_o0s0', 'a_o0s1', 'a_o1s0', 'a_o1s1')
cols_ca_0 <- c('a_o0s0_prop', 'a_o1s0_prop')
cols_ca_1 <- c('a_o0s1_prop', 'a_o1s1_prop')

dt_full <- generate_data_high_dimensional(
    n, m, p_t, p_x, p_c, corr_c, corr_x, sigma, 
    intercept_o, base_o, impact_t_o, 
    intercept_s, base_s, impact_t_s, 
    intercept_ca, intercept_s_ca, intercept_o_ca,
    s_only_0 = FALSE
) 

# dt <- dt_full[, .(
#   t, a, y, s, a_o0s0_prop, a_o0s1_prop, a_o1s0_prop, a_o1s1_prop
# )]

dt <- dt_full[, .(
  t, a, y, s, a_o0s0_prop, a_o1s0_prop
)]

X <- dt_full[, .SD, .SDcols = (names(dt_full) %like% "x")]
dt <- cbind(
  dt,
  X  
)

col_x <- select_cov(dt, X)

result_mnar <- mnar_estimator(
  dt, col_t='t', cols_w=c(cols_ca_0, 's'), 
  estimator=NULL, n_sample=100,
  delta=TRUE
)

result_mnar_covariate <- mnar_estimator_covariate(
  dt, col_t='t', col_s='s', 
  cols_ca = cols_ca_0, 
  col_x = col_x
)

dt[, mean(a)]

true_o <- dt_full[, .(value=mean(o)), t]
delta <- true_o[t==1, value] - true_o[t==0, value]
results <- rbind(
  true_o,
  data.table(t="delta", value = delta)
) %>% merge(result_mnar) %>% 
  merge(result_mnar_covariate, all.x=TRUE)

results[t=='delta', abs(p_o_t_mnar-value)] %>% print
results[t=='delta', abs(p_o_t_mnar_cov-value)] %>% print














