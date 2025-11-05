set.seed(1234)

sims_n <- 200
increment <- 0.2
max_impact_t_o <- 2
log_file <- "log_sim_mnar_questions.txt"
parallel <- T
outcome_tolerance_percentage <- 0.025

cols_ca <- c('a_o0s0_prop', 'a_o1s0_prop', 'a_o0s1_prop', 'a_o1s1_prop')

ns <- c(2000, 5000, 10000)

impacts_t_o <- seq(0, max_impact_t_o, increment)

capture.output({
  suppressMessages({
    sfInit(parallel = parallel, cpus = 10, slaveOutfile=log_file)
    sfExportAll()
    sfLibrary(data.table)
    sfLibrary(magrittr)
    sfLibrary(MASS)
    sfLibrary(rpart)
    sfLibrary(emmeans)
    invisible(sfLibrary(snowfall))
  })
}, file = nullfile())

results_list <- sfLapply(impacts_t_o, function(impact_t_o){
# for(impact_t_o in impacts_t_o){
  results <- data.table()
  sfCat(paste("Impact ", impact_t_o), sep="\n")
  for(n in ns){
    for(sim_n in 1:sims_n){
      dt_full <- generate_data_high_dimensional(
        n, m, p_t, p_x, p_c, corr_c, corr_x, sigma, 
        intercept_o, base_o, impact_t_o, 
        intercept_s, base_s, impact_t_s, 
        intercept_ca, intercept_s_ca, intercept_o_ca,
        covariate_importance, add_impact_s_o, s_only_0 = FALSE
      ) 
      
      min_outcome_prop <- dt_full[a==1, y %>% table %>% prop.table %>% min]
      diff_outcome <- dt_full[a==1, uniqueN(y)] != 1
      if( (min_outcome_prop > outcome_tolerance_percentage) & diff_outcome){
        dt <- dt_full[, .SD, .SDcols = c('t', 'a', 'y', 's', cols_ca)]
        
        result_mnar <- mnar_estimator(
          dt, col_t='t', cols_w=c(cols_ca, 's'), 
          estimator=NULL, n_sample=100,
          delta=TRUE
        )
        
        alt_estimator <- function(dt_w, val_t){
          dt_w[(a==1) & (t==val_t), mean(y)]
        }
        result_mnar_alt <- mnar_estimator(
          dt, col_t='t', cols_w=c(cols_ca, 's'), apply_bounds=TRUE, 
          estimator=alt_estimator, n_sample=100, delta=TRUE
        )
        result_mnar_alt[t=="delta", p_o_t_mnar := result_mnar[t==1, p_o_t_mnar] - result_mnar[t==0, p_o_t_mnar]]
        result_mnar_alt[t=='logodds', p_o_t_mnar := NA]
        
        setnames(result_mnar_alt, 'p_o_t_mnar', 'p_o_t_mnar_alt')
        setnames(result_mnar_alt, 'p_o_t_lb', 'p_o_t_lb_alt')
        setnames(result_mnar_alt, 'p_o_t_ub', 'p_o_t_ub_alt')
        result_mnar_alt[, prop_bounds:=NULL]
        
        sampled <- dt_full[, .(sample_p_o=mean(o)), .(t=as.character(t))]
        sampled <- rbind(sampled, list(
          't' = 'delta', sample_p_o = sampled[t == 1, sample_p_o] - sampled[t == 0, sample_p_o]
        ))
        dt_sm <- copy(dt_full)
        dt_sm[, y :=o]
        logodds_ratio <- calc_logodds_linear_models(dt_sm, col_t='t', cols_w=c(cols_ca, 's'))
        sampled <- rbind(sampled, list(
          't' = 'logodds', sample_p_o = logodds_ratio
        ))
        
        observed <- dt[, .(observed_p_o=mean(y, na.rm = TRUE)), t] %>% 
          .[, t := as.character(t)] 
        observed <- rbind(observed, data.table(
          t = 'delta', observed_p_o = observed[t == 1, observed_p_o] - observed[t == 0, observed_p_o]
        ))
        missings <- dt[, .(missing=mean(is.na(y))), t] %>% 
          .[, t := as.character(t)]
        
        result <- sampled %>%
          merge(observed, by="t", all.x=TRUE) %>%
          merge(result_mnar) %>% 
          merge(result_mnar_alt, all.x=TRUE) %>% 
          merge(missings, by="t", all.x=TRUE)
        
        result$impact_t_o <- impact_t_o
        result$n <- n
        
        results <- rbind(results, result)
      }
    }
  }
  return(results)
# }
})
suppressMessages(sfStop())

results <- do.call(rbind, results_list)
results[t=='0', t:='control']
results[t=='1', t:='treatment']
results[t=='delta', t:='ATE']
results[t=='logodds', t:='avg-adj-logodds']

results %>%
  .[t %chin% c("control", "treatment", 'ATE')] %>%
  .[, .(
    mnar = p_o_t_mnar - sample_p_o,
    mnar_smoothed = p_o_t_mnar_alt - sample_p_o,
    observed = observed_p_o - sample_p_o,
    impact_t_o,
    t, n
  )] %>%
  melt(id.vars=c("impact_t_o", "t", "n"), variable.name = "method", value.name = "error") %>%
  .[, .(
    mean_error = mean(error),
    l_ci = quantile(error, probs=0.025, na.rm=TRUE),
    u_ci = quantile(error, probs=0.975, na.rm=TRUE)),
    .(method, impact_t_o, t, n)
  ] %>% 
  ggplot(aes(x = impact_t_o, y = mean_error, linetype = method, fill = method, group = method)) +
  geom_ribbon(aes(ymin = l_ci, ymax = u_ci), alpha = 0.4) +
  geom_line() +
  geom_line(aes(x = impact_t_o, y=u_ci,linetype = method), alpha=0.2) +
  geom_line(aes(x = impact_t_o, y=l_ci,linetype = method), alpha=0.2) +
  facet_grid(n~t) +
  ylab("Bias") + 
  xlab("Impact from T to O") + 
  scale_linetype_manual(
    values = c("mnar" = "solid", "mnar_smoothed" = "dashed", "observed" = "dotted"),
    labels = c(
      "mnar" = "MNAR estimator", 
      "mnar_smoothed" = "Smoothed MNAR \n estimator", 
      "observed" = "Observed")
  ) +
  scale_fill_manual(
    values = c(
      "mnar" = "grey80",
      "mnar_smoothed" = "grey40",
      "observed" = "grey60"
    ),
    labels = c(
      "mnar" = "MNAR estimator", 
      "mnar_smoothed" = "Smoothed MNAR \n estimator", 
      "observed" = "Observed"
    )
  ) +
  scale_alpha(guide = "none") +
  guides(
    fill = guide_legend(
      title = "method",
      override.aes = list(
        linetype = c("solid", "dashed", "dotted"),
        fill = c("grey80", "grey40", "grey60")
      )
    )
  )  ->
  res_plot_1

results %>%
  .[!(t %chin% c("control", "treatment", 'ATE'))] %>%
  .[, .(
    mnar = p_o_t_mnar - sample_p_o,
    mnar_smoothed = p_o_t_mnar_alt - sample_p_o,
    observed = observed_p_o - sample_p_o,
    impact_t_o,
    t, n
  )] %>%
  melt(id.vars=c("impact_t_o", "t", "n"), variable.name = "method", value.name = "error") %>%
  .[, .(
    mean_error = mean(error),
    l_ci = quantile(error, probs=0.025, na.rm=TRUE),
    u_ci = quantile(error, probs=0.975, na.rm=TRUE)),
    .(method, impact_t_o, t, n)
  ] %>% 
  ggplot(aes(x = impact_t_o, y = mean_error, linetype = method, fill = method, group = method)) +
  geom_ribbon(aes(ymin = l_ci, ymax = u_ci), alpha = 0.4) +
  geom_line() +
  geom_line(aes(x = impact_t_o, y=u_ci,linetype = method), alpha=0.2) +
  geom_line(aes(x = impact_t_o, y=l_ci,linetype = method), alpha=0.2) +
  facet_grid(n~t) +
  ylab("Bias") + 
  xlab("Impact from T to O") + 
  scale_linetype_manual(
    values = c("mnar" = "solid", "mnar_smoothed" = "dashed", "observed" = "dotted"),
    labels = c(
      "mnar" = "MNAR estimator", 
      "mnar_smoothed" = "Smoothed MNAR \n estimator", 
      "observed" = "Observed")
  ) +
  scale_fill_manual(
    values = c(
      "mnar" = "grey80",
      "mnar_smoothed" = "grey40",
      "observed" = "grey60"
    ),
    labels = c(
      "mnar" = "MNAR estimator", 
      "mnar_smoothed" = "Smoothed MNAR \n estimator", 
      "observed" = "Observed"
    )
  ) +
  scale_alpha(guide = "none") +
  guides(
    fill = guide_legend(
      title = "method",
      override.aes = list(
        linetype = c("solid", "dashed", "dotted"),
        fill = c("grey80", "grey40", "grey60")
      )
    )
  )  ->
  res_plot_2

results %>%
  .[t %chin% c("control", "treatment")] %>%
  .[, .(
    bounds_range = p_o_t_ub - p_o_t_lb,
    impact_t_o,
    t, n=as.character(n)
  )] %>%
  melt(id.vars=c("impact_t_o", "t", "n"), variable.name = "quantity", value.name = "proportion") %>%
  .[, .(
    proportion = mean(proportion, na.rm=TRUE)),
    .(quantity, impact_t_o, n)
  ] %>%
  ggplot(aes(x=impact_t_o, y=proportion, linetype=n, group=n)) +
  geom_line() +
  ylab("Bound Range") + 
  xlab("Impact from T to O") + 
  ylim(c(0, NA))  ->
  res_plot_bounds

results %>%
  .[t %chin% c("control", "treatment")] %>%
  .[, .(
    prop_missing = missing,
    impact_t_o,
    t, n=as.character(n)
  )] %>%
  melt(id.vars=c("impact_t_o", "t", "n"), variable.name = "quantity", value.name = "proportion") %>%
  .[, .(
    proportion = mean(proportion, na.rm=TRUE)),
    .(quantity, impact_t_o, t, n)
  ] %>%
  ggplot(aes(x=impact_t_o, y=proportion, linetype=n, group=n)) +
  geom_line() +
  ylab("Proportion of missing data") + 
  xlab("Impact from T to O") + 
  ylim(c(0, NA)) +
  facet_grid(.~t, scales="free_y") ->
  res_plot_missing

ggsave(
  plot_name, 
  ((res_plot_1 + res_plot_2 + plot_layout(widths = c(2, 1))) / (res_plot_bounds + res_plot_missing)) + 
    plot_layout(guides = "collect") & 
    theme(legend.position = "bottom"),
  width=12, height=8
)
