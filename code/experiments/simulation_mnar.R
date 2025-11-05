set.seed(1234)

sims_n <- 200
increment <- 0.05
max_impact_t_o <- 0.4
n_sample <- 100
log_file <- "log_sim_mnar.txt"
parallel <- T

impacts_t_o <- seq(0, max_impact_t_o, increment)

capture.output({
  suppressMessages({
    sfInit(parallel = parallel, cpus = 10, slaveOutfile=log_file)
    sfExportAll()
    sfLibrary(data.table)
    sfLibrary(magrittr)
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
      data_info <- generate_data(
        base_s=base_s,
        impact_t_s=impact_t_s,
        impact_x_s=impact_x_s,
        impact_c_s=impact_c_s,
        base_o=base_o,
        impact_t_o=impact_t_o,
        impact_x_o=impact_x_o,
        impact_c_o=impact_c_o,
        add_impact_s_o=add_impact_s_o, 
        base_m=base_m,
        impact_o_m=impact_o_m,
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
      
      # result_mar <- mar_estimator(dt, col_t='t', cols_w=c("x"))
      # result_mar <- rbind(result_mar, list(
      #   't' = 'delta', p_o_t_mar = result_mar[t == 1, p_o_t_mar] - result_mar[t == 0, p_o_t_mar]
      # ))
      
      result_mnar <- mnar_estimator(
        dt, col_t='t', cols_w=c('x', 's'), apply_bounds=TRUE, 
        estimator=NULL, delta=TRUE
      )
      
      alt_estimator <- function(dt_w, val_t){
        dt_w[(a==1) & (t==val_t), mean(y)]
      }
      result_mnar_alt <- mnar_estimator(
        dt, col_t='t', cols_w=c('x', 's'), apply_bounds=TRUE, 
        estimator=alt_estimator, n_sample=n_sample, delta=TRUE
      )
      result_mnar_alt[t=="delta", p_o_t_mnar := result_mnar[t==1, p_o_t_mnar] - result_mnar[t==0, p_o_t_mnar]]
      result_mnar_alt[t=='logodds', p_o_t_mnar := NA]
      
      setnames(result_mnar_alt, 'p_o_t_mnar', 'p_o_t_mnar_alt')
      result_mnar_alt[, p_o_t_lb:=NULL][, p_o_t_ub:=NULL][, prop_bounds:=NULL]
      
      delta <- p_o_t - p_o_nt
      exp_logodds <- log((p_o_nt/(1-p_o_nt)) / (p_o_t/(1-p_o_t)))
      real <- data.table(t=c(1, 0, 'delta', 'logodds'), p_o_expected = c(p_o_t, p_o_nt, delta, exp_logodds))
      
      sampled <- dt_full[, .(sample_p_o=mean(o)), .(t=as.character(t))]
      sampled <- rbind(sampled, list(
        't' = 'delta', sample_p_o = sampled[t == 1, sample_p_o] - sampled[t == 0, sample_p_o]
      ))
      # sample_p_o_t <- sampled[t == 1, sample_p_o]
      # sample_p_o_nt <- sampled[t == 0, sample_p_o]
      # logodds_ratio <- log((sample_p_o_nt/(1-sample_p_o_nt)) / (sample_p_o_t/(1-sample_p_o_t)))
      dt_sm <- copy(data_info$dt)
      dt_sm[, y :=o]
      logodds_ratio <- calc_logodds_linear_models(dt_sm, col_t='t', cols_w=c('x', 's'))
      sampled <- rbind(sampled, list(
        't' = 'logodds', sample_p_o = logodds_ratio
      ))
      
      observed <- dt_full[, .(observed_p_o=mean(y, na.rm = TRUE)), .(t=as.character(t))]
      observed <- rbind(observed, list(
        't' = 'delta', observed_p_o = observed[t == 1, observed_p_o] - observed[t == 0, observed_p_o]
      ))
      # observed_p_o_t <- observed[t == 1, observed_p_o]
      # observed_p_o_nt <- observed[t == 0, observed_p_o]
      # logodds_ratio <- log((observed_p_o_nt/(1-observed_p_o_nt)) / (observed_p_o_t/(1-observed_p_o_t)))
      logodds_ratio <- NA
      observed <- rbind(observed, list(
        't' = 'logodds', observed_p_o = logodds_ratio
      ))
      
      result <- real %>%
        merge(sampled, by="t") %>%
        merge(observed, by="t") %>%
        merge(dt_full[, .(missing=mean(is.na(y))), .(t=as.character(t))], by="t", all.x=TRUE) %>%
        # merge(result_mar[, t := as.character(t)], by="t") %>%
        merge(result_mnar, by="t") %>%
        merge(result_mnar_alt, by="t")
      result$impact_t_o <- impact_t_o
      result$n <- n
      
      results <- rbind(results, result)
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
    # mar = abs(p_o_t_mar - sample_p_o),
    observed = observed_p_o - sample_p_o,
    impact_t_o,
    t, n
  )] %>%
  melt(id.vars=c("impact_t_o", "t", "n"), variable.name = "method", value.name = "error") %>%
  .[, .(
    mean_error = mean(error, na.rm=TRUE),
    l_ci = quantile(error, probs=0.025, na.rm=TRUE),
    u_ci = quantile(error, probs=0.975, na.rm=TRUE)), 
    .(method, impact_t_o, t, n)
  ] %>%
  ggplot(aes(x = impact_t_o, y = mean_error, linetype = method, fill = method, group = method)) +
  geom_ribbon(aes(ymin = l_ci, ymax = u_ci), alpha = 0.4) +
  geom_line() +
  geom_line(aes(x = impact_t_o, y=u_ci,linetype = method), alpha=0.2) +
  geom_line(aes(x = impact_t_o, y=l_ci,linetype = method), alpha=0.2) +
  facet_grid(n ~ t) +
  ylab("Bias") + 
  xlab("Impact from T to O") + 
  scale_linetype_manual(
    values = c(
      "mnar" = "solid", 
      "mnar_smoothed" = "dashed", 
      # "mar" = "twodash", 
      "observed" = "dotted"
    ),
    labels = c(
      "mnar" = "MNAR estimator", 
      "mnar_smoothed" = "Smoothed MNAR \n estimator", 
      # "mar" = "MAR estimator", 
      "observed" = "Observed"
    )
  ) +
  scale_fill_manual(
    values = c(
      "mnar" = "grey80",
      "mnar_smoothed" = "grey40",
      # "mar" = "grey50",
      "observed" = "grey60"
    ),
    labels = c(
      "mnar" = "MNAR estimator", 
      "mnar_smoothed" = "Smoothed MNAR \n estimator", 
      # "mar" = "MAR estimator", 
      "observed" = "Observed"
    )
  ) +
  scale_alpha(guide = "none") +
  guides(
    fill = guide_legend(
      title = "method",
      override.aes = list(
        # linetype = c("solid", "dashed", "twodash", "dotted"),
        # fill = c("grey80", "grey40", "grey50", "grey60")
        linetype = c("solid", "dashed", "dotted"),
        fill = c("grey80", "grey40", "grey60")
      )
    )
  ) ->
  res_plot_1

results %>%
  .[!(t %chin% c("control", "treatment", 'ATE'))] %>%
  .[, .(
    mnar = p_o_t_mnar - sample_p_o,
    mnar_smoothed = p_o_t_mnar_alt - sample_p_o,
    # mar = abs(p_o_t_mar - sample_p_o),
    observed = observed_p_o - sample_p_o,
    impact_t_o,
    t, n
  )] %>%
  melt(id.vars=c("impact_t_o", "t", "n"), variable.name = "method", value.name = "error") %>%
  .[, .(
    mean_error = mean(error, na.rm=TRUE),
    l_ci = quantile(error, probs=0.025, na.rm=TRUE),
    u_ci = quantile(error, probs=0.975, na.rm=TRUE)), 
    .(method, impact_t_o, t, n)
  ] %>%
  ggplot(aes(x = impact_t_o, y = mean_error, linetype = method, fill = method, group = method)) +
  geom_ribbon(aes(ymin = l_ci, ymax = u_ci), alpha = 0.4) +
  geom_line() +
  geom_line(aes(x = impact_t_o, y=u_ci,linetype = method), alpha=0.2) +
  geom_line(aes(x = impact_t_o, y=l_ci,linetype = method), alpha=0.2) +
  facet_grid(n ~ t, scales="free_y") +
  ylab("Bias") + 
  xlab("Impact from T to O") + 
  scale_linetype_manual(
    values = c(
      "mnar" = "solid", 
      "mnar_smoothed" = "dashed", 
      # "mar" = "twodash", 
      "observed" = "dotted"
    ),
    labels = c(
      "mnar" = "MNAR estimator", 
      "mnar_smoothed" = "Smoothed MNAR \n estimator", 
      # "mar" = "MAR estimator", 
      "observed" = "Observed"
    )
  ) +
  scale_fill_manual(
    values = c(
      "mnar" = "grey80",
      "mnar_smoothed" = "grey40",
      # "mar" = "grey50",
      "observed" = "grey60"
    ),
    labels = c(
      "mnar" = "MNAR estimator", 
      "mnar_smoothed" = "Smoothed MNAR \n estimator", 
      # "mar" = "MAR estimator", 
      "observed" = "Observed"
    )
  ) +
  scale_alpha(guide = "none") +
  guides(
    fill = guide_legend(
      title = "method",
      override.aes = list(
        # linetype = c("solid", "dashed", "twodash", "dotted"),
        # fill = c("grey80", "grey40", "grey50", "grey60")
        linetype = c("solid", "dashed", "dotted"),
        fill = c("grey80", "grey40", "grey60")
      )
    )
  ) ->
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

