library(tidyverse)
load("auxiliary/score-test-simulation-results.Rdata")

results_agg <- 
  results %>%
  gather("rate","reject", starts_with("reject_")) %>%
  mutate(rate = str_sub(rate, -3, -1)) %>%
  group_by(studies, n_factor, mean_effect, sd_effect, p_thresholds, p_RR, type, info, prior_mass, rate) %>%
  summarise(
    replicates = n(),
    reps = sum(reps),
    pct_all_sig = mean(pct_all_sig),
    reject = weighted.mean(reject, w = 1 - pct_NA),
    pct_NA = mean(pct_NA)
  ) %>%
  ungroup() %>%
  select(-n_factor, -p_thresholds, -p_RR, -reps) %>%
  mutate(type_info = paste(type, info, prior_mass, sep = "-"))


# percentage of NA results

results_agg %>%
  filter(rate == "050", type != "subscore") %>%
  ggplot(aes(mean_effect, pct_NA, color = type_info, linetype = type_info)) + 
  geom_point() + geom_line() + 
  expand_limits(y = 0) + 
  facet_grid(studies ~ sd_effect, scales = "free_y", labeller = "label_both") + 
  theme_light()


plot_rejection_rates <- function(dat, scales = "free_y") {
  
  rate <- as.numeric(unique(dat$rate)) / 1000
  
  ggplot(dat, aes(mean_effect, reject, color = type_info, linetype = type_info, shape = type_info)) + 
    geom_point() + geom_line() + 
    scale_color_brewer(type = "qual", palette = 6) + 
    expand_limits(y = 0) + 
    geom_hline(yintercept = rate) + 
    facet_grid(studies ~ sd_effect, scales = scales, labeller = "label_both") + 
    labs(x = "Mean effect size", y = "Rejection rate",
         color = "", linetype = "", shape = "") + 
    theme_light() +
    theme(
      strip.text = element_text(color = "black"),
      legend.position = "bottom"
    )
  
}

# rejection rates at alpha = .025

results_agg %>%
  filter(rate == "025", info == "expected") %>%
  plot_rejection_rates()

# rejection rates at alpha = .05

results_agg %>%
  filter(rate == "050") %>%
  plot_rejection_rates()

results_agg %>%
  filter(rate == "050", info == "expected") %>%
  plot_rejection_rates()

results_agg %>%
  filter(rate == "050", info == "expected", type == "robust") %>%
  plot_rejection_rates()

# rejection rates at alpha = .10

results_agg %>%
  filter(rate == "100") %>%
  plot_rejection_rates()

results_agg %>%
  filter(rate == "100", info == "expected") %>%
  plot_rejection_rates()
