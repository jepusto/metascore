library(metafor)
source("paper-3PSM-and-excess-significance-test/simulation-functions.R")
source("auxiliary/score-CR.R")

dat <- r_SMD(studies = 80, mean_effect = 0.2, sd_effect = 0.1,
             n_sim = n_beta(20,100,1,1.5), 
             p_thresholds = .025, p_RR = 1)

dat$studyid <- sample(LETTERS, size = nrow(dat), replace = TRUE)
table(dat$studyid)

# fixed effect meta-analysis
rma(yi = g, vi = Va, data = dat, method = "FE")

# weighted least squares meta-analysis (FE point estimate, REML tau-sq estimate)
meta_WLS <- rma(yi = g, vi = Va, weights = 1 / Va, data = dat, method = "REML")
meta_WLS
score_CR(mod = meta_WLS, cluster = dat$studyid, alpha = .025)

# random effects meta-analysis
meta_RML <- rma(yi = g, vi = Va, data = dat, method = "REML")
meta_RML
score_CR(mod = meta_RML, cluster = dat$studyid, alpha = .025)
