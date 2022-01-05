###Calculates valuation and preference parameters for participant data (Sample 1 here as an example, 
#as described in detail in the Tsypes et al., in press, Journal of Abnormal Psychology publication)
##For more details, see Tsypes, A., Szanto, K., Bridge, J. A., Brown, V. M., Keilp, J. G., & Dombrovski, A. Y. (in press). 
#Delayed reward valuation in high-lethality suicidal behavior: Myopic Preference or Inconsistent Valuation? 
#Journal of Abnormal Psychology.
# rstan ----
library(rstan)
library(rstanarm)
library(bayestestR)
library(effectsize)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

###Models with participant data
####### Read in tall-format delay discounting data, build multi-level models
setwd("~/OneDrive/discounting/data")
load('non_afsp_subs_long.Rda')

###change reference group for GLMs to high lethality
non_afsp_subs_long$lethgrp_ref_hl <- relevel(non_afsp_subs_long$groupLeth, ref = 'HL')

####Preference and valuation calculations with participant data, sample 1 in the JoAP paper
# example code for calculation of MCQ item-level k/logk_sc (i.e., k item indifference) that is used in the MLM below 
df <- as_tibble(df)  %>%
  dplyr::rename(immMag = Immediate.magnitude,
                delayMag = Delayed.magnitude,
                delay = Lenth.of.delay,
                choice = Response,
                item = Item) %>%
  mutate(magRatio = delayMag/immMag,
         logMagRatio = log(magRatio),
         k = (magRatio-1)/delay,
         logk = log(k),
         logk_sc = scale(logk))

####general model without covariates 
mb2<-stan_glmer(choice ~ logk_sc * lethgrp_ref_hl + (1|ID), family = binomial, non_afsp_subs_long, chains = 4, iter = 4000)
summary(mb2)
post1a <- describe_posterior(mb2, centrality = 'median', test = c('p_direction', 'p_significance'))
post2a <- insight::get_parameters(mb2)
###calculate preference parameter (i.e., discount rate) for each clinical group: k = -intercept/beta for each group reflecting that
#group’s discount rate
post2a$kHC <- -(post2a$lethgrp_ref_hlHC/(post2a$logk_sc + post2a$`logk_sc:lethgrp_ref_hlHC`))
post2a$kDEP <- -(post2a$lethgrp_ref_hlDEP/(post2a$logk_sc + post2a$`logk_sc:lethgrp_ref_hlDEP`))
post2a$kIDE <- -(post2a$lethgrp_ref_hlIDE/(post2a$logk_sc + post2a$`logk_sc:lethgrp_ref_hlIDE`))
post2a$kLL <-  -(post2a$lethgrp_ref_hlLL/(post2a$logk_sc + post2a$`logk_sc:lethgrp_ref_hlLL`))

#test
describe_posterior(post2a)

# plot consistency estimates
post2alc <- pivot_longer(post2a,values_to = "coefficient", cols = c("logk_sc:lethgrp_ref_hlHC", "logk_sc:lethgrp_ref_hlDEP",                                                                    "logk_sc:lethgrp_ref_hlIDE", "logk_sc:lethgrp_ref_hlLL"), names_to = "contrast")
ggplot(post2alc) + geom_density(aes(x = coefficient, color = contrast, ..scaled..)) + xlim(-3,3)
str(post2alc)





