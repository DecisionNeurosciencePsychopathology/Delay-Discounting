###Calculates valuation and preference parameters for participant data (Sample 1 here as an example)
##For more details, see Tsypes, A., Szanto, K., Bridge, J. A., Brown, V. M., Keilp, J. G., & Dombrovski, A. Y. (in press). 
#Delayed reward valuation in high-lethality suicidal behavior: Myopic Preference or Inconsistent Valuation? Journal of Abnormal Psychology.

###k item indifference: Discount rates that correspond to the indifference point for a given item consisting
#of a smaller immediate and a larger delayed reward options (Method, Equation 1). logk_sc is k item indifference
#scaled and log-transformed

##valuation consistency (beta coefficients in mb2 below): A degree to which choices are sensitive to the relative expected values of the
#available options. This parameter is analogous to the inverse temperature of the softmax choice function (Method, Equation 4).

###k subject: Subject-level discount rate that captures individual preferences for smaller, immediate, versus larger, delayed rewards.
#k subject in mb2 below = -intercept/beta for each group

library(rstan)
library(rstanarm)
library(bayestestR)
library(effectsize)
library(ggplot2)
library(wesanderson)
library(dplyr)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

###Models with participant data
####### Read in tall-format delay discounting data, build multi-level models
setwd("~/OneDrive/discounting/data")
load('non_afsp_subs_long.Rda')
###change reference group for GLMs to high lethality
non_afsp_subs_long$lethgrp_ref_hl <- relevel(non_afsp_subs_long$groupLeth, ref = 'HL')

# example code for calculation and log transformation and scaling of k item indifference used in the mb2 model below 
# choice (now or later), immMag, delayMag, delay, and item (number) are MCQ item characteristics
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

####Preference and valuation calculations with participant data, sample 1 in the JoAP paper (Table 5)
####general model without covariates 
mb2<-stan_glmer(choice ~ logk_sc * lethgrp_ref_hl + (1|ID), family = binomial, non_afsp_subs_long, chains = 4, iter = 4000)
summary(mb2)
post1a <- describe_posterior(mb2, centrality = 'median', test = c('p_direction', 'p_significance'))
post2a <- insight::get_parameters(mb2)
###calculate preference parameter/discount rate (i.e., k subject) for each clinical group: k subject = -intercept/beta for each group reflecting that
#group’s discount rate
post2a$kHC <- -(post2a$lethgrp_ref_hlHC/(post2a$logk_sc + post2a$`logk_sc:lethgrp_ref_hlHC`))
post2a$kDEP <- -(post2a$lethgrp_ref_hlDEP/(post2a$logk_sc + post2a$`logk_sc:lethgrp_ref_hlDEP`))
post2a$kIDE <- -(post2a$lethgrp_ref_hlIDE/(post2a$logk_sc + post2a$`logk_sc:lethgrp_ref_hlIDE`))
post2a$kLL <-  -(post2a$lethgrp_ref_hlLL/(post2a$logk_sc + post2a$`logk_sc:lethgrp_ref_hlLL`))

#test
describe_posterior(post2a)

##Sample 1
# plot consistency estimates, high lethality attempters as a reference group (Figure 2, plot F in the JoAP)
pal = wes_palette("Zissou1", 4, type = "continuous")
post2alc <- pivot_longer(post2a,values_to = "coefficient", cols = c("logk_sc:lethgrp_ref_hlHC", "logk_sc:lethgrp_ref_hlDEP",
                                                                    "logk_sc:lethgrp_ref_hlIDE", "logk_sc:lethgrp_ref_hlLL"), names_to = "contrast")
post2alc <- dplyr::rename (post2alc, Group = "contrast")
post2alc_1<-post2alc %>% mutate(Group = case_when(
  Group == 'logk_sc:lethgrp_ref_hlHC' ~ 'Controls',
  Group == 'logk_sc:lethgrp_ref_hlDEP' ~ 'MDD',
  Group == 'logk_sc:lethgrp_ref_hlIDE' ~ 'SI + MDD',
  Group == 'logk_sc:lethgrp_ref_hlLL' ~ "LL SA + MDD"
))

post2alc_1$Group <- factor(post2alc_1$Group, levels = c("Controls", "MDD", "SI + MDD", "LL SA + MDD"))

pal = wes_palette("Zissou1", type = "discrete")
f<-ggplot(post2alc_1) + geom_density(aes(x = coefficient, color = Group)) + xlim(-3,3) + labs (y='Posterior Density', x='valuation consistency vs. HL SA + MDD') +
  scale_color_manual(values = pal) + 
  theme(panel.grid.major = element_line(colour = "grey45"), 
        panel.grid.minor = element_line(colour = "grey45"), 
        panel.background = element_rect(fill = 'grey40'),aspect.ratio=1)  + guides(color=guide_legend(override.aes = list(size=3))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray")

##plot
f








