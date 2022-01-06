###As described in Tsypes et al., in press, Journal of Abnormal Psychology,
#"Delay discounting in suicidal behavior: Myopic preference or inconsistent valuation?"
#Code for the simulation: "We simulated the behavior of 800 pseudosubjects on the delay discounting task with the true
#hyperbolic discount rates of  3.95 (Group 1, n = 100) and of
#5.64 (Group 2, n = 100) with and without the noise injected into
#the delayed value (noise levels: 0, .1, .33, .67). The two discount
#rates were based on actual estimated means of the two participant
#groups from our prior study (that is, low lethality suicide attempters
#versus controls; (Dombrovski et al., 2011). To manipulate valuation
#consistency, we injected noise into delayed, but not
#immediate option, values because only the computation of delayed
#option values is cognitively demanding. We then recovered
#between-group differences in true hyperbolic discount rates and
#levels of noise via the MLM and using two-stage estimation
#method (Kirby et al., 1999) and compared the test statistics generated
#through this recovery. For comparison purposes, mean coefficients
#from these models were converted into z statistics through
#division by corresponding standard deviations."

##k_sub: one's individual discount rate in two-stage approaches to estimating individual discount rates and choice consistencies (p.4 JoAP for more details)
##log_k_sub: log-transformed individual discount rate 
##consistency: proportion of choices consistent with the most likely k_sub, 
#in two-stage approaches to estimating individual discount rates and choice consistencies (p.4 JoaP for more details)

library(bbmle)
library(lme4)
library(tidyr)
library(corrplot)
library(dplyr)
library(tidyverse)
library(psych)
library(corrplot)
library(lme4)
library(ggpubr)
library(car)
library(readxl)
library(compareGroups)
library(haven)
library(emmeans)

num_subjects=100 #number of subjects per group
# set up choices from MCQ
MCQ_options=matrix(data=c(
  54,55,117,
  55,75,61,
  19,25,53,
  84.5,85,155,
  14,25,19,
  47,50,160,
  15,35,13,
  25,60,14,
  78,80,162,
  40,55,62,
  34.75,35,186,
  67,75,119,
  34,35,186,
  27,50,21,
  69,85,91,
  49,60,89,
  80,85,157,
  24,35,29,
  33,80,14,
  28,30,179,
  34,50,30,
  25,30,80,
  41,75,20,
  54,60,111,
  54,80,30,
  22,25,136,
  59.75,60,109,
  34.5,35,186,
  84,85,150,
  59.5,60,108),nrow=30,byrow=TRUE) 
MCQ_options=as.data.frame(MCQ_options)
names(MCQ_options)=c('imm_reward','delay_reward','delay')

#compute discount rate for indifference point
MCQ_options$k_values=(MCQ_options$delay_reward/MCQ_options$imm_reward-1)/MCQ_options$delay
MCQ_options$logk=log(MCQ_options$k_values)
MCQ_options$logk_sc=scale(MCQ_options$logk)

###assigning values--means of controls and LL from 2011 paper
discounted_values_grp1=MCQ_options$delay_reward*
  (1/(1+exp(-3.95)*MCQ_options$delay))

discounted_values_grp2=MCQ_options$delay_reward*
  (1/(1+exp(-5.64)*MCQ_options$delay))

noises <- c(0, .33, .67, .1)
####Groups 1 and 2, noise 0
MCQ_choices1=array(data=NA,dim=c(num_subjects*2,30))
for (s in 1:num_subjects) {
  for (q in 1:30) {
    #group 1 randomly picks a value, mean is imm reward and SD is based on noise times that reward
    noisy_imm_reward1_1=rnorm(1,MCQ_options$imm_reward[q],0*MCQ_options$imm_reward[q])
    noisy_discounted_value1_1=rnorm(1,discounted_values_grp1[q],noises[1]*discounted_values_grp1[q]) 
    
    MCQ_choices1[s,q]=ifelse(noisy_imm_reward1_1>noisy_discounted_value1_1,0,1) #0 for imm; 1 for delayed (if first argument is true, it assigns the value of the second; if immreward is higher, then choice 1; if discounte d value is higher it's 0)
    
    #group 2 
    noisy_imm_reward2_1=rnorm(1,MCQ_options$imm_reward[q],0*MCQ_options$imm_reward[q])
    noisy_discounted_value2_1=rnorm(1,discounted_values_grp2[q],noises[1]*discounted_values_grp2[q])
    
    MCQ_choices1[(num_subjects+s),q]=ifelse(noisy_imm_reward2_1>noisy_discounted_value2_1,0,1) 
  }
}
MCQ_choices1=data.frame(MCQ_choices1)
names(MCQ_choices1)=paste0('Q',seq(1,30,by=1))
MCQ_choices_long1=pivot_longer(MCQ_choices1,cols=names(MCQ_choices1),names_to='Question')
MCQ_choices_long1$Question=as.numeric(sub('^Q','',MCQ_choices_long1$Question))
MCQ_choices_long1$k=MCQ_options$k_values[MCQ_choices_long1$Question]
MCQ_choices_long1$logk_sc=MCQ_options$logk_sc[MCQ_choices_long1$Question]
MCQ_choices_long1$group=as.factor(c(rep(1,num_subjects*30),rep(2,num_subjects*30)))
MCQ_choices_long1$ID=c(rep(1:(2*num_subjects),each=30))

######Kirby's for noise level 0
MCQ_choices_long1_renamed<-MCQ_choices_long1 %>% 
  dplyr::rename(choice = value)
df<-MCQ_choices_long1_renamed

df <- df %>% arrange(ID, k)
ks <- unique(df$k)
ids <- unique(df$ID)
df$consistency <- NA
df$k_sub <- NA
df$max_consistency <- NA
for (id in ids) {
  for (k in ks) {
    df$consistency[df$ID==id & df$k==k] = (sum(df$ID==id & df$k<k & df$choice==0, na.rm = T) + sum(df$ID==id & df$k>k & df$choice==1, na.rm = T))/(sum(!is.na(df$choice[df$ID==id]))-1)
  }
  best <- df %>% filter(ID==id & consistency == max(consistency[ID==id])) %>% select(k, consistency)
  df$k_sub[df$ID==id] <- geometric.mean(best$k)
  df$max_consistency[df$ID==id] <- max(best$consistency)
}
df$log_k_sub = log(df$k_sub)

sub_df1 <- df %>% select(ID, k_sub, log_k_sub, max_consistency, group) %>% unique()


####Groups 1 and 2, noise 0.33
MCQ_choices2=array(data=NA,dim=c(num_subjects*2,30))
for (s in 1:num_subjects) {
  for (q in 1:30) {
    #group 1 randomly picks a value, mean is imm reward and SD is based on noise times that reward
    noisy_imm_reward1_2=rnorm(1,MCQ_options$imm_reward[q],0*MCQ_options$imm_reward[q])
    noisy_discounted_value1_2=rnorm(1,discounted_values_grp1[q],noises[2]*discounted_values_grp1[q]) 
    
    MCQ_choices2[s,q]=ifelse(noisy_imm_reward1_2>noisy_discounted_value1_2,0,1) 
    #group 2 
    noisy_imm_reward2_2=rnorm(1,MCQ_options$imm_reward[q],0*MCQ_options$imm_reward[q])
    noisy_discounted_value2_2=rnorm(1,discounted_values_grp2[q],noises[2]*discounted_values_grp2[q])
    
    MCQ_choices2[(num_subjects+s),q]=ifelse(noisy_imm_reward2_2>noisy_discounted_value2_2,0,1) 
  }
}
MCQ_choices2=data.frame(MCQ_choices2)
names(MCQ_choices2)=paste0('Q',seq(1,30,by=1))
MCQ_choices_long2=pivot_longer(MCQ_choices2,cols=names(MCQ_choices2),names_to='Question')
MCQ_choices_long2$Question=as.numeric(sub('^Q','',MCQ_choices_long2$Question))
MCQ_choices_long2$k=MCQ_options$k_values[MCQ_choices_long2$Question]
MCQ_choices_long2$logk_sc=MCQ_options$logk_sc[MCQ_choices_long2$Question]
MCQ_choices_long2$group=as.factor(c(rep(1,num_subjects*30),rep(2,num_subjects*30)))
MCQ_choices_long2$ID=c(rep(1:(2*num_subjects),each=30))

######Kirby's for noise level noises[2]
MCQ_choices_long2_renamed<-MCQ_choices_long2 %>% 
  dplyr::rename(choice = value)
df<-MCQ_choices_long2_renamed

df <- df %>% arrange(ID, k)
ks <- unique(df$k)
ids <- unique(df$ID)
df$consistency <- NA
df$k_sub <- NA
df$max_consistency <- NA
for (id in ids) {
  for (k in ks) {
    df$consistency[df$ID==id & df$k==k] = (sum(df$ID==id & df$k<k & df$choice==0, na.rm = T) + sum(df$ID==id & df$k>k & df$choice==1, na.rm = T))/(sum(!is.na(df$choice[df$ID==id]))-1)
  }
  best <- df %>% filter(ID==id & consistency == max(consistency[ID==id])) %>% select(k, consistency)
  df$k_sub[df$ID==id] <- geometric.mean(best$k)
  df$max_consistency[df$ID==id] <- max(best$consistency)
}
df$log_k_sub = log(df$k_sub)
sub_df2 <- df %>% select(ID, k_sub, log_k_sub, max_consistency, group) %>% unique()

####Groups 1 and 2, noise noises[3]
MCQ_choices3=array(data=NA,dim=c(num_subjects*2,30))
for (s in 1:num_subjects) {
  for (q in 1:30) {
    #group 1 randomly picks a value, mean is imm reward and SD is based on noise times that reward
    noisy_imm_reward1_3=rnorm(1,MCQ_options$imm_reward[q],0*MCQ_options$imm_reward[q])
    noisy_discounted_value1_3=rnorm(1,discounted_values_grp1[q],noises[3]*discounted_values_grp1[q]) 
    
    MCQ_choices3[s,q]=ifelse(noisy_imm_reward1_3>noisy_discounted_value1_3,0,1) 
    #group 2 
    noisy_imm_reward2_3=rnorm(1,MCQ_options$imm_reward[q],0*MCQ_options$imm_reward[q])
    noisy_discounted_value2_3=rnorm(1,discounted_values_grp2[q],noises[3]*discounted_values_grp2[q])
    
    MCQ_choices3[(num_subjects+s),q]=ifelse(noisy_imm_reward2_3>noisy_discounted_value2_3,0,1) 
  }
}
MCQ_choices3=data.frame(MCQ_choices3)
names(MCQ_choices3)=paste0('Q',seq(1,30,by=1))
MCQ_choices_long3=pivot_longer(MCQ_choices3,cols=names(MCQ_choices3),names_to='Question')
MCQ_choices_long3$Question=as.numeric(sub('^Q','',MCQ_choices_long3$Question))
MCQ_choices_long3$k=MCQ_options$k_values[MCQ_choices_long3$Question]
MCQ_choices_long3$logk_sc=MCQ_options$logk_sc[MCQ_choices_long3$Question]
MCQ_choices_long3$group=as.factor(c(rep(1,num_subjects*30),rep(2,num_subjects*30)))
MCQ_choices_long3$ID=c(rep(1:(2*num_subjects),each=30))

######Kirby's for noise noises[3]
MCQ_choices_long3_renamed<-MCQ_choices_long3 %>% 
  dplyr::rename(choice = value)
df<-MCQ_choices_long3_renamed
df <- df %>% arrange(ID, k)
ks <- unique(df$k)
ids <- unique(df$ID)
df$consistency <- NA
df$k_sub <- NA
df$max_consistency <- NA
for (id in ids) {
  for (k in ks) {
    df$consistency[df$ID==id & df$k==k] = (sum(df$ID==id & df$k<k & df$choice==0, na.rm = T) + sum(df$ID==id & df$k>k & df$choice==1, na.rm = T))/(sum(!is.na(df$choice[df$ID==id]))-1)
  }
  best <- df %>% filter(ID==id & consistency == max(consistency[ID==id])) %>% select(k, consistency)
  df$k_sub[df$ID==id] <- geometric.mean(best$k)
  df$max_consistency[df$ID==id] <- max(best$consistency)
}
df$log_k_sub = log(df$k_sub)
sub_df3 <- df %>% select(ID, k_sub, log_k_sub, max_consistency, group) %>% unique()


####Groups 1 and 2, noise noises 4
MCQ_choices4=array(data=NA,dim=c(num_subjects*2,30))
for (s in 1:num_subjects) {
  for (q in 1:30) {
    #group 1 randomly picks a value, mean is imm reward and SD is based on noise times that reward
    noisy_imm_reward1_4=rnorm(1,MCQ_options$imm_reward[q],0*MCQ_options$imm_reward[q])
    noisy_discounted_value1_4=rnorm(1,discounted_values_grp1[q],noises[4]*discounted_values_grp1[q]) 
    
    MCQ_choices4[s,q]=ifelse(noisy_imm_reward1_4>noisy_discounted_value1_4,0,1) 
    #group 2 
    noisy_imm_reward2_4=rnorm(1,MCQ_options$imm_reward[q],0*MCQ_options$imm_reward[q])
    noisy_discounted_value2_4=rnorm(1,discounted_values_grp2[q],noises[4]*discounted_values_grp2[q])
    
    MCQ_choices4[(num_subjects+s),q]=ifelse(noisy_imm_reward2_4>noisy_discounted_value2_4,0,1) 
  }
}
MCQ_choices4=data.frame(MCQ_choices4)
names(MCQ_choices4)=paste0('Q',seq(1,30,by=1))
MCQ_choices_long4=pivot_longer(MCQ_choices4,cols=names(MCQ_choices4),names_to='Question')
MCQ_choices_long4$Question=as.numeric(sub('^Q','',MCQ_choices_long4$Question))
MCQ_choices_long4$k=MCQ_options$k_values[MCQ_choices_long4$Question]
MCQ_choices_long4$logk_sc=MCQ_options$logk_sc[MCQ_choices_long4$Question]
MCQ_choices_long4$group=as.factor(c(rep(1,num_subjects*30),rep(2,num_subjects*30)))
MCQ_choices_long4$ID=c(rep(1:(2*num_subjects),each=30))

######Kirby's for noise level 0.33
MCQ_choices_long4_renamed<-MCQ_choices_long4 %>% 
  dplyr::rename(choice = value)
df<-MCQ_choices_long4_renamed
df <- df %>% arrange(ID, k)
ks <- unique(df$k)
ids <- unique(df$ID)
df$consistency <- NA
df$k_sub <- NA
df$max_consistency <- NA
for (id in ids) {
  for (k in ks) {
    df$consistency[df$ID==id & df$k==k] = (sum(df$ID==id & df$k<k & df$choice==0, na.rm = T) + sum(df$ID==id & df$k>k & df$choice==1, na.rm = T))/(sum(!is.na(df$choice[df$ID==id]))-1)
  }
  best <- df %>% filter(ID==id & consistency == max(consistency[ID==id])) %>% select(k, consistency)
  df$k_sub[df$ID==id] <- geometric.mean(best$k)
  df$max_consistency[df$ID==id] <- max(best$consistency)
}
df$log_k_sub = log(df$k_sub)
sub_df4 <- df %>% select(ID, k_sub, log_k_sub, max_consistency, group) %>% unique()


MCQ_choices_long1_renamed$noise <- noises[1]
MCQ_choices_long2_renamed$noise <- noises[2]
MCQ_choices_long2_renamed$ID <- MCQ_choices_long2_renamed$ID + 200
MCQ_choices_long3_renamed$noise <- noises[3]
MCQ_choices_long3_renamed$ID <- MCQ_choices_long3_renamed$ID + 400
MCQ_choices_long4_renamed$noise <- noises[4]
MCQ_choices_long4_renamed$ID <- MCQ_choices_long4_renamed$ID + 600

ldf <- rbind(MCQ_choices_long1_renamed,
             MCQ_choices_long2_renamed,
             MCQ_choices_long3_renamed,
             MCQ_choices_long4_renamed)
# ldf <- MCQ_choices_long_combined
ldf <- ldf %>% mutate(log_k_true = case_when(
  group==1 ~'-3.95',
  group==2 ~ '-5.64'
), noise = as.factor(noise))


# df <- SubjectLevelCombinedWithRE
sub_df1$noise <- noises[1]
sub_df2$noise <- noises[2]
sub_df3$noise <- noises[3]
sub_df4$noise <- noises[4]
df <- rbind(sub_df1, sub_df2, sub_df3, sub_df4)
setwd('~/code/MCQ/')

# s <- load('Workspace_combined.RData')
df <- df %>% mutate(log_k_true = case_when(
  group==1 ~'-3.95',
  group==2 ~ '-5.64'
), noise = as.factor(noise))

###MLM to calculate valuation and preference parameters for simulated data

mb1 <- stan_glmer(choice ~ logk_sc * noise + logk_sc * log_k_true +
                    (1|ID), family = binomial, ldf, chains=4, iter = 4000)
summary(mb1)

post1 <- describe_posterior(mb1, centrality = 'median', test = c('p_direction', 'p_significance'))
post2 <- insight::get_parameters(mb1)
post2$k <- -((post2$`log_k_true-5.64`)/(post2$logk_sc + post2$`logk_sc:log_k_true-5.64`))
post2$k01 <- -((post2$noise0.1)/(post2$logk_sc + post2$`logk_sc:noise0.1`))
post2$k033 <- -((post2$noise0.33)/(post2$logk_sc + post2$`logk_sc:noise0.33`))
post2$k067 <- -((post2$noise0.67)/(post2$logk_sc + post2$`logk_sc:noise0.67`))

#test
describe_posterior(post2)
describe(post2)
rstan::summary(post2)

###plot the data to examine posteriors
ggplot(post2, aes(x = k)) +
  geom_density(fill = "orange") + xlim(-2,2)
ggplot(post2, aes(x = k01)) +
  geom_density(fill = "orange") + xlim(-2,2)
ggplot(post2, aes(x = k033)) +
  geom_density(fill = "orange") + xlim(-2,2)
ggplot(post2, aes(x = k067)) +
  geom_density(fill = "orange") + xlim(-2,2)

# KIRBY RECOVERING SUBJECT'S K
mk1 <- stan_lm(log_k_sub ~ log_k_true * noise, data = df, prior=NULL, chains=4, iter=4000)
summary(mk1)
describe_posterior(mk1)
post_mk1 <- describe_posterior(mk1, centrality = 'median', test = c('p_direction', 'p_significance'))
post2_mk1 <- insight::get_parameters(mk1)
describe(post2_mk1)

# KIRBY RECOVERING CONSISTENCY
mk2 <- stan_lm(max_consistency ~ log_k_true * noise, data = df, prior=NULL, chains=4, iter=4000)
summary(mk2)
describe_posterior(mk2)
post_mk2 <- describe_posterior(mk2, centrality = 'median', test = c('p_direction', 'p_significance'))
post2_mk2 <- insight::get_parameters(mk2)
describe(post2_mk2)

###Figure 1 in in the JoAP:
###Flipped the sign of log_k_true below (When the true k is lower, 
#Kirby estimates a lower recovered k whereas the MLM predicts a delayed choice, 
#which depends on how delayed/immediate are coded.)
A<- c(-20.75, -2.5, -11.5, -13.33)
B<- c("True discount rate", "Consistency, \n 10% noise", "Consistency, \n 33% noise ", "Consistency, \n 67% noise")
d <- as_tibble(A, names = "statistic")
d$method <- 'MLM'
d$predictor <- B
###converting t to z for lm predicting log_k_sub z=t*sqrt(N)/sqrt(df[error]).
C<-c(-18.7, -1,-4, -9)
# the vector of Kirby/lm statistic should have (1) effect of true k on recovered k and (2-4) effects of true noise levels on consistencies
D<-c("True discount rate", "Consistency, \n 10% noise", "Consistency, \n 33% noise ", "Consistency, \n 67% noise")
# barplot(C, names.arg=D, ylab="z score")
c <- as_tibble(C, names = "statistic")
c$method = 'Two-stage estimation'
c$predictor <- D
k <- rbind(d,c)
k1<-k%>% 
  dplyr::rename(
    Method = method)
ggplot(k1, aes(predictor, value, color=Method)) + geom_point(stat="identity", size = 3.5) + xlab("Predictor") + ylab("z statistic")


