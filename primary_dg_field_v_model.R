#hypothesis test similarity in slopes of field and simulated data
#set working directory (not necessary if using as part of a project)
require(tidyverse)
require(rms)
require(nlme)
#simulated dataset
simulated_dg <- read.csv('farber_kingsover_primary_dg.csv',header = T)
#empirical dataset
empirical_dg <- read.csv('island_field_local.csv',header = T)

#best fit c ordinary least squares
best.fit.c<-function(x, y,iterations)
  {c<-c(); pval<-c(); cval<-c(); rsq<-c()
  for(j in iterations)
    {l.bestc<-log10(x + j)
    tmp.mean<-lm(y ~ l.bestc)
    f <- summary(tmp.mean)$fstatistic
    p <- pf(f[1], f[2] , f[3], lower.tail = F)
    r <-summary(tmp.mean)$r.squared
    attributes(p) <- NULL
    pval<-append(pval, p)
    cval<-append(cval, j)
    rsq<-append(rsq, r)}
  candr<-matrix(c(cval, rsq),ncol=2)
  candr.sorted<-candr[order(candr[,2],decreasing=T),]
  candp<-matrix(c(cval, pval),ncol=2)
  candp.sorted<-candp[order(candp[,2],decreasing=F),]
  c.mat<-candr.sorted[1,1]
  l.dist<-log10(x + candr.sorted[1,1])
  print(c.mat)
  return(l.dist)}

#subsets and log transformations
local_sim <- simulated_dg %>% filter(scale == 'local' & m >=0) %>% 
  mutate(log_lesions = log10(lesions)) %>% 
  mutate(log_m = best.fit.c(m,log_lesions,seq(.01,2, by = .01))) #zeros removed

local_emp <- empirical_dg %>% filter(scale == 'local' & m >=0) 
loc_mean<-local_emp %>% group_by(m) %>% summarise(mean_lesions = mean(lesions))
local_emp<-cbind(local_emp,loc_mean[,'mean_lesions']) %>% select(m,mean_lesions)

local_emp_no0<-local_emp[1:42,1:2]
local_emp_no0<-rbind(local_emp_no0,apply(local_emp[43:44,],MARGIN=2,mean))
local_emp_no0<-rbind(local_emp_no0,apply(local_emp[45:46,],MARGIN=2,mean))
local_emp_no0<-rbind(local_emp_no0,apply(local_emp[47:48,],MARGIN=2,mean))
local_emp_no0<-rbind(local_emp_no0,apply(local_emp[49:50,],MARGIN=2,mean))
local_emp_no0<-rbind(local_emp_no0,apply(local_emp[51:52,],MARGIN=2,mean))
local_emp_no0<-rbind(local_emp_no0,local_emp[53:54,])
local_emp_no0<-rbind(local_emp_no0,apply(local_emp[55:56,],MARGIN=2,mean))
local_emp_no0<-rbind(local_emp_no0,apply(local_emp[57:58,],MARGIN=2,mean))
local_emp_no0<-rbind(local_emp_no0,local_emp[59,])
local_emp_no0<-rbind(local_emp_no0,apply(local_emp[60:61,],MARGIN=2,mean))

local_emp_no0 <- local_emp_no0 %>% mutate(log_lesions = log10(mean_lesions)) %>% 
  mutate(log_m = best.fit.c(m,log_lesions,seq(.01,2, by = .01)))


#ordinary linear models  
loc_sim_lm <- lm(local_sim$log_lesions ~ local_sim$log_m)
loc_emp_lm <- lm(local_emp_no0$log_lesions ~ local_emp_no0$log_m)
print(summary(loc_emp_lm))
print(summary(loc_sim_lm))
#lesson: removing zeros will not work.
#plan: subset only the madras 2012 data from local_emp

#dummy vector of 0s for field data and 1s for simulated data
dummy.var <- c(rep(0, length(local_sim[,1])), rep(1, length(local_emp_no0[,1])))
y.new<-c(local_sim$log_lesions,local_emp_no0$log_lesions)
x.new<-c(local_sim$log_m,local_emp_no0$log_m)

#example: gls.mod3 <- gls(y.new~x.new*dummy.var, weights=varIdent(form=~1|dummy.var))

gls.tmp <- gls(y.new~x.new*dummy.var, weights=varIdent(form=~1|dummy.var))
summary(gls.tmp)
#do the same as above for the kingsolver stuff

regional_sim <- simulated_dg %>% filter(scale == 'regional' & m >=0) %>% 
  mutate(log_lesions = log10(lesions)) %>%
  mutate(log_m = best.fit.c(m,log_lesions,seq(.762,150, by = .1)))

regional_emp <- empirical_dg  %>% 
  filter(scale == 'island' & m >=0) %>% mutate(log_lesions = log10(lesions)) %>% 
  mutate(log_m = best.fit.c(m,log_lesions,seq(.762,150, by = .1)))

reg_sim_lm <- lm(regional_sim$log_lesions ~ regional_sim$log_m)
reg_emp_lm <- lm(regional_emp$log_lesions ~ regional_emp$log_m)
print(summary(reg_emp_lm))
print(summary(reg_sim_lm))

#lesson: removing zeros will not work.
#plan: subset only the madras 2012 data from local_emp
  
#dummy vector of 0s for field data and 1s for simulated data
dummy.var <- c(rep(0, length(local_sim[,1])), rep(1, length(local_emp_no0[,1])))
y.new<-c(local_sim$log_lesions,local_emp_no0$log_lesions)
x.new<-c(local_sim$log_m,local_emp_no0$log_m)
           