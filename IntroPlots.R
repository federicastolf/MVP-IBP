rm(list=ls())

library(latex2exp)
library(tidyverse)
load("data/fungi_binary.Rdata")

####---------------------###### histogram ####--------------------######

m1 = apply(fungi, 2, mean)
fungi_prevalence = data.frame(Species = seq(1, NCOL(fungi)), Prevalence = m1)

sz_t = 22
pl1 = ggplot(fungi_prevalence) + 
  geom_histogram(aes(x = Prevalence), binwidth = 0.02) + 
  xlab("Mean prevalence") + 
  ylab("Count") + 
  theme_light()+
  theme(axis.title.x = element_text(size = 22),axis.text.x = element_text(size = 22),
        axis.title.y = element_text(size = 22), axis.text.y = element_text(size = 22))
pl1

####---------------------###### common species ####--------------------######

# empirical common species
Demp = seq(0,1,by=0.01)
common_emp = rep(0, length(Demp))
for(i in 1:length(Demp)){
  common_emp[i] = sum(m1>Demp[i])
}

mean(rowSums(fungi))
alpha = 70
common_PIBP = alpha*exp(-1/2-qnorm(Demp)) # common species for MVP-IBP
common_IBP = -alpha*log(Demp) # common species for IBP

dataC = cbind.data.frame(Demp, common_emp, common_PIBP, common_IBP)
dataC_long <- dataC %>% gather(key = "type", value = "value", -Demp)
dataC_long$type = as.factor(dataC_long$type)

fcn = ggplot(dataC_long, aes(x = Demp, y = value)) + 
  geom_line(aes(color = type, linetype=type), size=1.8) +
  scale_linetype_manual(values = c("longdash", "solid", "solid"),
                        labels = c("empirical","IBP", "MVP-IBP")) +
  scale_color_manual(values=c("steelblue", "black", "darkred"),
                     labels = c("empirical","IBP", "MVP-IBP")) +
  
  ylab("Common features") +
  xlab(TeX("$\\Delta$")) +
  theme_light() +
  theme(legend.position = "top", axis.title.y = element_text(size=21),
        axis.title.x = element_text(size=24), legend.text = element_text(size=20),
        axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),
        legend.title=element_blank(),  legend.box.spacing = unit(0.2,"line"),
        strip.text = element_text(size = 10, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"))
fcn

