
library(Rcpp) 
library(RcppArmadillo)
library(tidyverse)
require(phyloseq)
library(psadd)
library(dplyr)
library(latex2exp) 

rm(list=ls()) 
Rcpp::sourceCpp("TRACEh.cpp")
source("plot_functions.R")

#-------# load data #----#
load("data/fungi_binary.Rdata")
load("Data/X_fungi.Rdata")
load("Data/tree_fungi.Rdata")


#----------------# set parameters #-----------------#
max_it = 100
epsilon = 0.0001
m = 15
nmcmc = 200
burnin = 50
a_alpha = 17.5
b_alpha = 0.25
eps_MH = 0.2
prior_var = 10
eta0_rho = 0 
nu0_rho = 0.01 
gamma0_rho = 0.01 # a_omega
lambda0_sq_rho = 0.01 # b_omega
truncP = 200

#------------------------# fit TRACE for fungi data #--------------------------#
fit_TRACE = TRACEh(fungi, X_fungi,eta0_rho, nu0_rho, gamma0_rho, lambda0_sq_rho, 
                      a_alpha, b_alpha,prior_var,max_it, epsilon, m, nmcmc, burnin, 
                      eps_MH, truncP) 


#####------------------------#### Plots #####------------------------####

#--------# Krona wheels plots #-----#

n = nrow(X_fungi)
# constuct a dataframe with site as factor
df1 = cbind.data.frame(id=c(1:n), X_fungi[,2:4])
colnames(df1)[2:4]=c("1","3","4")
covariate_fungi = df1 %>% 
  pivot_longer(-id) %>%
  group_by(id) %>%
  slice(which.max(as.integer(factor(name))*value))%>%
  mutate(name = if_else(value == 0, '2',name), value= NULL)
colnames(covariate_fungi)[2] = "site"
covariate_fungi$site = as.factor(covariate_fungi$site)

# matrix with taxonomic tree
tax_matAB = tree_fungi %>% filter(Phylum %in% c("Ascomycota","Basidiomycota"))
tax_matAB = as.matrix(tax_matAB)

phylo_fungiAB = build_phylo(covariate_fungi, fit_TRACE, tax_matAB, X_fungi)

# interactive krona plots
# type "no_cov" instead of "site" for the taxonomic composition in all samples
plot_krona(phylo_fungiAB, "fungiAB", "site", trim=T)

#--------# sample-specific species richness #-----#

sprich = sprichnness_site(fit_TRACE, X_fungi)
plot_spr = cbind.data.frame(covariate_fungi$site, X_fungi[,5], sprich)
colnames(plot_spr) = c("site", "Week","sp")
plotsr = plot_spr %>% distinct()


figsp = ggplot(plot_spr, aes(x=as.factor(site), y=sp, colour = Week))+ 
  geom_point(size=5.5) +
  scale_color_gradient(low="yellow", high="darkred") +
  xlab("Site") + ylab(TeX("$n_i$")) +
  theme_light() + theme( axis.title=element_text(size=12),
                         axis.text.x = element_text(size=12),
                         axis.text.y = element_text(size=11))
figsp 

# ggsave(filename = "spSite.png", plot=figsp  ,  width = 8, height = 4)

