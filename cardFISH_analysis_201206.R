setwd("~/Dropbox/Projects/078_cardFISH/analysis/")

#load data
d <- read.csv2("cardFISH_counts_timepoint_experiment_JW1.csv")

# rename columns
colnames(d) <- c("Dish", "Sample", "Height", "Timepoint", "Archaea", "Bacteria", "NOB")
head(d)

# load package
require("brms")
require("lme4")

# fit model using lme4 (standard mixed effect model)
mod_lme4 <- lmer(Archaea ~ as.factor(Timepoint)*as.factor(Height) + (1|Dish/Sample), data = d)
summary(mod_lme4)

# in this case, we are testing whether (1) the time points differ, (2) the heights differ,
# and (3) the time differences change depending on height. The random effect part
# tells the model to fit a separate intercept for each dish and each sample, accounting
# for the fact that samples are grouped by dish. This should account for any
# effects of pseudoreplication within dishes.

# Note that we have wrapped time and height into as.factor() commands, since we want
# these to be treated as categorical variables.


# lme4 doesn't give p-values for fixed effects, so we need to compare different
# subsets of models to test hypotheses. E.g.:

mod_lme4_a1 <- lmer(Archaea ~ as.factor(Timepoint) + as.factor(Height) + (1|Dish/Sample), data = d)
mod_lme4_a2 <- lmer(Archaea ~ as.factor(Timepoint) + (1|Dish/Sample), data = d)
mod_lme4_a3 <- lmer(Archaea ~ 1 + (1|Dish/Sample), data = d)

# BUT note that none of these converge -- we really don't have enough data for a 
# standard freqentist framework. So, instead we can try a Bayesian approach with brms.

# mod_brm <- brm(Archaea ~ as.factor(Timepoint)*as.factor(Height) + (1|Dish/Sample), data = d)

# Note that syntax is the same as for lme4. BUT it will take much longer to run.
# When it finishes, it may complain about a "divergent transition". So long as
# there were only one or two, it's not a big problem. But, we can play around with
# step sizes to get rid of it if we wanted to.

mod_brm <- brm(Archaea ~ as.factor(Timepoint)*as.factor(Height) + (1|Dish/Sample), data = d, control =
                 list(adapt_delta = 0.98, max_treedepth = 15))

# Now, the warning should be gone. If it comes back, just set adapt_delta closer to 1.
# Next, we can check out some diagnostics

summary(mod_brm)
plot(mod_brm) # n.b. you need to press "enter" to show all plots

# R-hat is very close to 1 for all tests, and the plots look good.
# The plots show the distribution of each parameter, and traces the "chains" from
# the MCMC optimization. These are basically just independent replicates of the fitting
# process. Ideally, we want all four chains to look pretty similar, as this
# suggests that we have found a robust optimum.

# The paremters show, respectively:
# b_Intercept: overall model intercept, which we can more or less ignore
# b_as.factor(Timepoint2): difference between Timepoint1 and Timepoint 2 in hight class 1 - i.e. the time effect in the first height class
# b_as.factor(Height2): difference between Height1 and Height2, within Timepoint1 - i.e. the height effect in the first time period
# b_as.factor(Timepoint2):b_as.factor(Height2): interaction between time and height. You can either interpret this as the
# change in the time effect in height class 2 vs. height class 1, or the change in the height effect in time class 2 vs. time class 1.
# Both interpretations are equally valid, and can't really be separated, since the joint treatments are unseparable.
# sd_Dish_Intercept: Random effect -- can be ignored for our purposes. Shows variance in intercepts associated with dish identity (bigger numbers imply larger differences among dishes)
# sd_Sample_Intercept: Same as above, but for the variation among samples within a dish
# sigma: residual error term (i.e. what you would use to calculate deviance, residual squared error, etc)

# Finally, we can go ahead and test each of our parameters. The summary command shows 95% CI's -- we can calculate these by had too.
# First, we extract the estimates from the chains
parameter_estimates <- posterior_samples(mod_brm)

# Then we calculate 95% CI's
quantile(parameter_estimates$b_as.factorTimepoint2, c(0.025, 0.975))
quantile(parameter_estimates$b_as.factorHeight2, c(0.025, 0.975))
quantile(parameter_estimates$`b_as.factorTimepoint2:as.factorHeight2`, c(0.025, 0.975))

# At least for my analysis, I find that the 95% CI for time does not include zero,
# the CI for the height factor DOES include, zero, and for the interaction term does not include zero.

# This means that at alpha = 0.05, the time effect is significant (i.e. within the first height class, time 1 and time 2 differ significantly)
# whereas the height effect is not (i.e. the difference between height 1 and height 2 within time step 1 is not significantly different from zero).
# Lastly, the way that I would interpret the interaction would be that the time effect is significantly different in the height 2 class than in the height 1 class --
# that is, the time effect is significantly more *negative* in the height 2 region than in the height 1 region.

# We can also directly extract p-values for these as follows:
mean(parameter_estimates$b_as.factorTimepoint2 < 0)
mean(parameter_estimates$b_as.factorHeight2 > 0)
mean(parameter_estimates$`b_as.factorTimepoint2:as.factorHeight2` > 0)

# Note that the height effect is very close to p = 0.05.
# The code above is a bit hacky because we need to know what the dominant sign of each parameter is, in order
# to set the direction of the inequality. We can make this easier with a simple function

getpval <- function(x) {p = mean(x > 0); out = pmax(pmin(p, 1-p), 1/length(x)); if(all(x==0)) {out = 1}; out}

getpval(parameter_estimates$b_as.factorTimepoint2)
getpval(parameter_estimates$b_as.factorHeight2)
getpval(parameter_estimates$`b_as.factorTimepoint2:as.factorHeight2`)


# In publications, I would suggest just reporting the 95% CI for these three variables, along with the R-hats.
# You should be able to use the same code as above to fit and test the other two models.
# Breifly, this gives us:

mod_brm2 <- brm(Bacteria ~ as.factor(Timepoint)*as.factor(Height) + (1|Dish/Sample), data = d, control =
                 list(adapt_delta = 0.982, max_treedepth = 15))
summary(mod_brm2)
parameter_estimates2 <- posterior_samples(mod_brm2)
getpval(parameter_estimates2$b_as.factorTimepoint2)
getpval(parameter_estimates2$b_as.factorHeight2)
getpval(parameter_estimates2$`b_as.factorTimepoint2:as.factorHeight2`)
# so, no significant time or time*height interaction, but there is a significant height effect


mod_brm3 <- brm(NOB ~ as.factor(Timepoint)*as.factor(Height) + (1|Dish/Sample), data = d, control =
                  list(adapt_delta = 0.98, max_treedepth = 15))
summary(mod_brm3)
parameter_estimates3 <- posterior_samples(mod_brm3)
getpval(parameter_estimates3$b_as.factorTimepoint2)
getpval(parameter_estimates3$b_as.factorHeight2)
getpval(parameter_estimates3$`b_as.factorTimepoint2:as.factorHeight2`)
# so, no significant time effect for height 1, but a significant height effect, and a significant
# change in the time effect in height class 2.


########################################################
# Extended analyses, giving contrast for a, b, c, d, etc notation in bar plots.
########################################################

# Archaea
aH1T1 = parameter_estimates$b_Intercept
aH1T2 = parameter_estimates$b_Intercept + parameter_estimates$b_as.factorTimepoint2
aH2T1 = parameter_estimates$b_Intercept + parameter_estimates$b_as.factorHeight2
aH2T2 = parameter_estimates$b_Intercept + parameter_estimates$b_as.factorTimepoint2 + parameter_estimates$b_as.factorHeight2 + parameter_estimates$`b_as.factorTimepoint2:as.factorHeight2`

# check that contrasts were done correctly (values should be similar)
mean(aH1T1); mean(d$Archaea[d$Height == 1 & d$Timepoint == 1])
mean(aH1T2); mean(d$Archaea[d$Height == 1 & d$Timepoint == 2])
mean(aH2T1); mean(d$Archaea[d$Height == 2 & d$Timepoint == 1])
mean(aH2T2); mean(d$Archaea[d$Height == 2 & d$Timepoint == 2])
# looks good

# calculate significance of pairwise differences
contrmat = data.frame(aH1T1, aH1T2, aH2T1, aH2T2)
pvalmat = matrix(nrow=4, ncol=4)
rownames(pvalmat) = colnames(pvalmat) = colnames(contrmat)

for(i in 1:ncol(contrmat)) {
  for(j in i:ncol(contrmat)) {
    pvalmat[i,j] = getpval(contrmat[,i]-contrmat[,j])
  }
}


# Bacteria
bH1T1 = parameter_estimates2$b_Intercept
bH1T2 = parameter_estimates2$b_Intercept + parameter_estimates2$b_as.factorTimepoint2
bH2T1 = parameter_estimates2$b_Intercept + parameter_estimates2$b_as.factorHeight2
bH2T2 = parameter_estimates2$b_Intercept + parameter_estimates2$b_as.factorTimepoint2 + parameter_estimates2$b_as.factorHeight2 + parameter_estimates2$`b_as.factorTimepoint2:as.factorHeight2`

# check that contrasts were done correctly (values should be similar)
mean(bH1T1); mean(d$Bacteria[d$Height == 1 & d$Timepoint == 1])
mean(bH1T2); mean(d$Bacteria[d$Height == 1 & d$Timepoint == 2])
mean(bH2T1); mean(d$Bacteria[d$Height == 2 & d$Timepoint == 1])
mean(bH2T2); mean(d$Bacteria[d$Height == 2 & d$Timepoint == 2])
# looks good

# calculate significance of pairwise differences
contrmat2 = data.frame(bH1T1, bH1T2, bH2T1, bH2T2)
pvalmat2 = matrix(nrow=4, ncol=4)
rownames(pvalmat2) = colnames(pvalmat2) = colnames(contrmat2)

for(i in 1:ncol(contrmat2)) {
  for(j in i:ncol(contrmat2)) {
    pvalmat2[i,j] = getpval(contrmat2[,i]-contrmat2[,j])
  }
}



# NOB
nH1T1 = parameter_estimates3$b_Intercept
nH1T2 = parameter_estimates3$b_Intercept + parameter_estimates3$b_as.factorTimepoint2
nH2T1 = parameter_estimates3$b_Intercept + parameter_estimates3$b_as.factorHeight2
nH2T2 = parameter_estimates3$b_Intercept + parameter_estimates3$b_as.factorTimepoint2 + parameter_estimates3$b_as.factorHeight2 + parameter_estimates3$`b_as.factorTimepoint2:as.factorHeight2`

# check that contrasts were done correctly (values should be similar)
mean(nH1T1); mean(d$NOB[d$Height == 1 & d$Timepoint == 1])
mean(nH1T2); mean(d$NOB[d$Height == 1 & d$Timepoint == 2])
mean(nH2T1); mean(d$NOB[d$Height == 2 & d$Timepoint == 1])
mean(nH2T2); mean(d$NOB[d$Height == 2 & d$Timepoint == 2])
# looks good

# calculate significance of pairwise differences
contrmat3 = data.frame(nH1T1, nH1T2, nH2T1, nH2T2)
pvalmat3 = matrix(nrow=4, ncol=4)
rownames(pvalmat3) = colnames(pvalmat3) = colnames(contrmat3)

for(i in 1:ncol(contrmat3)) {
  for(j in i:ncol(contrmat3)) {
    pvalmat3[i,j] = getpval(contrmat3[,i]-contrmat3[,j])
  }
}



pvalmat  # Archaea
pvalmat2 # Bacteria
pvalmat3 # NOB




########################################################
# Check model fit
########################################################
summary(mod_brm)
summary(mod_brm2)
summary(mod_brm3)


save.image(file = "cardFISH.Rdata")


