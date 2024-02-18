# Modified by Alie Fink-Skular from check_gambling_behav_stim.r script from Ignacio Saez
# Last modified: 09/06/2023
# Inputs preprocessed data file from preprocess_gamble_stim_AF.r script

library(ggplot2)
library(reshape)
library(plyr)
library(scales)
library(stats)

# SPECIFY PROJECT DIRECTORY, SUBJECT INFO, DATA DIRECTORIES, SAVE DIRECTORIES


# Quick examination of gamble proportion by win probability, for individual subjects


################# BASELINE BEHAVIOR ################
behav <- read.csv(baseline, stringsAsFactors=FALSE)

win.prob.seq <- seq(0,1,0.1)
gamble.prop <- 0
choice <- behav$choice.class
win.prob <- behav$win.prob
# Proportion of gamble choices by win prob
for(j in 1:length(win.prob.seq)){
  gamble.prop[j] <- sum(choice[win.prob == (j-1)/10] == 'Gamble') / # (j-1)/10 is the win probability
    (sum(choice[win.prob == (j-1)/10] == 'Gamble') + sum(choice[win.prob == (j-1)/10] == 'Safebet')) # normalized by total number of trials with that winprob
}

# Fit
choices = as.numeric(behav$choice.class == 'Gamble')
g=glm(choices~win.prob,family=binomial) # run a logistic regression model (in this case, generalized linear model with logit link)

################# STIM BEHAVIOR ################
behav.stim <- read.csv(stim, stringsAsFactors=FALSE)
win.prob.seq <- seq(0,1,0.1)
gamble.prop.stim <- 0
choice.stim <- behav.stim$choice.class
win.prob.stim <- behav.stim$win.prob
# Proportion of gamble choices by win prob
for(j in 1:length(win.prob.seq)){
  gamble.prop.stim[j] <- sum(choice.stim[win.prob.stim == (j-1)/10] == 'Gamble') / # (j-1)/10 is the win probability
    (sum(choice.stim[win.prob.stim == (j-1)/10] == 'Gamble') + sum(choice.stim[win.prob.stim == (j-1)/10] == 'Safebet'))
}

# Fit
choices.stim = as.numeric(behav.stim$choice.class == 'Gamble')
g.stim=glm(choices.stim~win.prob.stim,family=binomial) # run a logistic regression model (in this case, generalized linear model with logit link)



################### Plotting & save fig ####################

pdf = TRUE
if(pdf) {
  pdf(paste0(save.path,subj,'_stim_behav.pdf'),height=3.5,width=3.5)
}
par(mar=c(4,4,1,1))
plot(win.prob.seq,gamble.prop,las=1,bty='l',ylab='Proportion of gambles',xlab='Gamble win probability',ylim=c(0,1),pch=16)
points(jitter(win.prob,1),jitter(choices,0.1),col=alpha('black',.3),pch=16,cex=.5)
curve(predict(g,data.frame(win.prob=x),type="resp"),add=TRUE,col='black')
points(win.prob.seq,gamble.prop,pch=16)

points(win.prob.seq,gamble.prop.stim,pch=16,col='dodgerblue3',lwd=2)
points(jitter(win.prob.stim,1),jitter(choices.stim,0.1),col=alpha('dodgerblue3',.3),pch=16,cex=.5)
curve(predict(g.stim,data.frame(win.prob.stim=x),type="resp"),add=TRUE,col='dodgerblue3')
points(win.prob.seq,gamble.prop,pch=16)

legend(c('Baseline','Stim'),x='topleft',bty='n',col=c('black','dodgerblue3'),pch=c(16,16),lwd=1)
if(pdf) {
  dev.off()
}




stim.curve_fit = g.stim$fitted.values
plot(win.prob.stim,stim.curve_fit)

base.curve_fit = g$fitted.values
plot(win.prob,base.curve_fit)

fisher.test(base.curve_fit,stim.curve_fit,simulate.p.value=TRUE)

base.predict = predict(g,data.frame(win.prob),type="resp")
stim.predict = predict(g.stim,data.frame(win.prob),type="resp")

fisher.test(base.predict,stim.predict,simulate.p.value=TRUE)
