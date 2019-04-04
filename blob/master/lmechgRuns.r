#read the raw data
source("simData.dump")

#look at the raw data
plot(groupedData(FVC~Age|Subj,data=d))

# Standard LME run, with fixed breakpoints set at 18
dat18 <- setbreakpoint(d,18,rep(1,dim(d)[1]),"Age","Subj")
m.18 <- lme(FVC~I(Lo*Age)+Lo+I(Hi*Age)+Hi-1,data=dat18,random=~I(Lo*Age)+Lo+I(Hi*Age)+Hi-1|Subj,method="ML")

# Augmented lme with changepoints at 12 15 18 21 24 using m.18 for initial parmameters
m.mix15.r6.b5 <- lmeChgPt(fixed = FVC ~ I(Lo*Age)+Lo+I(Hi*Age)+Hi-1,
                  data = d,
                  random = ~ I(Lo*Age)+Lo+I(Hi*Age)+Hi-1|Subj,
                  brkvar="Age",
                  bkpts=c(12,15,18,21,24),
                  initProbs=c(.2,.2,.2,.2,.2),
                  initParams=list(fixed=fixed.effects(m.18),
                                  random=VarCorr(m.18)),
                  initChs=rep(3,50),
                  nRepl=6,
                  max.fail=1,
                  max.iter=8,
                  thresh=1e-3,
                  verbose=TRUE)

# It is best to use summary to view results; the result from lmeChgPt contains "replicated" data, so 
# inference needs to be adjusted - that's what summary function does.
summary(m.mix15.r6.b5) 
