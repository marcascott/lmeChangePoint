require(nlme)
require(mvtnorm)

# =============================================================================
# Fit a change point lme model
# =============================================================================
lmeChgPt <- function(fixed, data, random, brkvar, bkpts, initProbs, initParams, initChs,
                     nRepl=10, max.fail=5, max.iter=5, thresh=1e-3, opt=F, verbose=T,...)
{

   # ====================================================
   # process calling parameters
   # ====================================================

   # match the call for storage in class list
   mcall <- match.call()
   
   # check for fixed effects formula, try to get from groupedData object if missing
   if (missing(fixed)) {
     # if data is grouped data get formula
   }
   
   # check for random effects formula, set to intercept only if missing
   if (missing(random)) {
       random=~1
   }
   
   # breakpoints must be specified
   if (missing(bkpts)) {
       stop('You must supply an array of change points!')
   }
   
   # if initprob is missing assume uniform distribution
   if (missing(initProbs)) {
       initProbs <- vector("numeric",length(bkpts))
       initProbs <- rep(1/length(bkpts),3)
   }
   
   # get the information on the formulas
   fixedterms   <- attr(terms(formula(fixed)),"variables")
   fixedlabels  <- attr(terms(formula(fixed)),"term.labels")
   randomterms  <- attr(terms(formula(random)),"variables")
   randomlabels <- attr(terms(formula(random)),"term.labels")
   
   # get the grouping variable
   reStructRandom <- reStruct(random, REML = FALSE, data = NULL)
   groups <- attr(terms(getGroupsFormula(reStructRandom)),"term.labels")[[1]]

   #get groups if missing and data is groupedData object
   if (missing(groups)) {
      if (inherits(data,"groupedData")) {
         groups <- attr(terms(getGroupsFormula(data)),"term.labels")[[1]]
      } else {
         stop("Grouping must be specified in random effects formula")
      }
   }
   
   # brk var set to first independendant var
   if (missing(brkvar)) {
       brkvar <- fixedlabels[[2]]
   }

   # if initParams are missing run an LME using fixed and random specified
   if (missing(initParams)) {
       # run a model using midpoint prob if odd, random of  midpoints if even
   }

   # set up formulae for model matrix calls
   fixmat <- fixed
   ranmat <- formula(reStructRandom)[[1]]

   #update the data frame to include a second ID that stays fixed
   data <- cbind(data,dupid=data[,groups])

   # ====================================================
   # begin fitting algorithm
   # ====================================================

   #EM housekeeping initializations
   n        <- length(initChs) # count one per subject
   nComp    <- length(initProbs)
   fcalls   <- matrix(NA,n,nComp)
   p        <- matrix(NA,n,nComp)
   ct       <- 1
   max.iter <- max.iter + 1  # allows for "0-th" startup parameters & lliks to be stored in rec

   #record of each iteration stored in rec
   rec <- list(mixProbs=matrix(NA,max.iter,nComp),params = vector("list",length=nComp),
               loglik = rep(NA,max.iter),
               mllik=rep(NA,max.iter),delta=rep(NA,max.iter))

   #initialize rec
   tmp.ps           <- matrix(0,n,nComp)
   tmp.ps[,initChs] <- 1 # assume that the choice was made "definitively"
   mixProbs         <- apply(tmp.ps,2,mean)

   tmp <- lmeChgPt.llik(data=data, fixed=fixmat, random=ranmat, mixProbs=mixProbs,
                        bkpts=bkpts, ps=tmp.ps, params=initParams,
                        breakvar=brkvar, groupvar=groups, verbose=T)

   rec$delta[ct]           <- delta <- Inf
   rec$mixProbs[ct,]      <- initProbs
   rec$params[[ct]]       <- initParams
   llik <- rec$loglik[ct] <- tmp$llik
   rec$mllik[ct]          <- tmp$mllik

   # information display
   if (verbose) {
      cat("In EM Algorithm, iteration ",ct," logliks= ",round(llik,4),
          " delta= ", round(rec$delta[ct],4), "\n")
      cat("Marginal llik = ",round(rec$mllik[ct],4),"\n")
      cat("mixProbs = (",round(mixProbs,4),")\n")
      cat("Current parameters:\n") #display something intelligent using params (a list)
      cat("Fixed: ", round(rec$params[[ct]]$fixed,5),'\n')
      print(rec$params[[ct]]$random)
   }

   # final set up for EM iteration
   newParams <- initParams
   chs       <- initChs
   mixProbs  <- initProbs

   # EM iteration loop
   while (ct < max.iter && delta > thresh) {

      # just for safety to prevent memory overflow
      gc()

      # improvment success count
      ct <- ct + 1

      # E-STEP: - - - - - - - - - - - - - - - - - -
      if (verbose) cat("E-Step...  ")

      for (j in 1:nComp) {
          #generate a separate dataframe for each component
          data.j     <- setbreakpoint(data=data, breakpoints=bkpts, breakchoice=rep(j,n),
                                      breakvar=brkvar, groupvar=groups)
          fcalls[,j] <- dlme(data=data.j, fixed=fixmat, random=ranmat,
                            params=newParams,groupvar=groups)
                            # dlme uses lmeParams to set the mean effs and covar struct,
                            # calls rmvnorm ...
      } # for j

      for (i in 1:n) {
         for (j in 1:nComp) {
            p[i,j] <- mixProbs[j]*fcalls[i,j]/sum(mixProbs*fcalls[i,]) #new probs of mbrshp
         } # for i
      } # for j

      if (verbose) cat("done\n")

      # M-STEP: - - - - - - - - - - - - - - - - - -
      mixProbs <- apply(p,2,mean)
      if (verbose) {
         cat("M-Step...  mixProbs = ",round(mixProbs,4),
             "\nRecord classification probs = ",
         round(apply((t(apply(p,1,rank))==nComp),2,mean),4)," ...  ")
      } # if verbose

      attempt <- 0
      delta0  <- -Inf
      mllik   <- -Inf

      while (attempt < max.fail & delta0<0) {

         # REPLICATION HERE:
         for (j in 1:nRepl) {

            for (i in 1:n) chs[i] <- sample(1:nComp,1,prob=p[i,])

            tmp        <- setbreakpoint(data=data.j, breakpoints=bkpts, breakchoice=chs,
                                        breakvar=brkvar, groupvar=groups)
            nms0       <- levels(tmp[,groups])
            nms        <- sapply(nms0,"paste",as.character(j-1),sep="")
            names(nms) <- NULL
            cts        <- table(tmp[,groups])
            levels(tmp[,groups]) <- nms
            if (j == 1) {
               data2 <- tmp
            } else {
               curLvls <- levels(data2)
               levels(data2) <- c(curLvls,levels(tmp[,groups]))
               data2 <- rbind(data2,tmp)
            }
         } # end for i

         lme.fit <- lme(fixed, data2, random=random, method="ML")

         # SAVE THIS FOR POSSIBLE USE LATER ON control=list(natUnconstrained = F))
         tmpParams <- list(fixed=fixed.effects(lme.fit), random=VarCorr(lme.fit))

         #compute new loglik; update rec, compute delta...
         if (verbose) cat("Computing Marginal Log Likelihood...\n")
         tmp <- lmeChgPt.llik(data,fixmat,ranmat,mixProbs,bkpts,p,tmpParams,
                              brkvar,groups,verbose)

         delta0 <- tmp$mllik-rec$mllik[ct-1]

         if (tmp$mllik > mllik) {
            mllik <- rec$mllik[ct] <- tmp$mllik
            if (verbose) cat("Marginal llik = ",round(mllik,4),"\n")
            llik <- rec$loglik[ct] <- tmp$llik
            newParams <- tmpParams
         }
         attempt <- attempt + 1
         
      } # end while

      if (verbose) cat("done...")

      if (delta0 < 0) {
         cat("Failed to find a better fit.\n")
         #roll back record
         rec$params[[ct]]  <- rec$params[[ct-1]]
         rec$mixProbs[ct,] <- rec$mixProbs[ct-1,]
         rec$loglik[ct]    <- rec$loglik[ct-1]
         rec$mllik[ct]     <- rec$mllik[ct-1]
         rec$delta[ct]     <- rec$delta[ct-1]
      } else {
         cat("Successfully found improved fit.\n")
         #update rec, compute delta...
         rec$params[[ct]] <- list(fixed=fixed.effects(lme.fit), random=VarCorr(lme.fit))
         # NOTE: does this get error
         # Struct-such as sigma^2_epsilon
         rec$mixProbs[ct,] <- mixProbs
         if (ct>1) rec$delta[ct] <- delta <- rec$mllik[ct]-rec$mllik[ct-1]
         ##ADD## RGN 12/22/01
         last.data  <- data2
         last.model <- lme.fit
      }

      if (verbose) {
          cat("In EM Algorithm, iteration ",ct," logliks= ",round(llik,4),
              " delta= ", round(rec$delta[ct],4), "\n")
          cat("Marginal llik = ",round(rec$mllik[ct],4),"\n")
          cat("mixProbs = (",round(mixProbs,4),")\n")
          cat("Current parameters:\n") #display something intelligent using params (a list)
          cat("Fixed: ", round(rec$params[[ct]]$fixed,5),'\n')
          print(rec$params[[ct]]$random)
      }

   } # end of EM iteration loop
   
   if (opt) {
       # run optimization code to get standard errors
       # reset opt parameter to results
   }
   
   # set up result object
   res <- list(call=mcall,
               group=groups,
               changepts=bkpts,
               initprob=initProbs,
               history=rec,
               pMembshp = p,
               mixProbs=mixProbs,
               recClasif=apply((t(apply(p,1,rank))==nComp),2,mean),
               params=rec$params[[ct]],
               last.data=last.data,
               last.model=last.model,
               optim=opt)
   class(res) <- "lmechg"

   # return the resulting object
   return(res)
   
}




# ==============================================================================
# LOG LIKELIHOOD ESTIMATION
# ==============================================================================
lmeChgPt.llik <- function(data,fixed,random,mixProbs,bkpts,ps,params,breakvar,
                          groupvar,verbose=T) {

    #EM housekeeping initializations
    nComp  <- length(mixProbs)
    subjs  <- dim(ps)[1]
    n      <- dim(data)[1]
    fcalls <- matrix(NA,subjs,nComp)

    #generate a separate dataframe for each component
    for (j in 1:nComp) {
        datafull   <- setbreakpoint(data=data, breakpoints=bkpts, breakchoice=rep(j,n),
                                    breakvar=breakvar, groupvar=groupvar)
        fcalls[,j] <- dlme(data=datafull, fixed=fixed, random=random,
                           params=params, groupvar=groupvar)
    }

    # return the results
    list(mllik=sum(log(apply(t(fcalls)*mixProbs,2,sum))),
         llik=sum(log(fcalls)*ps))

}




# =============================================================================
# Multivariate Normal Sample
# =============================================================================
dlme <- function(data,fixed,random,params,groupvar)
{

    # get the model frames we need
    X        <- model.matrix(fixed,data)
    meanV    <- X%*%params$fixed
    covParms <- convVarCorr(params$random)
    Z0       <- model.matrix(random,data)
    
    # set local variables
    tmp      <- getGroupCounts(data=data,fixed=fixed,groupvar=groupvar)
    lns      <- tmp$lns
    obs      <- tmp$obs
    n        <- length(lns)
    cutpts   <- c(0,cumsum(lns))
    r        <- rep(0,n)

    # compute e-hat, then find coresponding dmvnorm
    for (i in 1:n) {
       rng  <- (cutpts[i]+1):cutpts[i+1]
       Z    <- Z0[rng,]
       R    <- diag(lns[i])*covParms$var
       obs0 <- obs[rng]
       mean0<- as.vector(meanV[rng])
       if (lns[i]==1) Z <- t(as.matrix(Z))  #this is necessary b/c S converts matrix to vector
       r[i] <- dmvnorm(obs0,mean0,Z%*%covParms$cov%*%t(Z)+R)
    }

    return(r)

}




# checked 11/19/03 RGN
# =============================================================================
# Set change point indicator variables
#      NOTE: assumes indicators already exist in the data.frame
# =============================================================================
setbreakpoint <- function ( data, breakpoints, breakchoice, breakvar, groupvar)
{

   # get a local data.frame with a break point element
   tempdata <- data
   datlen   <- dim(tempdata)[1]
   
   # get a list of subjects factor names
   groups <- levels(tempdata[,groupvar])

   # get a limit for the loop
   looplim  <- length(groups)

   # if tempdata$brk does not exist create it
   if (is.null(tempdata$brk)) {
       tempdata <- cbind(tempdata,brk=rep(NA,datlen))
   } else {
       tempdata$brk <- rep(NA,datlen)
   }

   # loop through all subjects
   for (i in 1:looplim) {

      # get subject name and probability
      s <- groups[i]

      # set breakpoint to a randomly selected value
      tempdata[tempdata[,groupvar]==s,"brk"] <- breakpoints[breakchoice[i]]

   }

   # if Lo and Hi don't exist then create them
   if (is.null(tempdata$Lo)) {
      tempdata <- cbind(tempdata,Lo=NA)
   }
   if (is.null(tempdata$Hi)) {
      tempdata <- cbind(tempdata,Hi=NA)
   }
   
   # now set them base on the age
   tempdata$Lo <- as.numeric(tempdata[,breakvar]<tempdata[,"brk"])
   tempdata$Hi <- 1-tempdata$Lo

   # return the data frame
   tempdata

}




# checked 11/19/03 RGN
# =============================================================================
# Get the count of the number of subjects
# =============================================================================
getGroupCounts <- function(data,fixed,groupvar)
{

    # data is a data frame
    # fixed is the fixed effects formula
    
    # get the counts per subject from a split data frame
    cts <- lapply(split(data[,groupvar],data[,groupvar]),length)

    # get the response values
    resp <- data[,paste(attr(terms(formula(fixed)),"variables")[[2]])]

    # return the count vector & the reponse vector
    #   list(lns = cts, obs=eval(parse(text=resp)))
    list(lns = cts, obs = resp)
    
}




# unchanged 11/19/03 RGN
# =============================================================================
# Convert Variance Covariance Matrix
# ==============================================================================
convVarCorr <- function (varcor)
{

    # Get the VarCorr structure
    varlist <- varcor

    # maxlen is the number of entries in variance column (= 1 resid + terms in model)
    maxlen  <- length(varlist[,1])
    numvar  <- maxlen-1

    # get the names
    n <- dimnames(varlist)[[1]]

    # get the variance array and sd array
    v <- as.numeric(varlist[,1])
    names(v) <- n
    s <- as.numeric(varlist[1:numvar,2])

    # create a starting covariance matrix (with 1/2 on diag)
    cv <- diag(rep(0.5,numvar))
    rownames(cv) <- n[1:numvar]
    colnames(cv) <- n[1:numvar]

    # now the ugly part to extract the lower triangular correlation data
    for (i in 1:(numvar-1)) {
        cv[(i+1):numvar,i] <- as.numeric(varlist[(i+1):numvar,i+2])
    }

    # convert lower triangular to complete correlation matrix (diag now 1)
    cv <- cv+t(cv)

    # convert to covariance matrix
    cv <- cv*(s%*%t(s))

    # return the result
    list(cov=cv, var=v[maxlen])

}




# =============================================================================
# Support function for par parameter convert
# ==============================================================================
lowVectToMat <- function(v)
{

#
#   m <- matrix(0,4,4)
#   m[1,1] <- v[1]
#   m[2,2] <- v[2]
#   m[3,3] <- v[3]
#   m[4,4] <- v[4]
#   m[2,1] <- v[5]
#   m[3,1] <- v[6]
#   m[4,1] <- v[7]
#   m[3,2] <- v[8]
#   m[4,2] <- v[9]
#   m[4,3] <- v[10]
#
#   diag(m) <- exp(diag(m))
#   m <- m %*% t(m)
#   m
#

}




# =============================================================================
# Convert paramters to vector for optimization
# ==============================================================================
convParamToPar <- function (params,doChol=TRUE)
{

#
#    # get the values
#    t      <- convVarCorr(params$random)
#    varcov <- t$cov
#    resid  <- t$var
#    fixed  <- params$fixed
#
#    # remove the strings
#    rownames(varcov) <- NULL
#    names(fixed)     <- NULL
#    names(resid)     <- NULL
#
#    if (doChol==TRUE) {
#        varcov <- chol(varcov)
#        diag(varcov) <- log(diag(varcov))
#        varcov <- t(varcov)
#    }
#
#    varcov <- as.vector(varcov[lower.tri(varcov,diag=TRUE)])
#
#    #create the vector (fixed,resid,covariance)
#    c(fixed,resid,varcov[1],varcov[5],varcov[8],varcov[10],
#      varcov[2],varcov[3],varcov[4],varcov[6],varcov[7],
#      varcov[9])

}




# =============================================================================
# Convert optimization vector to param structure
# ==============================================================================
convParToParam <- function(par)
{

#    list(fixed=par[1:5],var=par[6],random=lowVectToMat(par[7:16]))

}




# =============================================================================
# Calculate the marginal log likelihood
# ==============================================================================
mllik <- function(par,fixed,random,data,mixProbs,bkpts)
{

#
#    # get the variance covariance matrix
#    params <- convParToParam(par)
#
#
#    #EM housekeeping initializations
#    nComp  <- length(mixProbs)
#    n      <- dim(data)[1]
#    subjs  <- length(levels(data$Subj))
#    fcalls <- matrix(NA,subjs,nComp)
#
#    #generate a separate dataframe for each component
#    for (j in 1:nComp) {
#        datafull   <- setbreakpoint3(data=data,bkpt=bkpts,choice=rep(j,n))
#        fcalls[,j] <- moddlme(data=datafull,fixed=fixed,random=random,params=params)
#    }
#
#    # return the results
#    list(mllik=sum(log(apply(t(fcalls)*mixProbs,2,sum))))
#

}




# =============================================================================
# Multivariate normal likelihoood for data
# ==============================================================================
moddlme <- function(data,fixed,random,params)
{

#
#    # get model matrices
#    X  <- createdatamatrix(data,fixed)
#    Z0 <- createdatamatrix(data,fixed)[,-5] #MODIFY *** FOR Z matrix!!!
#
#    # get parameters from modified parameter list
#    meanV    <- X%*%params$fixed
#    covParms <- params$random
#    pvar     <- params$var
#
#    # set up locals
#    tmp    <- getSubjCts(data,fixed,random)
#    lns    <- tmp$lns
#    obs    <- tmp$obs
#    dfn    <- length(lns)
#    cutpts <- c(0,cumsum(lns))
#    r      <- rep(0,dfn)
#
#    # calculate
#    for (i in 1:dfn) {
#      rng   <- (cutpts[i]+1):cutpts[i+1]
#      Z     <- Z0[rng,]
#      R     <- diag(lns[i])*pvar
#      obs0  <- obs[rng]
#      mean0 <- as.vector(meanV[rng])
#      if (lns[i]==1) {
#        Z   <- t(as.matrix(Z))  #this is necessary b/c S converts matrix to vector
#      }
#      r[i]  <- dmvnorm(obs0,mean0,Z%*%covParms%*%t(Z)+R)
#    }
#    r
#

}




# ******************************************************************************
# SUPPORT ROUTINES TO BE IMPLEMENTED FOR USER OUTPUT
# ******************************************************************************



# OK 11/25/03 RGN
# ==============================================================================
# Summary function for class lmechg (to not rely on lme output)
# ==============================================================================

# need to change to use optim results if available for standard errors

summary.lmechg <- function(obj)
{
   # @ @ @ @ need to modify to return an object @ @ @ @
   # to report
   # if optimization is run use those standard errors etc.
   callparm <- as.list(obj$call)
   nRepl <- callparm$nRepl
   fix  <- obj$last.model$coefficients$fixed
   ran  <- VarCorr(obj$last.model)
   lenf <- length(fix)
   lenr <- dim(ran)[1]-1
   nmf  <- names(fix)
   fixnm<- abbreviate(nmf,14)
   vf   <- obj$last.model$varFix
   se   <- sqrt(diag(vf)*nRepl)
   obs  <- dim(obj$last.model$groups)[1]
   grp  <- dim(unique(obj$last.model$groups))[1]
   df   <- (obs-grp)/nRepl-(lenf-1)
   tval <- fix/se
   sig  <- 2*(1-pt(abs(tval),df))
   sd   <- sqrt(diag(vf))
   vf   <- vf/sd%*%t(sd)
   ll   <- -2*getLogLik.lmechg(obj)
   aic  <- getAIC.lmechg(obj)
   bic  <- getBIC.lmechg(obj)
   
   # configuration
   out1 <- data.frame(data=as.character(callparm$data),
                      brkvar = as.character(callparm$brkvar),
                      replic = as.integer(nRepl),
                      iter   = as.integer(callparm$max.iter),
                      obs    = as.integer(obs/nRepl),
                      grp    = as.integer(grp/nRepl))

   # Goodness of Fit
   out2 <-  data.frame(LL2=ll, AIC=aic, BIC=bic)
                
   # random effects
   out3 <- list(form   = as.character(callparm)[[4]],
                varcor = ran)
                
   # fixed effects
   out4 <- list(form   = as.character(callparm)[[2]],
                result = data.frame( fixnam = fixnm,
                                     fix    = fix,
                                     se     = se,
                                     df     = df,
                                     tval   = tval,
                                     sig    = sig),
                corr   = vf)

   # change points
   out5 <- data.frame(names=as.character(eval(callparm$bkpts)),
                      probs=obj$mixProbs)

   # notes
   out6 <- ifelse(obj$optim==F,
                  "NOTE: Standard errors estimated for full data problem",
                  "NOTE: Standard errors estimated by optimization")

   # configure output list objecgt
   res <- list(config=out1,gfit=out2,random=out3,fixed=out4,bkpt=out5,note=out6)
   class(res) <- c("summary.lmechg")
   
   # return the summary object
   return(res)
   
}




# OK 11/25/03 RGN
# ==============================================================================
# Print the summary object passed (usually implicitly (generic override)
# ==============================================================================
print.summary.lmechg <- function(obj)
{

   # title
   cat("VARIABLE CHANGEPOINT AUGMENTED LINEAR MIXED EFFECTS MODEL\n\n")

   # Configuration of model run - method, nRepl etc
   dat <- unlist(obj$config)
   cat(sprintf("   Data        : %-15s     Break Variable  : %-15s\n",
               as.character(obj$config[1,1]),as.character(obj$config[1,2])))
   cat(sprintf("   Replicates  : %4d                Iterations      : %3d\n",
               dat[3],dat[4]))
   cat(sprintf("   Observations: %4d                Groups          : %3d\n",
               dat[5],dat[6]))

   # goodness of fit
   cat("\n   GOODNESS OF FIT\n")
   dat <- unlist(obj$gfit)
   cat(sprintf("%28s %10s %10s\n","-2LL","AIC","BIC"))
   cat(sprintf("%28.3f %10.3f %10.3f\n",dat[1],dat[2],dat[3]))

   # random effect
   cat(paste("\n   RANDOM EFFECTS  ",as.character(obj$random$form),sep=""))
   print.VarCorr.lmechg(obj$random$varcor)

   # fixed effects
   tmp   <- obj$fixed$result
   fixnm <- abbreviate(tmp$fixnam,17)
   cat(paste("\n   FIXED EFFECTS  ",as.character(obj$fixed$form),sep=""))
   cat(sprintf("\n   %15s %10s %10s %6s %8s %8s \n","Variable","Coeff","SE","DF",
               "T","Sig"))
   for (i in 1:dim(tmp)[1]) {
      cat(sprintf(" %17s %10.4f %10.4f %7.0f %8.3f %8.4f \n",fixnm[i],tmp[i,2],
                  tmp[i,3],tmp[i,4],tmp[i,5],tmp[i,6]))
   }

   # correlations of fixed effects
   cat("\n       Correlation")
   vf <- obj$fixed$corr
   mnames <- dimnames(vf)[[1]]
   fixnm <- abbreviate(mnames,8)
   for (i in 1:(dim(vf)[1]-1)) {
       cat(sprintf(" %10s",fixnm[i]))
   }
   for (i in 2:dim(vf)[1]) {
      cat(sprintf("\n      %12s",fixnm[i]))
      for (j in 1:(i-1)) {
          cat(sprintf(" %10.3f",vf[i,j]))
      }
   }

   #   mixing probs
   bkpts <- unlist(obj$bkpt$names)
   probs <- unlist(obj$bkpt$probs)
   len   <- length(bkpts)
   cat("\n\n   MIXING PROBABTILITES\n")
   cat("         Break Points")
   for (i in 1:len) {
       cat(sprintf("%8s ",as.character(bkpts[i])))
   }
   cat("\n")
   cat("         Probability ")
   for (i in 1:len) {
       cat(sprintf("%9.3f",probs[i]))
   }

   cat(paste("\n\n   ",obj$note,"\n\n",sep=""))

}




# OK 11/25/03 RGN
# =============================================================================
# VarCorr Random effects retrieval
# ==============================================================================
VarCorr.lmechg <- function(obj,corr=F)
{

   # Get correlation matrix or covariance matrix from last.model
   if (corr) {
      out <- convVarCorr(VarCorr(obj$last.model))
   } else {
      out <- VarCorr(obj$last.model)
   }

   # set the class and return it
   class(out) <- "VarCorr.lmechg"
   out

}



# OK 11/25/03 RGN
# =============================================================================
# Print Variance Covariance Matrix  - normally by default (generic override)
# ==============================================================================
print.VarCorr.lmechg <- function (x)
{

   # Set up the names and headers
   rown <- abbreviate(dimnames(x)[[1]],15)
   coln <- abbreviate(dimnames(x)[[2]],11)
   cat("\n                  ")
   cat(sprintf("%11s",coln[1]))
   cat(sprintf("%11s",coln[2]))
   cat(sprintf("%11s",coln[3]))

   # output each row of the matrix
   for (i in 1:length(rown)) {
       # row name
       cat("\n   ")
       cat(sprintf("%15s",rown[i]))
       # print the entries on the row
       for (j in 1:length(rown)) {
           if (nchar(x[i,j])==0) {
              cat(sprintf("%-11s"," "))
           } else {
              # avoid coercion errors
              tmp <- suppressWarnings(as.numeric(x[i,j]))
              if (is.na(tmp)){
                 cat(sprintf("%11s",x[i,j]))
              } else {
                 cat(sprintf("%11.4f",tmp))
              }
           }
       }
   }

   cat("\n")

}




# =============================================================================
# Plot function for class lmechg (to not rely on lme output)
# ==============================================================================
plot.lmechg <- function(obj)
{

   #options - residual analysis
   #        - lattice plot with fixed / random / changepoint
   cat("\nNot yet implemented!\n\n")

}




# =============================================================================
# augPred function for class lmechg (to not rely on lme output)
# ==============================================================================
augPred.lmechg <- function(obj,bkpt="exp",bkat)
{

   # need to create interpolated values like lme

   callparm <- as.list(obj$call)
   primary  <- as.character(callparm$brkvar)
   response <- attr(terms(as.formula(callparm$fixed)),"variables")[[2]]
   nRepl    <- callparm$nRepl
   obs      <- dim(obj$last.model$groups)[1]
   obs      <- obs/nRepl
   grp      <- dim(unique(obj$last.model$groups))[1]
   group    <- obj$group
   dat      <- as.data.frame(obj$last.data[1:obs,])
   brks     <- eval(callparm$bkpts)
   mprob    <- obj$mixProbs

   # get breakpoint for prediction
   type <- pmatch(bkpt,c("exponent","fixed"))
   if (type==1) {
       bkat <- sum(brks*mprob)
   } else {
       if (missing(bkat)) {
          stop("Fixed breakpoint value must be specified!")
       }
   }
   
   # for plot function to produce effects lines
   dat  <- setbreakpoint(dat,c(bkat),rep(1,grp),primary,group)
   mm   <- model.matrix(as.formula(callparm$fixed),dat)
   resp <- mm %*% fixef(obj$last.model)
   
   out <- data.frame(dat[,primary],dat[,group],resp,rep("predicted",obs))
   
   names(out) <- c(primary,".group",response,".type")
   class(out) <- c("augPred","data.frame")
   return(out)
   
}




# =============================================================================
# Predict function for class lmechg (to not rely on lme output)
# ==============================================================================
predict.lmechg <- function(obj,bkpt="exp",bkat)
{

   # under construction
   # to return data.frame with 2 columns
   # predicted for fixed effects
   # predicted for blups
   # how to handle breakpoints to be decided (average, fixed, diff estimates?)
   callparm <- as.list(obj$call)
   primary  <- as.character(callparm$brkvar)
   nRepl    <- callparm$nRepl
   obs      <- dim(obj$last.model$groups)[1]
   obs      <- obs/nRepl
   grp      <- dim(unique(obj$last.model$groups))[1]
   group    <- obj$group
   dat      <- as.data.frame(obj$last.data[1:obs,])
   brks     <- eval(callparm$bkpts)
   mprob    <- obj$mixProbs

   # get breakpoint for prediction
   type <- pmatch(bkpt,c("exponent","fixed"))
   if (type==1) {
       bkat <- sum(brks*mprob)
   } else {
       if (missing(bkat)) {
          stop("Fixed breakpoint value must be specified!")
       }
   }

   # for plot function to produce effects lines
   dat  <- setbreakpoint(dat,c(bkat),rep(1,grp),primary,group)
   mm   <- model.matrix(as.formula(callparm$fixed),dat)
   resp <- mm %*% fixef(obj$last.model)

   out <- array(data=resp,dimnames=list(as.character(dat[,group])))

   return(out)

}



# OK 11/25/03 RGN
# =============================================================================
# Fixed effects retrieval
# ==============================================================================
fixef.lmechg <- function(obj)
{

   # get the fixed effects from the last model run
   out <- fixef(obj$last.model)

   # set names to shorter version
   names(out) <- abbreviate(names(out),10)

   # return fixed effects
   return(out)

}




# =============================================================================
# Random effects retrieval
# ==============================================================================
ranef.lmechg <- function(obj)
{

   # under construction
   callparm <- as.list(obj$call)
   nRepl    <- callparm$nRepl
   groups   <- obj$group

   tmp  <- obj$last.data
   tmp  <- split.data.frame(tmp,tmp$dupid)
   tmpl <- length(tmp)

   brkval <- array(NA, c(tmpl,nRepl))

   for (i in 1:tmpl) {
      tmp2 <- split.data.frame(tmp[[i]],tmp[[i]][,groups])
      for (j in 1:length(tmp2)) {
         tmp3 <- tmp2[[j]]
         brkval[i,j] <- tmp3[1,"brk"]
      }
   }

   # get random effect coefficients
   # ? how to combine replicates ?
   cat("\nNot yet implemented!\n\n")
   brkval
}




# OK 11/25/03 RGN
# =============================================================================
# get Log Likelihood
# ==============================================================================
getLogLik.lmechg <- function(obj)
{

   # get the marginal log likelhood for the last iteration
   ll   <- obj$history$mllik
   ll   <- ll[length(ll)]
   return(ll)

}




# OK 11/25/03 RGN
# =============================================================================
# get number of parameters for model
# ==============================================================================
getNumParm.lmechg <- function(obj)
{

   # get the number of fixed, random and breakpoint parameters
   lenf <- length(obj$last.model$coefficients$fixed)
   lenr <- dim(VarCorr(obj$last.model))[1]-1
   lenb <- length(eval(as.list(obj$call)$bkpts))
   
   # calculate the effective parameter number
   parms <- lenf + lenb+(lenr+lenr*lenr)/2
   return(parms)
   
}



# OK 11/25/03 RGN
# =============================================================================
# get AIC
# ==============================================================================
getAIC.lmechg <- function(obj)
{

   # claculate the AIC from model params and -2LL
   aic  <- -2*getLogLik.lmechg(obj)+2*getNumParm.lmechg(obj)
   names(aic) <- "AIC"
   return(aic)

}




# OK 11/25/03 RGN
# =============================================================================
# get BIC - first declare a generic then the lmechg function
# ==============================================================================
getBIC.lmechg <- function(obj)
{

   # calculate the BIC from sample size, model params and -2LL
   obs  <- dim(obj$last.model$groups)[1]
   bic  <- -2*getLogLik.lmechg(obj)+log10(obs)*getNumParm.lmechg(obj)
   names(bic) <- "BIC"
   return(bic)

}



# OK 11/25/03 RGN
# =============================================================================
# ANOVA function to compare two (or more) models
# ==============================================================================
anova.lmechg <- function(obj,...,test=T)
{

   ancall <- sys.call()
   ancall$verbose <- ancall$test <- NULL

   # initially only evaluate comparison between modles
   if (length(list(...))==0) {
      stop("Only a single model specified, please list at least 2 models.")
   }

   # all models must be lmechg objects
   object <- list(obj, ...)
   termsClass <- unlist(lapply(object, data.class))
   if (!all(termsClass=="lmechg")) {
      stop("Objects must inherit from class lmechg")
   }

   # get degrees of freedom for all models
   dfModel <- unlist(lapply(object, function(x)
      {
         nRepl <- as.list(x$call)$nRepl
         obs   <- dim(x$last.model$groups)[1]
         grp  <- dim(unique(x$last.model$groups))[1]
         lenf <- length(x$last.model$coefficients$fixed)
         df   <- (obs-grp)/nRepl-(lenf-1)
      }))
      
   # get log Liklihood, AIC and BIC for all models
   logLik <- unlist(lapply(object, getLogLik.lmechg))
   logLik <- logLik*-2
   AIC    <- unlist(lapply(object, function(x) getAIC.lmechg(x)))
   BIC    <- unlist(lapply(object, function(x) getBIC.lmechg(x)))

   modval <- data.frame(Model=c(1:length(unlist(AIC))), df=dfModel, "LL2"=logLik,
                        AIC=AIC, BIC=BIC)
   row.names(modval) <- unlist(lapply(as.list(ancall[-1]),deparse))

   # if test is to be done process it
   if (test) {

      # get df and don't test if zero
      dfDiff <- c(abs(diff(dfModel)))
      dfDiff <- sapply(dfDiff, function(x) ifelse(x==0,NA,x))

      # calculate likelihood ratio tests
      LLDiff <- c(diff(logLik))
      pval   <- 1-pchisq(abs(LLDiff),dfDiff)

      # set up test names
      nm <- array()
      for (i in 1:length(unlist(AIC))) {
         if (i==1) {
            nm[i] <- " "
         } else {
            nm[i] <- paste(i-1," vs ",i,sep="")
         }
      }
      
      # add test to output data frame
      modval <- data.frame(modval, test=nm, df2=c(NA,dfDiff), Ratio=c(NA,LLDiff),
                           pval=c(NA,pval))
   }
   
   # return an object that will be assigned or printed
   class(modval) <- c("anova.lmechg","data.frame")
   return(modval)

}



# OK 11/25/03 RGN
# =============================================================================
# Print the ANOVA object (usually by default call) (generic override)
# ============================================================================
print.anova.lmechg <- function(modval)
{

    # get the dimensions and names
    rlen <- dim(modval)[1]
    clen <- dim(modval)[2]
    cnm <- colnames(modval)
    rnm <- rownames(modval)

    # output the common field headers
    cat(sprintf("%15s %2s %4s %6s %6s  %6s"," ",cnm[1],cnm[2],cnm[3],cnm[4],cnm[5]))
    
    # add test headers if done
    if (clen>5) {
       tname <- as.character(levels(modval[,"test"]))
       cat(sprintf(" %7s  %3s %6s %6s",cnm[6],cnm[7],cnm[8],cnm[9]))
    }
    cat("\n")
    
    # print the model rows
    for (i in 1:rlen) {
        # get the parameters on the row
        dd <- unlist(modval[i,])
        # output the common columns
        cat(sprintf("%17s  %2.0f %4.0f %6.2f  %6.2f  %6.2f",
            rnm[i],dd[1],dd[2],dd[3],dd[4],dd[5]))
        # if test done output the tests results
        if (clen>5) {
           if (!is.na(dd[7])) {
              cat(sprintf("%8s %3.0f %6.2f %6.4f",tname[i],dd[7],dd[8],dd[9]))
           }
        }
        # go to next row on page
        cat("\n")
    } # for

}




# =============================================================================
# ? others to add
# =============================================================================

