
###############################################################################
# JMLE:                                                                       #
#                                                                             #
# Joint maximum likelihood estimation of item parameters and attribute        #
# profiles in cognitive diagnostic models.                                    #
#                                                                             # 
# Inputs:                                                                     #
# (1) Y: a matrix of binary responses (1=correct, 0=incorrect). Rows          #
#        represent persons, and columns represent items.                      #
# (2) Q: the Q-matrix of the test. Rows represent items, and colums represent #
#        attributes.                                                          #
# (3) model: currently has three options, "DINA", "DINO", and "NIDA".         #
# (4) conv.crit.par: the critical value for the change in item parameter      #
#                    values to determine convergence                          #      
# (5) conv.crit.att: the critical value for the percentage of attributes that #
#                    are changed to determine convergence                     #
#                                                                             #
# Outputs:                                                                    #
# (1) item parameters estimates and standard errors for each item             #
# (2) attribute profiles for each examinee                                    #
#                                                                             #
###############################################################################


JMLE <- function(Y, Q, model="DINA", NP.method="Weighted", conv.crit.par=0.001, conv.crit.att=0.01, max.ite=100) {

  #####
  # 1 #
  ##### Check dimension consistency and convert data to the right formats 
    
  check <- NULL
  check <- CheckInput(Y, Q)  
  if (!is.null(check)) return(warning(check))
  
  Y <- as.matrix(Y)
  Q <- as.matrix(Q)
  
  #####
  # 2 #
  ##### Step 1: Initial model-free estimation of the ability pattern 
  
  if (model %in% c("DINA", "NIDA", "GNIDA"))
  {
    gate="AND"
  } else if (model == "DINO")
  {
    gate="OR"
  } else
  {
    return(warning("Model specification is not valid."))
  }
  
  alpha.est.NP <- AlphaNP(Y, Q, gate=gate, method=NP.method, wg=1, ws=1)
  alpha.est <- alpha.est.NP$alpha.est
  loss.matrix <- alpha.est.NP$loss.matrix
  
  #####
  # 3 #
  ##### Step 2: Iterative MLE estimation for item parameters and ability pattern
  
  d.par <- d.att <- d.undefined <- 1
  ite <- 0
  
  while ((max(d.par) > conv.crit.par || d.att > conv.crit.att || d.undefined > 0) & ite < max.ite)
  {
    ite <- ite + 1
    #cat(paste(paste("Iteration:", ite), "\n"))
    nitem <- dim(Y)[2]
   
    # MLE of item parameters
    
    if (model == "DINA" || model == "DINO") 
    {
      temp <- ParMLE(Y, Q, alpha.est, model)
      par.est <- list(slip=temp$slip, guess=temp$guess, se.slip=temp$se.slip, se.guess=temp$se.guess)
      undefined.flag <- 1 - ((1 - is.infinite(par.est$slip)) * (1 - is.infinite(par.est$guess)) 
                             * (1 - is.nan(par.est$slip)) * (1 - is.nan(par.est$guess)))
    } else if (model == "NIDA") 
    {
      temp <- ParMLE(Y, Q, alpha.est, model)
      par.est <- list(slip=temp$slip, guess=temp$guess, se.slip=temp$se.slip, se.guess=temp$se.guess)
      undefined.flag <- rep(0, nitem)
    } else 
    {
      return(warning("Model specification is not valid."))
    }
       
    # MLE of alpha
    
    alpha.out.MLE <- AlphaMLE(Y, Q, par.est, model, undefined.flag)
    alpha.est <- alpha.out.MLE$alpha.est
    loglike.matrix <- alpha.out.MLE$loglike.matrix
    
    # Compute changes  
    
    if (ite > 1) 
    {
      d.par <- abs(unlist(par.est) - unlist(par.est.old))
      d.par <- d.par[(1 - is.infinite(d.par)) * (1 - is.na(d.par))]
      d.att <- mean(rowSums(abs(alpha.est - alpha.est.old)) != 0)
      d.undefined <- sum(abs(undefined.flag - undefined.flag.old))
    }
    
    par.est.old <- par.est
    alpha.est.old <- alpha.est
    undefined.flag.old <- undefined.flag
  }
  
  conv <- "Convergence criteria met."
  if (ite == max.ite) conv <- "Maximum iteration reached."

  output <- list(alpha.est=alpha.est, par.est=par.est, undefined.flag=undefined.flag, convergence=conv, n.ite=ite, loglike.matrix=loglike.matrix, NPloss.matrix=loss.matrix, alpha.est.NP=alpha.est.NP$alpha.est, NP.method=NP.method, model=model)
  class(output) <- "JMLE"
  return(output)
  
}