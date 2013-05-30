
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
# (4) NP.method: "Hamming", the plain hamming distance method;                #
#             "Weighted", the hamming distance weighted by inversed variance  #
#             "Penalized", the hamming distance weighted by inversed variance #
#                          and specified penalizing weights for guess and slip#
# Additional input for the "penalized" NP.method:                             #
# (5) wg = weight assigned to guess                                           #
# (6) ws = weight assigned to slip                                            #

# (7) conv.crit.par: the critical value for the change in item parameter      #
#                    values to determine convergence                          #      
# (8) conv.crit.att: the critical value for the percentage of attributes that #
#                    are changed to determine convergence                     #
#                                                                             #
# Outputs:                                                                    #
# (1) item parameters estimates and standard errors for each item             #
# (2) attribute profiles for each examinee                                    #
#                                                                             #
###############################################################################


JMLE <- function(Y, Q, model=c("DINA", "DINO", "NIDA"), NP.method=c("Weighted", "Hamming", "Penalized"), wg=1, ws=1, conv.crit.par=0.001, conv.crit.att=0.01, max.ite=100) {

  #####
  # 1 #
  ##### Check dimension consistency and convert data to the right formats 
    
  Y <- as.matrix(Y)
  Q <- as.matrix(Q)
  check <- NULL
  check <- CheckInput(Y, Q)  
  if (!is.null(check)) return(warning(check))
  
  model <- match.arg(model)
  NP.method <- match.arg(NP.method)
  
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
  
  while ((max(d.par) > conv.crit.par || d.undefined > 0 || (d.att > conv.crit.att & max(d.par) > 0)) & ite < max.ite)
  {
    ite <- ite + 1
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
    
    # Compute loglikelihood
    
    loglike <- 0
    for (i in 1:ncol(loglike.matrix)){
      loglike <- loglike + loglike.matrix[alpha.out.MLE$est.class[i],i]
    }
    
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
    
    cat(sprintf("Iteration: %d, loglike = %.5f, diff.par = %.5f, diff.att = %.5f, diff.undefined = %.5f", ite, loglike, max(d.par), d.att, d.undefined))
    cat("\n")
#   cat(format(unlist(par.est))); cat("\n")
#   cat(format(alpha.est[1:5,])); cat("\n")
#   cat(format(alpha.out.MLE$n.tie)); cat("\n")
  }
  
  conv <- "Convergence criteria met."
  if (ite == max.ite) conv <- "Maximum iteration reached."
  if (d.att > conv.crit.att) conv <- "Stationary item parameters reached. Estimated attribute profiles subject to random sampling."
    
  output <- list(alpha.est=alpha.est, par.est=par.est, n.tie=alpha.out.MLE$n.tie, undefined.flag=undefined.flag, 
                 loglike=loglike, convergence=conv, n.ite=ite, loglike.matrix=loglike.matrix, 
                 est.class=alpha.out.MLE$est.class, NP.loss.matrix=loss.matrix, 
                 NP.alpha.est=alpha.est.NP$alpha.est, NP.method=NP.method, 
                 NP.est.class=alpha.est.NP$est.class, pattern=alpha.est.NP$pattern, model=model)
  class(output) <- "JMLE"
  return(output)
  
}