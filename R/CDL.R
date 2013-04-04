
###############################################################################
# CDL:                                                                        #
#                                                                             #
# Compute likelihood for cognitive diagnostic models.                         #
#                                                                             #
# Input:                                                                      #
# (1) Y: a vector of binary responses (1=correct, 0=incorrect).               #
# (2) Q: the Q-matrix of the test. Rows represent items and columns represent #
#        attributes.                                                          #
# (3) par: a list of parameters.                                              #
#          DINA&DINO --- par$slip: a vector of slip parameters for each item  #
#                   par$guess: a vector of guessing parameters for each item  #
#          NIDA --- par$slip: a vector of slip parameters for each attributes #
#                   par$guess: a vector of slip parameters for each attributes#
# (4) alpha: a vector of attribute profile.                                   #
# (5) model: "DINA", "DINO", "NIDA"                                           #
#                                                                             #
# Output:                                                                     #
# (1) loglike: the log likelihood for the given data.                         #
#                                                                             #
###############################################################################

CDL <- function(Y, Q, par, alpha, model="DINA", undefined.flag){
  
  nitem <- length(Y)
  natt <- dim(Q)[2]
  if (is.null(undefined.flag)) undefined.flag <- rep(0, nitem)
  
  par$slip[par$slip == 0] <- 0.001
  par$guess[par$guess == 0] <- 0.001
  par$slip[par$slip == 1] <- 0.999
  par$guess[par$guess == 1] <- 0.999
  
  if (model == "DINA") 
  {
    ita <- apply(matrix(alpha, byrow=T, nrow=nitem, ncol=natt) ^ Q, 1, FUN=prod)
    select <- (undefined.flag == 0)
    loglike.vec <- (Y[select] * ita[select] * log(1-par$slip[select]) 
                   + (1 - Y[select]) * ita[select] * log(par$slip[select]) 
                   + Y[select] * (1 - ita[select]) * log(par$guess[select]) 
                   + (1 - Y[select]) * (1 - ita[select]) * log(1 - par$guess[select]))     
    loglike <- sum(loglike.vec)   

  } else if (model == "DINO")
  {
    omega <- 1 - apply(matrix((1 - alpha), byrow=T, nrow=nitem, ncol=natt) ^ Q, 1, FUN=prod)
    select <- (undefined.flag == 0)
    loglike.vec <- (Y[select] * omega[select] * log(1-par$slip[select]) 
                    + (1 - Y[select]) * omega[select] * log(par$slip[select]) 
                    + Y[select] * (1 - omega[select]) * log(par$guess[select]) 
                    + (1 - Y[select]) * (1 - omega[select]) * log(1 - par$guess[select])) 
    loglike <- sum(loglike.vec)   
  
  } else if (model == "NIDA")
  {
    loglike <- 0
    for (j in 1:nitem)
    {
      P <- prod(((1 - par$slip) ^ alpha * par$guess ^ (1 - alpha)) ^ Q[j, ])
      loglike <- loglike + Y[j] * log(P) + (1 - Y[j]) * log(1 - P)
    }
        
  } else
  {
    return(warning("Model specification is not valid."))
  }
  
  return(loglike)
}
