
###############################################################################
# AlphaMLE                                                                    #
#                                                                             #
# Maximum likelihood estimation of attribute profiles for cogitive diagnostic #
# models.                                                                     #
#                                                                             #
# Input:                                                                      #
# (1) Y: a matrix of binary responses (1=correct, 0=incorrect). Rows          #
#               represent persons and columns represent items.                #
# (2) Q: the Q-matrix of the test. Rows represent items and columns represent #
#        attributes.                                                          #
# (3) par: a list of parameters.                                              #
#          DINA & DINO --- par$slip: a vector of slip parameters for each item#
#                   par$guess: a vector of guessing parameters for each item  #
#          NIDA --- par$slip: a vector of slip parameters for each attributes #
#                   par$guess: a vector of slip parameters for each attributes#
# (4) model: "DINA", "DINO", "NIDA"                                           #
# (5) NP.method: "Hamming", "Weighted", "Panelized"                           #
# (6) undefined.flag: a binary vector indicating whether the parameters of    #
#                     each item are undefined (1=undefined, 0=defined).       #
#                                                                             #
# Output:                                                                     #
# (1) alpha.est: a matrix of estimated attribute profiles for all examinees   #
# (2) pattern: all attribute profiles in the search space.                    #
# (3) loglike.matrix: The values for the loss function. Rows represent        #
#                     candidate attribute profiles in the same order with the #
#                     pattern matrix; Columns represent different examinees.  #
#                                                                             #
###############################################################################

AlphaMLE <- function(Y, Q, par, model="DINA", undefined.flag=NULL) {
  
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
  ##### Estimation
  
  nperson <- dim(Y)[1]
  nitem <- dim(Q)[1]
  natt <- dim(Q)[2]
  pattern <-AlphaPermute(natt)
  loglike.matrix <- matrix(NA, dim(pattern)[1], nperson)  
  alpha.est <- matrix(NA, nperson, natt)
  
  for (i in 1:nperson)
  {
    loglike <- NULL
    
    for (j in 1:nrow(pattern))
    {
      loglike[j] <- CDL(Y[i, ], Q, par, pattern[j, ], model, undefined.flag)  
    }
    
    loglike.matrix[, i] <- loglike
    
    if (length(which(loglike == max(loglike))) == 1){
      alpha.est[i, ] <- pattern[which(loglike == max(loglike)), ]
    } else {
      alpha.est[i, ] <- pattern[sample(which(loglike == max(loglike)), 1), ]
    }     
  }  
  
  output <- list(alpha.est=alpha.est, pattern=pattern, loglike.matrix=loglike.matrix)
  class(output) <- "AlphaMLE"
  return(output)
}
