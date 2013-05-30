
###############################################################################
# ParMLE                                                                      #
#                                                                             #
# Maximum likelihood estimation of item parameters for cognitive diagnostic   #
# models.                                                                     #
#                                                                             #
# Input:                                                                      #
# (1) Y: a matrix of binary responses (1=correct, 0=incorrect). Rows          #
#               represent persons and columns represent items.                #
# (2) Q: the Q-matrix of the test. Rows represent items and columns represent #
#        attributes.                                                          #
# (3) alpha: examinee attribute profiles. Rows represent persons and columns  #
#            represent attributes.                                            #
# (4) model: currently has three options, "DINA", "DINO", and "NIDA".         #
#                                                                             #
# Output:                                                                     #
# (1) slip: a vector or matrix of slip parameters                             #
# (2) guess: a vector or matrix of guessing parameters                        #
# (3) se.slip: a vector of or matrix standard error for slip parameters.      #
# (4) se.guess: a vector of or matrix standard error for guessing parameters. #
#                                                                             #
###############################################################################


ParMLE <- function(Y, Q, alpha, model=c("DINA", "DINO", "NIDA"))
{
  #####
  # 1 #
  ##### Check dimension consistency and convert data to the right formats 
  
  Y <- as.matrix(Y)
  Q <- as.matrix(Q)
  check <- NULL
  check <- CheckInput(Y, Q)  
  if (!is.null(check)) return(warning(check))
  
  model <- match.arg(model)
  
  #####
  # 2 #
  ##### Estimation
  
  nitem <- dim(Y)[2]
  nperson <- dim(Y)[1]
  natt <- dim(Q)[2]
  
  if (model == "DINA")
  {
    slip <- se.slip <- guess <- se.guess <- matrix(NA, nitem, 1)
    
    for (i in 1:nitem)
    {
      ita <- NULL
      for (j in 1:nperson)
      {
        ita[j] <- prod(alpha[j, ] ^ Q[i, ])
      }
      
      slip[i] <- sum((1 - Y[ , i]) * ita) / sum(ita)
      se.slip[i] <- sqrt(slip[i] * (1 - slip[i]) / sum(ita))
      guess[i] <- sum(Y[ , i] * (1 - ita)) / sum(1 - ita)   
      se.guess[i] <- sqrt(guess[i] * (1 - guess[i]) / sum(1 - ita))
    }        
        
  } else if (model == "DINO")
  {
    slip <- se.slip <- guess <- se.guess <- matrix(NA, nitem, 1)
    
    for (i in 1:nitem)
    {
      omega <- NULL
      for (j in 1:nperson)
      {
        omega[j] <- 1 - prod((1 - alpha[j, ]) ^ Q[i, ])
      }
            
      slip[i] <- sum((1 - Y[ , i]) * omega) / sum(omega)
      se.slip[i] <- sqrt(slip[i] * (1 - slip[i]) / sum(omega))
      guess[i] <- sum(Y[ , i] * (1 - omega)) / sum(1 - omega)   
      se.guess[i] <- sqrt(guess[i] * (1 - guess[i]) / sum(1 - omega))
    }        
        
  } else if (model == "NIDA")
  {
    slip <- se.slip <- guess <- se.guess <- matrix(NA, natt, 1)

    Like <- function(par, alpha, Q, Y)
    {
      l <- 0
      for (i in 1:nitem)
      {
        for (j in 1:nperson)
        {
          P <- prod((1 - par[1:natt]) ^ (alpha[j, ] * Q[i, ]) * par[(natt + 1):(2 * natt)] ^ ((1 - alpha[j, ]) * Q[i, ]))
          l <- l + Y[j, i] * log(P) + (1 - Y[j, i]) * log(1 - P)
        }        
      }
      return(l)
    }
    
    derLike <- function(par, alpha, Q, Y)
    {      
      dls <- dlg <- rep(NA, natt)
      
      for (k in 1:natt)
      {
        dls[k] <- 0
        dlg[k] <- 0
        
        for (i in 1:nitem)
        {
          for (j in 1:nperson)
          {
            prod.temp <- prod((1 - par[1:natt]) ^ (alpha[j, ] * Q[i, ]) * par[(natt + 1):(2 * natt)] ^ ((1 - alpha[j, ]) * Q[i, ]))
            dls[k] <- dls[k] + alpha[j, k] * Q[i, k] * (Y[j, i] - prod.temp) / (1 - prod.temp)
            dlg[k] <- dlg[k] + (1 - alpha[j, k]) * Q[i, k] * (Y[j, i] - prod.temp) / (1 - prod.temp)
          }
        }
      }     
      return(c(dls, dlg))
    }
    
    p0 <- rep(0.3, 2 * natt)
    #ans <- BBoptim(par=p0, fn=Like, gr=NULL, alpha=alpha, Q=Q, Y=Y, lower=0, upper=1, control=list(maximize=TRUE))
    ans <- BBsolve(par=p0, fn=derLike, alpha=alpha, Q=Q, Y=Y)
    #ans <- dfsane(par=p0, fn=derLike, tol=1.e-04, alpha=alpha, Q=Q, Y=Y)
    slip <- ans$par[1:natt]
    guess <- ans$par[(natt + 1):(2 * natt)]
    
    slip[slip < 0] <- 0; slip[slip > 1] <- 1
    guess[guess < 0] <- 0; guess[guess > 1] <- 1
    
    for (k in 1:natt)
    {
      se.slip[k] <- sqrt(slip[k] * (1 - slip[k]) / (sum(alpha[, k]) * nitem))
      se.guess[k] <- sqrt(guess[k] * (1 - guess[k]) / ((nperson - sum(alpha[,k])) * nitem))
    } 
          
  } else
  {
    return(warning("Model specification is not valid."))
  } 
  
  output <- list(slip=slip, guess=guess, se.slip=se.slip, se.guess=se.guess, model=model)
  class(output) <- "ParMLE"
  return(output)
  
}

