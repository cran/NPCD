
###############################################################################
# AlphaNP:                                                                    #
#                                                                             #
# Estimate attribute profiles using the nonparametric approaches for the "AND #
# gate" or the "OR gate" cognitive diagnostic models.                         #
#                                                                             #
# Input:                                                                      #
# (1) Y: a matrix of binary responses (1=correct, 0=incorrect). Rows          #
#               represent persons and columns represent items.                #
# (2) Q: the Q-matrix of the test. Rows represent items and columns represent #
#        attributes.                                                          #
# (3) gate: "AND", the examinee needs to possess all related attributes       #
#                  in order to answer an item correctly                       #
#           "OR", the examinee needs to possess only one of the related       #
#                 attributes in order to answer an item correctly             #
# (4) method: "Hamming", the plain hamming distance method;                   #
#             "Weighted", the hamming distance weighted by inversed variance  #
#             "Panelized", the hamming distance weighted by inversed variance #
#                          and specified panelizing weights for guess and slip#
# Additional input for the "panelized" method:                                #
# (5) wg = weight assigned to guess in DINA model and DINO model              #
# (6) ws = weight assigned to slip in DINA model and DINO model               #
#                                                                             #
# Output:                                                                     #
# (1) alpha.est: estimated ability patterns                                   #
# (2) est.ideal: estimated ideal response to all the items for all the        #
#                 examinees                                                   #
# (3) n.ite: number of ties in the hamming distance. When we encounter ties,  # 
#            the best way is to radomly choose one. However, we want to       #
#           monitor the number of ties.                                       #
# (4) pattern: all attribute profiles in the search space.                    #
# (5) loss.matrix: The values for the loss function. Rows represent candidate #
#                  attribute profiles in the same order with the pattern      #
#                 matrix; Columns represent different examinees.              #
# For now, we report the estimated classes, estimated ideal response patterns #
# (mostly for the Q-refinement algorithm) and number of ties. We may want to  #
# report more later, such as the final RSS.                                   # 
#                                                                             #
###############################################################################

AlphaNP <- function(Y, Q, gate="AND", method="Hamming", wg=1, ws=1) {
  
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
  M <- 2 ^ natt
  
  # all possible alpha vectors
  
  pattern <- AlphaPermute(natt)
  
  # ideal responses for each possible alpha vector
  
  Ideal <- matrix(NA, M, nitem)
  for (m in 1:M) 
  {
    for (j in 1:nitem) 
    {
      if (gate == "AND")
      {
        u <- prod(pattern[m, ] ^ Q[j, ])
      } else if (gate == "OR")
      {
        u <- 1 - prod((1 - pattern[m, ]) ^ Q[j, ])
      } else 
      {
        return(warning("Gate specification not valid."))
      }
      Ideal[m, j] <- u
    }
  }
  
  if (method == "Hamming")
  {
    weight <- rep(1, nitem)
    ws <- wg <- 1
  } else if (method == "Weighted") 
  {
    p.bar <- apply(Y, 2, mean)
    weight <- 1 / (p.bar * (1 - p.bar))
    ws <- wg <- 1
  } else if (method == "Panelized") 
  {
    p.bar <- apply(Y, 2, mean)
    weight <- 1 / (p.bar * (1 - p.bar))    
    if (ws == wg) warning("Panelzing weights for guess and slip are the same --> equivalent with the \"Weighted\" method.")
  } else 
  {
    return(warning("Method specification not valid.")) 
  }
  
  loss.matrix <- matrix(NA, nrow=M, ncol=nperson)
  est.class <- NULL
  est.pattern <- NULL
  n.tie <-0
  
  for (i in 1:nperson)
  {
    Y.matrix <- matrix(rep(Y[i,], M), M, nitem, byrow=TRUE)
    loss <- apply(matrix(rep(weight, M), M, nitem, byrow=TRUE) * (wg * abs(Y.matrix - Ideal) * Y.matrix + ws * abs(Y.matrix - Ideal) * (1 - Y.matrix)), 1, sum)
    loss.matrix[,i] <- loss
    
    min.loss <- which(loss == min(loss))

    if (length(min.loss) != 1) 
    {
      n.tie <- n.tie + 1
      min.loss <- sample(min.loss, 1, prob=rep(1 / length(min.loss), length(min.loss))) 
    }
    
    est.class <- c(est.class, min.loss)
  }
  
  est.pattern <- pattern[est.class,]
  est.ideal <- Ideal[est.class,]
  output <- list(alpha.est=est.pattern, est.ideal=est.ideal, est.class=est.class, n.tie=n.tie, pattern=pattern, loss.matrix=loss.matrix, method=method)
  class(output) <- "AlphaNP"
  return(output)
}

