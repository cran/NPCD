###############################################################################
# print.NPCD:                                                                 #
#                                                                             #
# Define functions to print the various classes of outputs generated from the #
# functions in this package, including AlphaNP, AlphaMLE, JMLE, and Qrefine.  #                                                          #
#                                                                             #
###############################################################################

#install.packages("R.oo")
#library(R.oo)

################################################################################

# Print outputs for Nonparametric Alpha estimation
# inputs: x: the output from AlphaNP
# outputs: print the estimated examinee attribute profiles

setMethodS3("print", "AlphaNP", function(x, ...)  
{
 out <- as.matrix(x$alpha.est)
 rowname.tmp <- colname.tmp <- NULL
 for (i in 1:nrow(out)) rowname.tmp <- c(rowname.tmp, paste("Examinee", i))
 for (i in 1:ncol(out)) colname.tmp <- c(colname.tmp, paste("Attribute", i))
 rownames(out) <- rowname.tmp
 colnames(out) <- colname.tmp
 cat("The estimated examinee attribute profiles\n")
 cat(paste(paste("Method:", x$method), "\n"))
 print(out)
}
)

################################################################################

# Print outputs for MLE Alpha estimation
# inputs: x: the output from AlphaMLE
# outputs: print the estimated examinee attribute profiles

setMethodS3("print", "AlphaMLE", function(x, ...)  
{
  out <- as.matrix(x$alpha.est)
  rowname.tmp <- colname.tmp <- NULL
  for (i in 1:nrow(out)) rowname.tmp <- c(rowname.tmp, paste("Examinee", i))
  for (i in 1:ncol(out)) colname.tmp <- c(colname.tmp, paste("Attribute", i))
  rownames(out) <- rowname.tmp
  colnames(out) <- colname.tmp
  cat("The estimated examinee attribute profiles\n")
  cat("Method: conditional MLE\n")
  print(out)
}
)

################################################################################

# Print outputs for conditional MLE estimation of the item parameters
# inputs: x: the output from ParMLE
# outputs: print the estimated examinee attribute profiles

setMethodS3("print", "ParMLE", function(x, ...)  
{
  out <- as.matrix(cbind(x$slip, x$se.slip, x$guess, x$se.guess))
  rowname.tmp <- NULL
  if (x$model %in% c('DINA','DINO'))
  {
    for (i in 1:nrow(out)) rowname.tmp <- c(rowname.tmp, paste("Item", i))
  } else if (x$model == 'NIDA')
  {
    for (i in 1:nrow(out)) rowname.tmp <- c(rowname.tmp, paste("Attribute", i))
  }
  rownames(out) <- rowname.tmp
  colnames(out) <- c("slip","SE.slip","guess","SE.guess")
  cat("The estimated item parameters\n")
  cat(paste(paste("Model:", x$model), "\n"))
  cat("Method: conditional MLE\n")
  print(out)
  
}
)

################################################################################

# Print outputs for JMLE estimation of the examinee attribute profiles and item parameters
# inputs: x: the output from JMLE
# outputs: print the estimated examinee attribute profiles and the item parameters

setMethodS3("print", "JMLE", function(x, ...)  
{
  # Item parameters
  
  out <- cbind(x$par.est$slip, x$par.est$se.slip, x$par.est$guess, x$par.est$se.guess)
  rowname.tmp <- NULL
  if (x$model %in% c('DINA','DINO'))
  {
    for (i in 1:nrow(out)) rowname.tmp <- c(rowname.tmp, paste("Item", i))
  } else if (x$model == 'NIDA')
  {
    for (i in 1:nrow(out)) rowname.tmp <- c(rowname.tmp, paste("Attribute", i))
  }
  rownames(out) <- rowname.tmp
  colnames(out) <- c("slip","SE.slip","guess","SE.guess")
  cat(paste(paste("Model:", x$model), "\n"))
  cat("Method: JMLE\n")
  cat(paste(paste("Convergence:", x$conv), "\n"))
  cat(paste(paste("Number of iterations:", x$n.ite), "\n"))
  cat("The estimated item parameters\n")  
  print(out)
  
  # Examinee attribute profiles
  
  out <- as.matrix(x$alpha.est)
  rowname.tmp <- colname.tmp <- NULL
  for (i in 1:nrow(out)) rowname.tmp <- c(rowname.tmp, paste("Examinee", i))
  for (i in 1:ncol(out)) colname.tmp <- c(colname.tmp, paste("Attribute", i))
  rownames(out) <- rowname.tmp
  colnames(out) <- colname.tmp
  cat("\nThe estimated examinee attribute profiles\n")
  print(out)
}
)

################################################################################

# Print outputs for Q refinement
# inputs: x: the output from Qrefine
# outputs: print the refined Q-matrix and the modified entries

setMethodS3("print", "Qrefine", function(x, ...)  
{
  out <- as.matrix(x$modified.Q)
  rowname.tmp <- colname.tmp <- NULL
  for (i in 1:nrow(out)) rowname.tmp <- c(rowname.tmp, paste("Item", i))
  for (i in 1:ncol(out)) colname.tmp <- c(colname.tmp, paste("Attribute", i))
  rownames(out) <- rowname.tmp
  colnames(out) <- colname.tmp
  cat("The modified Q-matrix\n")
  print(out)
  
  cat("\nThe modified entries\n")
  print(x$modified.entries) 
}
)

