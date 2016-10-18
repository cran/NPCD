###############################################################################
# ModelFit:                                                                   #
#                                                                             #
# Generate model fit statistics for various classes of outputs generated from #
# the functions in this package, including ParMLE and JMLE. Currently supports#
# -2LL, AIC, and BIC.                                                         #
#                                                                             #
###############################################################################

#install.packages("R.oo")
#library(R.oo)

# inputs: x: the output from other functions
# outputs: model fit statistics

################################################################################

ModelFit <- function(x) {

	# Extract necessary information
	
	Q <- x$Q
	Y <- x$Y

	nitem <- nrow(Q)
	natt <- ncol(Q)
	nperson <- nrow(Y)
	pattern <- AlphaPermute(natt)
	npattern <- nrow(pattern)

	# Different treatments for classes
	
	if (class(x) %in% c("JMLE")) {
		model <- x$model
		loglike <- x$loglike
		
		if (model %in% c("DINA", "DINO")) npar <- 2 * nitem + 2 ^ natt - 1
		if (model == "NIDA") npar <- 2 * natt + 2 ^ natt - 1
		if (model == "GNIDA") npar <- 2 * natt * nitem + 2 ^ natt - 1
		if (model == "RRUM") npar <- nitem * (natt + 1) + 2 ^ natt - 1
		
	} else if (class(x) == "ParMLE") {
		model <- x$model
		par <- x
		est.class <- rep(NA, nperson)
		for (j in 1:npattern) est.class[apply(x$alpha, 1, function(m) all(m == pattern[j, ]))] <- j		
		
		if (model %in% c("DINA", "DINO")) 
    	{
      		undefined.flag <- 1 - ((1 - is.infinite(par$slip)) * (1 - is.infinite(par$guess)) 
                             * (1 - is.nan(par$slip)) * (1 - is.nan(par$guess)))
        } else {
        	undefined.flag <- rep(0, nitem)
        }
		
		loglike <- 0
		for (i in 1:nperson)
  		{
    		loglike <- loglike + CDL(Y[i, ], Q, par, pattern[est.class[i], ], model, undefined.flag)  
    	} 
    	
    	if (model %in% c("DINA", "DINO")) npar <- 2 * nitem
		if (model == "NIDA") npar <- 2 * natt
		if (model == "GNIDA") npar <- 2 * natt * nitem
		if (model == "RRUM") npar <- nitem * (natt + 1)
		
	} else {
		cat("Model fit statistics not appropriate for this class of object.\n")
	}
		
	AIC <- -2 * loglike + 2 * npar
	BIC <- -2 * loglike + npar * log(nperson)

	out <- list(AIC=AIC, BIC=BIC)

	return(out)
}