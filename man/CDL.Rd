\name{CDL}
\alias{CDL}

\title{Log-likelihood for cognitive diagnostic models}
\description{
  This function returns the log-likelihood of a particular examinee's responses to a set of cognitive diagnostic items. Currently supported cognitive diagnostic models include DINA model, DINO model, and NIDA model. This function is called by the \code{AlphaMLE} function and the \code{JMLE} function in the package.
}

\usage{
CDL(Y, Q, par, alpha, model = "DINA", undefined.flag)
}
\arguments{
  \item{Y}{A vector of binary examinee responses. 1=correct, 0=incorrect.
}
  \item{Q}{The Q-matrix of the test. Rows represent items and columns represent attributes. 1=attribute required by the item, 0=attribute not required by the item. 
}
  \item{par}{A list of parameters.                                            
          DINA & DINO --- \code{par$slip}: a vector of slipping parameters for each item;       
                   \code{par$guess}: a vector of guessing parameters for each item.  
          NIDA --- \code{par$slip}: a vector of slipping parameters for each attribute; 
                   \code{par$guess}: a vector of guessing parameters for each attribute.
}
  \item{alpha}{A vector of examinee ability profile. 1=examinee masters the attribute, 0=examinee does not master the attribute.
}
  \item{model}{Currently supports three models: \code{"DINA"}, \code{"DINO"}, and \code{"NIDA"}. The default is \code{"DINA"}.
}
  \item{undefined.flag}{A binary vector indicating whether the parameters of each item are undefined. 1=undefined, 0=defined.
}
}

\value{
  \item{loglike}{The log likelihood function value for the given data.}
}

\seealso{
  \code{\link{AlphaMLE}}, \code{\link{JMLE}}
}

\examples{
# Generate item and examinee profiles

nitem <- 4
Q <- rbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1), c(1, 1, 1))
alpha <- c(1, 0, 0)

# Generate DINA model-based response data

slip <- rep(0.1, nitem)
guess <- rep(0.1, nitem)
my.par <- list(slip=slip, guess=guess)

data <- NA
eta <- NA

for (i in 1:nitem) {
  eta[i] <- prod(alpha ^ Q[i, ])
  P <- (1 - slip[i]) ^ eta[i] * guess[i] ^ (1 - eta[i])
  u <- runif(1)
  data[i] <- as.numeric(u < P)
}

# Using the function to compute the log-likelihood of the given data

CDL(data, Q, my.par, alpha, model="DINA", undefined.flag=rep(0, nitem))
}
