\name{ParMLE}
\alias{ParMLE}

\title{Maximum likelihood estimation of item parameters for cognitive diagnostic models.}
\description{
  This function returns maximum likelihood estimates of item parameters for cognitive diagnostic models when examinee ability patterns are known. This function can either be used independently or called in the \code{\link{JMLE}} function. Currently supported cognitive diagnostic models include the DINA model, the DINO model, and the NIDA model.
}

\usage{
ParMLE(Y, Q, alpha, model = "DINA")
}
\arguments{
  \item{Y}{A matrix of binary responses. Rows represent persons and columns represent items. 1=correct, 0=incorrect.
}
  \item{Q}{The Q-matrix of the test. Rows represent items and columns represent attributes. 1=attribute required by the item, 0=attribute not required by the item.
}
  \item{alpha}{Examinee attribute profiles. Rows represent persons and columns represent attributes. 1=examinee masters the attribute, 0=examinee does not master the attribute.
}
  \item{model}{Currently support three models: \code{"DINA"}, \code{"DINO"}, and \code{"NIDA"}. The default is \code{"DINA"}.
}
}

\value{
  \item{slip}{A vector of the estimates of the slipping parameters.}
  \item{guess}{A vector of the estimates of the guessing parameters.}
  \item{se.slip}{A vector of the standard error for the slipping parameters.}
  \item{se.guess}{A vector of the standard error for the guessing parameters.}
}

\seealso{
  \code{\link{JMLE}}, \code{\link{print.ParMLE}}  
}

\examples{
# Generate item and examinee profiles

natt <- 3
nitem <- 4
nperson <- 5
Q <- rbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1), c(1, 1, 1))
alpha <- rbind(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0), c(0, 0, 1), c(1, 1, 1))

# Generate DINA model-based response data

slip <- c(0.1, 0.15, 0.2, 0.25)
guess <- c(0.1, 0.15, 0.2, 0.25)
my.par <- list(slip=slip, guess=guess)

data <- matrix(NA, nperson, nitem)
eta <- matrix(NA, nperson, nitem)

for (i in 1:nperson) {
  for (j in 1:nitem) {
  eta[i, j] <- prod(alpha[i,] ^ Q[j, ])
  P <- (1 - slip[j]) ^ eta[i, j] * guess[j] ^ (1 - eta[i, j])
  u <- runif(1)
  data[i, j] <- as.numeric(u < P)
  }
}

# Using the function to estimate item parameters

parMLE.result <- ParMLE(data, Q, alpha, model="DINA")
print(parMLE.result)  # Print the estimated item parameters and standard errors

}
