\name{ParMLE}
\alias{ParMLE}

\title{Maximum likelihood estimation of item parameters for cognitive diagnostic models.}
\description{
  This function returns maximum likelihood estimates of item parameters for cognitive diagnostic models when examinee ability patterns are known. This function can either be used independently or called in the \code{\link{JMLE}} function. Currently supported cognitive diagnostic models include the DINA model, the DINO model, the NIDA model, the G-NIDA model, and the R-RUM model.
}

\usage{
ParMLE(Y, Q, alpha, model = c("DINA", "DINO", "NIDA", "GNIDA", "RRUM"))
}
\arguments{
  \item{Y}{A matrix of binary responses. Rows represent persons and columns represent items. 1=correct, 0=incorrect.
}
  \item{Q}{The Q-matrix of the test. Rows represent items and columns represent attributes. 1=attribute required by the item, 0=attribute not required by the item.
}
  \item{alpha}{Examinee attribute profiles. Rows represent persons and columns represent attributes. 1=examinee masters the attribute, 0=examinee does not master the attribute.
}
  \item{model}{Currently support five models: \code{"DINA"}, \code{"DINO"}, \code{"NIDA"}, \code{"GNIDA"}, and \code{"RRUM"}. The default is \code{"DINA"}.
}
}

\value{
For the DINA, DINO, and NIDA models:
\item{slip}{a vector of slip parameters.}
\item{guess}{a vector of guessing parameters.}
\item{se.slip}{a vector of the standard errors for slip parameters.}
\item{se.guess}{a vector of the standard errors for guessing parameters.}
For the G-NIDA model: 
\item{slip}{a matrix (# items by # attributes) of slip parameters.}
\item{guess}{a matrix (# items by # attributes) of guessing parameters.}
\item{se.slip}{a matrix (# items by # attributes) of the standard errors for slip parameters.}
\item{se.guess}{a matrix (# items by # attributes) of the standard errors for guessing parameters.}
For the R-RUM model:
\item{pi}{a vector of pi parameters for each item.}
\item{r}{a matrix (# items by # attributes) of r parameters.}
\item{se.pi}{a vector of the standard errors for pi parameters.}
\item{se.r}{a matrix (# items by # attributes) of the standard errors for r parameters.}
Note that for the G-NIDA model and the R-RUM model, the item parameter estimates and standard errors are not available for the entries where the Q-matrix is 0.

Additionally, for all models:
\item{model}{The chosen model.}
\item{Q}{The Q-matrix of the test.}
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
ItemFit(parMLE.result)
ModelFit(parMLE.result)
}
