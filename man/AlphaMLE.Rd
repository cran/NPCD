\name{AlphaMLE}
\alias{AlphaMLE}

\title{Maximum likelihood estimation of attribute profile}
\description{
  This function returns the model-based Maximum likelihood estimator(s) of the cognitive diagnostic attribute profile(s). Currently supported cognitive diagnostic models include the DINA, DINO, NIDA, GNIDA, and R-RUM models.
}

\usage{
AlphaMLE(Y, Q, par, model = c("DINA", "DINO", "NIDA", "GNIDA", "RRUM"),
undefined.flag = NULL)
}

\arguments{
  \item{Y}{A matrix of binary responses. Rows represent persons and columns represent items. 1=correct, 0=incorrect.
}
  \item{Q}{The Q-matrix of the test. Rows represent items and columns represent attributes. 1=attribute required by the item, 0=attribute not required by the item.
}
  \item{par}{A list of parameters.                                            
          DINA & DINO --- \code{par$slip}: a vector of slipping parameters for each item;       
                   \code{par$guess}: a vector of guessing parameters for each item.  
          NIDA --- \code{par$slip}: a vector of slipping parameters for each attribute; 
                   \code{par$guess}: a vector of guessing parameters for each attribute.
          GNIDA --- \code{par$slip}: a matrix (items by attributes) of slipping parameters; 
                   \code{par$guess}: a matrix (items by attributes) of guessing parameters.
          RRUM --- \code{par$pi}: a vector of pi parameters for each item;
                   \code{par$r}: a matrix (items by attributes) of r parameters.                   
}
  \item{model}{Currently supports five models: \code{"DINA"}, \code{"DINO"}, \code{"NIDA"}, \code{"GNIDA"}, and \code{"RRUM"}. The default is \code{"DINA"}.
}
  \item{undefined.flag}{A binary vector indicating whether the parameters of each item are undefined. 1=undefined, 0=defined.  Generally, this argument is only needed in \code{"JMLE"} for the DINA and DINO models, where the data may generate undefined item parameters for some items.
}
}

\value{
  \item{alpha.est}{A matrix of estimated attribute profiles for all examinees. Rows represent persons and columns represent attributes. 1=examinee masters the attribute, 0=examinee does not master the attribute.}
  \item{est.class}{The class number (row index in \code{pattern}) for each person's attribute profile. It can also be used for locating the log-likelihood value in \code{loglike.matrix} for the estimated attribute profile for each person.}
  \item{n.tie}{Number of ties in the log-likelihood among the candidate attribute profiles for each person. When we encounter ties, one of the tied attribute profiles is randomly chosen.}
  \item{class.tie}{The class numbers (row index in \code{pattern}) of the tied patterns for each person.}
  \item{pattern}{All possible attribute profiles in the search space.}
  \item{loglike.matrix}{The matrix of the log-likelihood function values. Rows represent candidate attribute profiles in the same order with the pattern matrix; columns represent different examinees.}
}

\seealso{
\code{\link{AlphaNP}}, \code{\link{JMLE}}, \code{\link{print.AlphaMLE}}, \code{\link{plot.AlphaMLE}}
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

# Using the function to estimate examinee attribute profile

alpha.est.MLE <- AlphaMLE(data, Q, my.par, model="DINA", undefined.flag=NULL)

nperson <- 1   # Choose an examinee to investigate
print(alpha.est.MLE)   # Print the estimated examinee attribute profiles
plot(alpha.est.MLE, nperson)   # Plot the sorted log-likelihood function 
#of different attribute profiles for this examinee
ItemFit(alpha.est.MLE)
}
