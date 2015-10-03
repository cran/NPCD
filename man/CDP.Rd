\name{CDP}
\alias{CDP}

\title{Probability of correct response for cognitive diagnostic models}
\description{
  This function returns the model-predicted probability of correct response of one item for one person given the item parameters, Q vector, and alpha vector. Currently supported cognitive diagnostic models include the DINA model, DINO model, NIDA model, G-NIDA model, and R-RUM model. This function is called by the \code{ItemFit} function in the package.
}

\usage{
CDP(Q, par, alpha, model = c("DINA", "DINO", "NIDA", "GNIDA", "RRUM"))
}
\arguments{
  \item{Q}{The Q-vector of the item. Columns represent attributes. 1=attribute required by the item, 0=attribute not required by the item. 
}
  \item{par}{A list of parameters.                                            
          DINA & DINO --- \code{par$slip}: a scaler slip parameter for the item;       
                   \code{par$guess}: a scaler guessing parameter for the item.  
          NIDA --- \code{par$slip}: a vector of slip parameters for each attribute; 
                   \code{par$guess}: a vector of guessing parameters for each attribute.
          GNIDA --- \code{par$slip}: a vector of slip parameters for each attribute for the item; 
                   \code{par$guess}: a vector of guessing parameters for each attribute for the item.
          RRUM --- \code{par$pi}: a scaler pi parameter for the item;
                   \code{par$r}: a vector of r parameters for each attribute for the item.  
}
  \item{alpha}{A vector of examinee ability profile. 1=examinee masters the attribute, 0=examinee does not master the attribute.
}
  \item{model}{Currently supports five models: \code{"DINA"}, \code{"DINO"}, \code{"NIDA"}, \code{"GNIDA"}, and \code{"RRUM"}. The default is \code{"DINA"}.
}
}

\value{
  \item{P}{The probability of correct response for the item by the person.}
}

\examples{
# Generate item and examinee profiles

Q <- c(1, 0, 0)
alpha <- c(1, 0, 0)
slip <- 0.2
guess <- 0.1
my.par <- list(slip=slip, guess=guess)
CDP(Q, my.par, alpha, model="DINA")
}
