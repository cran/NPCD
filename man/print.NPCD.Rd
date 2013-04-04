\name{print.NPCD}
\alias{print.AlphaNP}
\alias{print.AlphaMLE}
\alias{print.ParMLE}
\alias{print.JMLE}
\alias{print.Qrefine}

\title{Print outputs generated from the functions in the package.}

\description{
  This function prints outputs generated from the functions in this package, including \code{\link{AlphaNP}}, \code{\link{AlphaMLE}}, \code{\link{ParMLE}}, \code{\link{JMLE}}, and\code{\link{Qrefine}}.
}

\usage{
\method{print}{AlphaNP}(x, \dots)
\method{print}{AlphaMLE}(x, \dots)
\method{print}{ParMLE}(x, \dots)
\method{print}{JMLE}(x, \dots)
\method{print}{Qrefine}(x, \dots)
}

\arguments{
  \item{x}{The output from the function (The list of all outputs).}
  \item{...}{Other arguments.}
}

\value{
  \item{AlphaNP}{The estimated examinee attribute profiles.}
  \item{AlphaMLE}{The estimated examinee attribute profiles.}
  \item{ParMLE}{The estimated item parameters and standard errors.}
  \item{JMLE}{The estimated item parameters and examinee attribute profiles.}
  \item{Qrefine}{The modified Q-matrix and the modified entries}
}

\seealso{
\code{\link{AlphaNP}}, \code{\link{AlphaMLE}}, \code{\link{ParMLE}}, \code{\link{JMLE}}, \code{\link{Qrefine}}
}

\examples{
# See examples in AlphaNP, AlphaMLE, ParMLE, JMLE, and Qrefine.
}
