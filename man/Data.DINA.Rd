\name{Data.DINA}
\alias{Data.DINA}
\docType{data}
\title{Example dataset generated using DINA model}
\description{This dataset is a \code{list} containing components generated using DINA model, which can be used to test out the estimation functions in the package.
}
\usage{data(Data.DINA)}
\format{The format is: List of 4 elements.}

\value{
\item{$ Q}{The Q matrix. Numerical, 20 (item) by 3 (attribute). Binary values: 1=attribute required by the item, 0=attribute not required by the item.}
\item{$ response}{The examinee responses to all items. Numerical, 300 (examinee) by 20 (item). Binary values: 1=correct, 0=incorrect.}
\item{$ true.alpha}{The true examinee ability profiles. Numerical, 300 (examinee) by 3 (attribute). Binary values: 1=examinee masters the attribute, 0=examinee does not master the attribute.}
\item{$ true.par $ slip}{The true slipping parameters for each item. Numerical, vector of 20.}
\item{$ true.par $ guess}{The true guessing parameters for each item. Numerical, vector of 20.} 
}

\keyword{datasets}

