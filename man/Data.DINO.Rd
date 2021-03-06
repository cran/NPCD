\name{Data.DINO}
\alias{Data.DINO}
\docType{data}
\title{Example dataset generated using DINO model}
\description{This dataset is a \code{list} containing components generated using DINO model, which can be used to test out the estimation functions in the package.

\code{$Q}: The Q matrix. Numerical, 20 (item) by 3 (attribute). Binary values: 1=attribute required by the item, 0=attribute not required by the item.

\code{$response}: The examinee responses to all items. Numerical, 300 (examinee) by 20 (item). Binary values: 1=correct, 0=incorrect.

\code{$true.alpha}: The true examinee ability profiles. Numerical, 300 (examinee) by 3 (attribute). Binary values: 1=examinee masters the attribute, 0=examinee does not master the attribute.

\code{$true.par$slip}: The true slipping parameters for each item. Numerical, vector of 20.

\code{$true.par$guess}: The true guessing parameters for each item. Numerical, vector of 20.
}
\usage{data(Data.DINO)}
\format{The format is: List of 4 elements.}
\keyword{datasets}

