\name{Data.Qrefine}
\alias{Data.Qrefine}
\docType{data}
\title{Example dataset used for the \code{Qrefine} function}
\description{This dataset is a \code{list} containing components that can be used for the \code{Qrefine} function. The response marix is generated using the DINA model.

\code{$mis.Q}: The mis-specified Q matrix. Numerical, 20 (item) by 3 (attribute). Binary values: 1=attribute required by the item, 0=attribute not required by the item.

\code{$true.Q}: The true Q-matrix. Numerical, 20 (item) by 3 (attribute). Binary values: 1=attribute required by the item, 0=attribute not required by the item.

\code{$response}: The examinee responses to all items, generated using the true Q-matrix. Numerical, 300 (examinee) by 20 (item). Binary values: 1=correct, 0=incorrect.

\code{$true.alpha}: The true examinee ability profiles. Numerical, 300 (examinee) by 3 (attribute). Binary values: 1=examinee masters the attribute, 0=examinee does not master the attribute.

\code{$true.par$slip}: The true slipping parameters for each item. Numerical, vector of 20.

\code{$true.par$guess}: The true guessing parameters for each item. Numerical, vector of 20.

}
\usage{data(Data.Qrefine)}
\format{The format is: List of 5 elements.}
\keyword{datasets}
