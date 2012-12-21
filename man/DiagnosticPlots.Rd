\name{DiagnosticPlots}
\alias{plot.AlphaNP}
\alias{plot.AlphaMLE}
\alias{plot.JMLE}
\alias{plot.Qrefine}

\title{Produce diagnostic plots}

\description{
  This function produces diagnostic plots of various outputs generated from the functions in this package, including \code{\link{AlphaNP}}, \code{\link{AlphaMLE}}, \code{\link{JMLE}}, and\code{\link{Qrefine}}.
}

\usage{
\method{plot}{AlphaNP}(x, nperson, cex.main, cex.legend, \dots)
\method{plot}{AlphaMLE}(x, nperson, cex.main, \dots)
\method{plot}{JMLE}(x, nperson, cex.main, \dots)
\method{plot}{Qrefine}(x, filename="Qrefine.plot.png", cex.main,
                 cex.lab, cex.axis, cex.legend, \dots)
}

\arguments{
  \item{x}{The output from the function (The list of all outputs).}
  \item{nperson}{The choice of examinee to be investigated.}
  \item{filename}{The filename of the plot, ending in ".png". Directory can be included; otherwise the plot will be saved in the working directory.}
  \item{cex.main}{A numerical value giving the amount by which title text should be magnified relative to the default. This starts as 1 when a device is opened.}
  \item{cex.axis}{The magnification to be used for axis annotation.}
  \item{cex.lab}{The magnification to be used for x and y labels.}
  \item{cex.legend}{The magnification to be used for the legend if there is one.}
  \item{...}{Other arguments.}
}

\value{
  \item{AlphaNP}{Bar plot of sorted loss function values (plain, weighted, or panelized Hamming distance) for each candidate attribute profile for the chosen examinee. The bar with the tilted shade is the estimated attribute profile for this examinee, which ideally should have the smallest loss function value.}
  \item{AlphaMLE}{Bar plot of sorted negative log-likelihood function values for each candidate attribute profile for the chosen examinee.The bar with the tilted shade is the estimated attribute profile for this examinee, which ideally should have the smallest value.}
  \item{JMLE}{A pair of barplots. Plot1: Bar plot of unsorted loss function values for each candidate attribute profile for the chosen examinee using the nonparametric estimation; Plot2: bar plot of unsorted negative log-likelihood function values for each candidate attribute profile for the chosen examinee in the last iteration of the JMLE estimation.The bar with the tilted shade is the estimated attribute profile for this examinee by each method, which ideally should have the smallest value.}
  \item{Qrefine}{A panel plot for all items with refined q-vectors. Each subplot shows the unsorted residual sum of squares (RSS) between the real response patterns and the ideal response patterns generated from each possible q-vector. The bar with the sparsely tilted shade is the initial q-vector for this item, and the bar with the densely tilted shade is the refined q-vector for this item, which ideally should have the smallest RSS value.}
}

\seealso{
\code{\link{AlphaNP}}, \code{\link{AlphaMLE}}, \code{\link{JMLE}}, \code{\link{Qrefine}}
}

\examples{
# See examples in AlphaNP, AlphaMLE, JMLE, and Qrefine.
}
