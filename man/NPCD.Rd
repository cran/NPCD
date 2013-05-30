\name{NPCD-package}

\alias{NPCD-package}
\docType{package}

\title{
  NPCD: The \R Package for Nonparametric Methods for Cognitive Diagnosis
}

\description{
  This package implements an array of nonparametric estimation methods for cognitive diagnostic models, including nonparametric estimation of examinee attribute profiles, joint maximum likelihood estimation (JMLE) of examinee attribute profiles and item parameters, and nonparametric refinement of the Q-matrix, as well as parametric algorithms such as maximum likelihood estimation (MLE) of examinee attribute profiles given item parameters and maximum likelihood estimation (MLE) of item parameters given examinee attribute profiles.
}

\details{
  \tabular{ll}{
    Package: \tab \code{NPCD}\cr
    Type: \tab Package\cr
    Version: \tab 1.0-2\cr
    Date: \tab 2013-04-04\cr
    License: \tab \acronym{LGPL} (>= 2.1)\cr
    Depends: \tab R (>= 2.14.2), BB, R.oo
  }
  
  Cognitive diagnostic models (CDM) are a group of latent class models specialized for cognitive diagnosis in educational measurement. Based on the responses to the items in a test, CDMs classify each examinee into mastery or non-mastery on each of a few attributes, producing an "attribute profile" for the examinee. For different tests, the attributes can represent different constructs (e.g., skills, traits, strategies). The attribute profile can provide helpful diagnostic information regarding the targeted constructs of the examinee. Some examples of CDMs include the deterministic inputs, noisy "AND" gate (DINA) (Junker & Sijtsma, 2001), the Deterministic Input, Noisy Output "OR" gate (DINO) (Templin & Henson, 2006), the NIDA model (Maris, 1999), the reduced reparametrized unified model (R-RUM) (Hartz, Roussos, Henson, & Templin, 2005), the log-linear CDM (Hen son, Templin, & Willse, 2009), the general diagnostic model (GDM) (Davier, 2010), and the generalized DINA model (G-DINA) (de la Torre, 2011).
  
Despite the different ways to model the relationship between examinees' attribute profiles and their answers to test items, most of the existing CDMs make use of a so-called Q-matrix, which specifies which attributes are required by each item. The Q-matrix and examinee attribute profile together produce the ideal response pattern to the items by the examinee. However, the real response pattern may not be exactly the same with the ideal response pattern due to other factors. Therefore, CDMs incorporate stochastic elements to account for random deviation from the ideal response pattern. For an acceptable model fit, the stochastic terms should allow some departure from the ideal response pattern, but not so much that the model becomes doubtful. And in such nice conditions, despite the distinctive stochastic terms and assumptions of different CDMs, according to the likelihood functions, the ideal response pattern should still be the most likely one among all possible response patterns. Therefore, classification of examinee attribute profiles based on the ideal response patterns can be effective without fitting the parametric models. Chiu and Douglas (In press) proposed a few nonparametric classification methods based on this idea. This \R package NPCD was built to carry out a group of methods that are all based on the above-mentioned nonparametric classification idea, including (1) the ones proposed in Chiu and Douglas (In press), (2) a joint maximum likelihood estimation (JMLE) approach that utilizes the nonparametric methods to give good initial classifications of examinee attribute profiles and then estimates model-based attribute profiles and item parameters simultaneously, and (3) a nonparametric approach to refine the Q-matrix based on the residual sum of squares (RSS) criterion (Chiu, in press).

The aim of the NPCD package is to provide an easy-to-use tool to analyze cognitive diagnostic test data using the aforementioned nonparametric methods as well as a few model-based methods for some CDMs. Currently the nonparametric methods in the package support both conjunctive and disjunctive models, and the parametric methods in the package support the DINA model, the DINO model, and the NIDA model. The package also contains demonstrative data generated from those models. Besides the numeric outputs, the package also provides diagnostic plots for users to look inside the algorithms. 
}

\author{
  Yi Zheng, Chia-Yi Chiu.
  
  Maintainer: Yi Zheng <yizheng1@illinois.edu>
}

\references{
  Chiu, C. Y. (2011). \emph{Flexible approaches to cognitive diagnosis: nonparametric methods and small sample techniques}. Invited session of cognitive diagnosis and item response theory at 2011 Joint Statistical Meeting.
  
  Chiu, C. Y., & Douglas, J. A. (In press). A model-free approach to cognitive diagnosis by proximity to ideal response patterns. \emph{Journal of Classification}.
  
  Chiu, C. Y. (In press). Statistical Refinement of the Q-matrix in Cognitive Diagnosis. \emph{Applied Psychological Measurement}. 
  
  Davier, M. (2010). A General Diagnostic Model Applied to Language Testing Data. \emph{British Journal of Mathematical and Statistical Psychology, 61}(2), 287-307.
  
  de la Torre, J. (2011). The Generalized DINA Model Framework. \emph{Psychometrika, 76}(2), 179-199.
  
  Hartz, S., Roussos, L., Henson, R., & Templin, J. (2005). The Fusion Model for Skill Diagnosis: Blending Theory with Practicality. \emph{Unpublished Manuscript}.
  
  Henson, R.A., Templin, J.L., & Willse, J.T. (2009). Defining a Family of Cognitive Diagnosis Models Using Log-Linear Models with Latent Variables. \emph{Psychometrika, 74}(2), 191-210.
  
  Junker, B., Sijtsma, K. (2001). Cognitive Assessment Models with Few Assumptions, and Connections with Nonparametric Item Response Theory. \emph{Applied Psychological Measurement, 25}(3), 258-272.
  
  Maris, E. (1999). Estimating Multiple Classification Latent Class Models. \emph{Psychometrika, 64}(2), 187-212.
  
  Templin, J.L,, Henson, R.A. (2006). Measurement of Psychological Disorders Using Cognitive Diagnosis Models. \emph{Psychological Methods, 11}(3), 287-305.
}

\keyword{package}