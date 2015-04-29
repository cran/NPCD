\name{NPCD-package}

\alias{NPCD-package}
\docType{package}

\title{
  NPCD: The \R Package for Nonparametric Methods for Cognitive Diagnosis
}

\description{
  This package implements an array of nonparametric and parametric estimation methods for cognitive diagnostic models, including nonparametric classification of examinee attribute profiles, joint maximum likelihood estimation (JMLE) of examinee attribute profiles and item parameters, and nonparametric refinement of the Q-matrix, as well as conditional maximum likelihood estimation (CMLE) of examinee attribute profiles given item parameters and CMLE of item parameters given examinee attribute profiles.
}

\details{
  \tabular{ll}{
    Package: \tab \code{NPCD}\cr
    Type: \tab Package\cr
    Version: \tab 1.0-8\cr
    Date: \tab 2015-4-29\cr
    License: \tab \acronym{LGPL} (>= 2.1)\cr
    Depends: \tab R (>= 2.14.2), BB, R.oo
  }
  
  Cognitive diagnostic models (CDM) are a group of latent class models specialized for cognitive diagnosis in educational testing. Based on the responses to the items in a test, CDMs classify each examinee into mastery or non-mastery on each of a few attributes, producing an "attribute profile" for the examinee. For different tests, the attributes can represent different constructs (e.g., skills, traits, strategies). An attribute profile can provide helpful diagnostic information of the examinee regarding the targeted constructs. Some examples of CDMs include the deterministic inputs, noisy "AND" gate (DINA) model (Junker & Sijtsma, 2001), the Deterministic Input, Noisy Output "OR" gate (DINO) model (Templin & Henson, 2006), the NIDA model (Maris, 1999), the reduced reparametrized unified model (R-RUM) (Hartz, Roussos, Henson, & Templin, 2005), the log-linear CDM (Hen son, Templin, & Willse, 2009), the general diagnostic model (GDM) (Davier, 2010), and the generalized DINA model (G-DINA) (de la Torre, 2011).
  
Despite the different ways to model the relationship between examinees' attribute profiles and their answers to test items, most of the existing CDMs make use of a so-called Q-matrix (Tatsuoka, 1985), which specifies which attributes are required by each item. The Q-matrix and the examinee attribute profile together produce the ideal response pattern to the items by an examinee. However, the real response pattern may not be exactly the same with the ideal response pattern due to other factors. Therefore, CDMs incorporate stochastic elements to account for random deviations from the ideal response pattern. For an acceptable model fit, the stochastic terms should allow some departure from the ideal response pattern, but not so much that the model becomes doubtful. And in such nice conditions, despite the distinctive stochastic terms and assumptions of different CDMs, according to the likelihood functions, the ideal response pattern should still be the most likely one among all possible response patterns. Therefore, classification of examinee attribute profiles based on the ideal response patterns can be effective without fitting the parametric models. Chiu and Douglas (2013) proposed a few nonparametric classification methods based on this idea. This \R package NPCD was built to carry out a group of methods that are all based on the above-mentioned nonparametric classification idea, including (1) the nonparametric classification methods proposed in Chiu and Douglas (2013), (2) a joint maximum likelihood estimation (JMLE) approach that utilizes the nonparametric classification methods to give good initial classifications of examinee attribute profiles and then estimates model-based attribute profiles and item parameters simultaneously, and (3) a nonparametric approach to refine the Q-matrix based on the residual sum of squares (RSS) criterion (Chiu, 2013).

Joint maximum likelihood estimation (JMLE) is a parameteric approach to fit a CDM to response data and estimate both the item parameters and examinee attribute profiles. Two other common methods are the expectation maximization (EM) algorithm and the Markov chain Monte Carlo (MCMC) technique. Despite JMLE's attractive advantages of mathematical simplicity and computational rapidity especially for more complex models, so far JMLE has not been successfully implemented for CDMs due to the potential inconsistency resulted from arbitrary initial values. In this package, we try to overcome this limitation by taking the classification results obtained from the nonparametric classification methods (Chiu & Douglas, 2013) as the initial values. Simulation results have shown JMLE is efficient and accurate for the models incorporated in this package.

In most practices, the Q-matrix of a test is constructed from the judgments of content experts. The misspecification of the Q-matrix can impair the estimation of the model parameters as well as classification of examinee attribute profiles. This package implements the Q-matrix refinement method developed by Chiu (2013), which is also based on the aforementioned nonparametric classification methods (Chiu & Douglas, 2013). This Q-matrix refinement method corrects potential misspecified entries of the Q-matrix through comparisons of the residual sum of squares computed from the observed and the ideal item responses. Interested users can refer to Chiu (2013) for more details.

The aim of the NPCD package is to provide an easy-to-use tool to analyze cognitive diagnostic test data using the aforementioned nonparametric methods as well as a few model-based methods for some CDMs. Currently the nonparametric methods in the package support both conjunctive and disjunctive models, and the parametric methods in the package support the DINA model, the DINO model, the NIDA model, the G-NIDA model, and the R-RUM model. Note that the G-NIDA model and the R-RUM model are equivalent based on reparametrization; nevertheless, we think it may still provide some convenicence to the users to include both of them us an option. The package also contains demonstrative data generated from some of these models. Besides the numeric outputs, the package also provides diagnostic plots for users to look inside the algorithms. 
}

\author{
  Yi Zheng (Arizona State University) and Chia-Yi Chiu (Rutgers, the State University of New Jersey)
  
  Maintainer: Yi Zheng <yi.isabel.zheng@asu.edu>
}

\references{
  Chiu, C. Y. (2011). \emph{Flexible approaches to cognitive diagnosis: nonparametric methods and small sample techniques}. Invited session of cognitive diagnosis and item response theory at 2011 Joint Statistical Meeting.
  
  Chiu, C. Y., & Douglas, J. A. (2013). A nonparametric approach to cognitive diagnosis by proximity to ideal response patterns. \emph{Journal of Classification 30}(2), 225-250.
  
  Chiu, C. Y. (2013). Statistical Refinement of the Q-matrix in Cognitive Diagnosis. \emph{Applied Psychological Measurement, 37}(8), 598-618. 
  
  Davier, M. (2010). A General Diagnostic Model Applied to Language Testing Data. \emph{British Journal of Mathematical and Statistical Psychology, 61}(2), 287-307.
  
  de la Torre, J. (2011). The Generalized DINA Model Framework. \emph{Psychometrika, 76}(2), 179-199.
  
  Hartz, S., Roussos, L., Henson, R., & Templin, J. (2005). The Fusion Model for Skill Diagnosis: Blending Theory with Practicality. \emph{Unpublished Manuscript}.
  
  Henson, R.A., Templin, J.L., & Willse, J.T. (2009). Defining a Family of Cognitive Diagnosis Models Using Log-Linear Models with Latent Variables. \emph{Psychometrika, 74}(2), 191-210.
  
  Junker, B., Sijtsma, K. (2001). Cognitive Assessment Models with Few Assumptions, and Connections with Nonparametric Item Response Theory. \emph{Applied Psychological Measurement, 25}(3), 258-272.
  
  Maris, E. (1999). Estimating Multiple Classification Latent Class Models. \emph{Psychometrika, 64}(2), 187-212.
  
  Tatsuoka, K. K. (1985). A probabilistic model for diagnosing misconceptions by the pattern classification approach. \emph{Journal of Educational and Behavioral Statistics, 10}, 55-73.

  Templin, J.L,, Henson, R.A. (2006). Measurement of Psychological Disorders Using Cognitive Diagnosis Models. \emph{Psychological Methods, 11}(3), 287-305.
}

\keyword{package}