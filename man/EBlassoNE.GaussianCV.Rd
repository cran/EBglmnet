\name{EBlassoNE.GaussianCV}
\alias{EBlassoNE.GaussianCV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Cross Validation (CV) Function to Determine Hyperparameter of the Empirical Bayes Lasso Algorithm for Linear Model with Normal + Exponential Prior Distribution}

\description{Hyperparameter controls degree of shrinkage, and is obtained via Cross Validation (CV). This program calculates the maximum lambda
			that allows one non-zero basis; and performs a search down to 0.0001*lambda_max at even steps. (nStep = 20 by default) \cr\cr
			EBlasso-NE is a special case of EBEN with alpha = 1.  This function provides the optimal lambda (lambda_opt), and users should call the 
			function EBelasticNet.Gaussian with parameter (alpha = 1, lambda_opt) to obtain variable selection and effect estimation.
			}
\usage{
EBlassoNE.GaussianCV(BASIS, Target, nFolds, nStep, Epis = "no",foldId)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BASIS}{sample matrix; rows correspond to samples, columns correspond to features}
  \item{Target}{Response each individual}
  \item{nFolds}{number of n-fold cv}
   \item{nStep}{number of steps to be examined from lambda_max to lambda_min}
  \item{Epis}{"yes" or "no" for including two-way interactions}
  \item{foldId}{random assign samples to different folds(see example)}
}

\details{If Epis="yes", the program adds two-way interaction K*(K-1)/2 more columns to BASIS}
\value{
  \item{CrossValidation}{col1: hyperparameter; col2: loglikelihood mean; standard ERROR of nfold mean log likelihood}
  \item{Lambda_optimal}{the optimal hyperparameter as computed}

}

\references{
Huang, A., Xu, S., and Cai, X. (2013). Empirical Bayesian elastic net for multiple quantitative trait locus mapping. submitted.
}
\author{Anhui Huang; Dept of Electrical and Computer Engineering, Univ of Miami, Coral Gables, FL}



%% ~Make other sections like Warning with \section{Warning }{....} ~



\examples{
library(EBglmnet)
data(BASIS)
data(y)
#reduce sample size to speed up the running time
n = 50;
k = 100;
BASIS = BASIS[1:n,1:k];
y  = y[1:n];
CV = EBelasticNet.GaussianCV(BASIS, y, nFolds = 3,Epis = "no")
# to obtain variable selection and effect estimation
alpha = 1;
lambda = CV$Lambda_optimal;
fit=EBelasticNet.Gaussian(BASIS, y,lambda,alpha, Epis = "no")
}

