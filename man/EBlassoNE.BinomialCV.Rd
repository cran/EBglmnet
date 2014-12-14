\name{EBlassoNE.BinomialCV}
\alias{EBlassoNE.BinomialCV}

\title{ Cross Validation (CV) Function to Determine Hyperparameter of the Empirical Bayes Lasso Algorithm for Binomial Model with Normal + Exponential Prior Distribution}
\description{Hyperparameter controls degree of shrinkage, and is obtained via Cross Validation (CV). This program calculates the maximum lambda
			that allows one non-zero basis; and performs a search down to 0.0001*lambda_max at even steps. (nStep = 20 by default) \cr\cr
			EBlasso-NE is a special case of EBEN with alpha = 1.  This function provides the optimal lambda (lambda_opt), and users should call the 
			function EBelasticNet.Binomial with parameter (alpha = 1, lambda_opt) to obtain variable selection and effect estimation.
			}
\usage{
EBlassoNE.BinomialCV(BASIS, Target, nFolds,nStep =20, Epis = "no")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BASIS}{sample matrix; rows correspond to samples, columns correspond to features}
  \item{Target}{Class label of each individual, TAKES VALUES OF 0 OR 1}
  \item{nFolds}{number of n-fold cv}
  \item{nStep}{number of steps to be examined from lambda_max to lambda_min}
  \item{Epis}{"yes" or "no" for including two-way interactions}
}
\details{If Epis="yes", the program adds two-way interaction K*(K-1)/2 more columns to BASIS}
\value{
  \item{CrossValidation}{col1: hyperparameter; col2: loglikelihood mean; standard ERROR of nfold mean log likelihood}
  \item{Lambda_optimal}{the optimal hyperparameter as computed}
}
\references{Huang A, Xu S, Cai X: Empirical Bayesian LASSO-logistic regression for multiple binary trait locus mapping. BMC genetics 2013, 14(1):5.}
\author{Anhui Huang; Dept of Electrical and Computer Engineering, Univ of Miami, Coral Gables, FL}



%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
library(EBglmnet)
data(BASISbinomial)
data(yBinomial)
#reduce sample size to speed up the running time
n = 50;
k = 100;
BASIS = BASISbinomial[1:n,1:k];
y  = yBinomial[1:n];
CV = EBlassoNE.BinomialCV(BASIS, y, nFolds = 3,nStep = 10,Epis = "no")
# to obtain variable selection and effect estimation
alpha = 1;
lambda = CV$Lambda_optimal;
fit = EBelasticNet.Binomial(BASIS, y,lambda,alpha, Epis = "no")
}

