\name{EBelasticNet.Binomial}
\alias{EBelasticNet.Binomial}
\title{The Empirical Bayesian Elastic Net (EBEN) Algorithm for Binomial Model with Normal-Gamma(NG) Prior Distribution}
\description{Generalized linear regression, normal-Gamma (NG) hierarchical prior for regression coefficients \cr
			Multinomial model is not available in current package}
\usage{EBelasticNet.Binomial(BASIS, Target, lambda, alpha,Epis = "no",verbose = 0)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BASIS}{sample matrix; rows correspond to samples, columns correspond to features}
  \item{Target}{Class label of each individual, TAKES VALUES OF 0 OR 1}
  \item{lambda}{Hyperparameter controls degree of shrinkage; can be obtained via Cross Validation; lambda>0}
  \item{alpha}{Hyperparameter controls degree of shrinkage; can be obtained via Cross Validation; 0<alpha<1}
  \item{Epis}{"yes" or "no" for including two-way interactions}
  \item{verbose}{0 or 1; 1: display message; 0 no message}
}
\details{If Epis="yes", the program adds two-way interaction K*(K-1)/2 more columns to BASIS;
		for memory efficient, the function pass n_effect to C. n_effect > n_true effects, which is 
		a rough guess on how many variables will be selected by the function.
		by providing a relative 'small' n_effect, the function will not allocate
		a large trunk of memory during computation. If you run into situations where the the model ran out of memory and needed to select
		more variables than the guess, please feel free to contact with me regarding to modifying the memory allocations.\cr\cr
		By default, only the diagonal elements of the posterior covariance matrix are provided in the output for t-test. In cases where the whole covariance 
		matrix is interested, please also contact with me for technique support.}
\value{
	\item{weight}{the none-zero regression coefficients: \cr 
	col1,col2 are the indices of the bases(main if equal); \cr
	col3: coefficent value; \cr
	col4: posterior variance;\cr
	col5: t-value;\cr
	col6: p-value}
	\item{logLikelihood}{log likelihood from the final regression coefficients}
	\item{WaldScore}{Wald Score}
	\item{Intercept}{Intercept}
	\item{lambda}{the hyperparameter; same as input lambda}
	\item{alpha}{the hyperparameter; same as input alpha}	
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
output = EBelasticNet.Binomial(BASIS, y,lambda = 1,alpha = 0.5, Epis = "no")
}
