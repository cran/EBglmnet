\name{cv.EBglmnet}
\alias{cv.EBglmnet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Cross Validation (CV) Function to Determine Hyperparameters of the EBglmnet Algorithms}
\description{The degree of shrinkage, or equivalently, the number of non-zero effects selected by EBglmnet are 
				controlled by the hyperparameters in the prior distribution, which can be obtained 
				via Cross Validation (CV). This function performs k-fold CV to obtain the optimal hyperparameters, and outputs the model fit results using the optimal parameters. Therefore, this function runs \code{EBglmnet} for (\code{k + 1}) times.
			}
\usage{
cv.EBglmnet(x, y, family=c("gaussian","binomial"),
		prior= c("lassoNEG","lasso","elastic net"), nfolds=5, 
		foldId, Epis = FALSE, group = FALSE, verbose = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{input matrix of dimension \code{n} x \code{p}; each row is an
  observation vector, and each column is a variable. When epistasis is considered, users do not need
  to create a giant matrix including both main and interaction terms. Instead, \code{x} should always be
  the matrix corresponding to the \code{p} main effects, and \code{cv.EBglmnet} will generate the interaction terms
  dynamically during running time.}
  \item{y}{response variable. Quantitative for \code{family="gaussian"}, and binary for 
  \code{family="binomial"}. For binary response variable, y can be a Boolean or numeric vector, or factor type
  array.}
  \item{family}{Model type taking values of "gaussian" (default) or "binomial". } 
  \item{prior}{Prior distribution to be used. Taking values of "lassoNEG"(default), "lasso", and "elastic net". 
  All priors will produce a sparse outcome of the regression coefficients, see Details for choosing priors. }  
  \item{nfolds}{number of n-fold CV. \code{nfolds} typically >=3. Although \code{nfolds}
    can be as large as the sample size (leave-one-out CV), it can be computationally intensively for large datasets. }
  \item{foldId}{an optional vector of values between 1 and \code{nfold}
    identifying which fold each observation is assigned to. If not supplied, each of the \code{n} samples will be 
	assigned to the \code{nfolds} randomly} 
  \item{Epis}{Boolean parameter for including two-way interactions. By default, \code{Epis = FALSE}. When
  \code{Epis = TRUE}, EBglmnet will take all pair-wise interaction effects into consideration. EBglmnet does not 
  create a giant matrix for all the \code{p(p+1)/2} effects. Instead, it dynamically allocates the memory for the 
  nonzero effects identified in the model, and read the corresponding variables from the original input matrix 
  \code{x}}
  \item{group}{Boolean parameter for \code{"Group EBlasso"} (currently only available for the \code{"lassoNEG"} 
  prior). This parameter is only valid when \code{Epis = TRUE}, 
	and is set to \code{FALSE} by default. When \code{Epis = TRUE} and \code{group = TRUE}, the hyperparameter
		controlling degree of shrinkage will be further scaled such that the scale hyperparameter for interaction terms
			is different with that of main effects by a factor of  \eqn{\sqrt{p(p-1)/2}}. 
			When \eqn{p} is large, eg., several thousands of genetic markers, 
			the total number of effects can easily be more than 10 millions, and \code{"Group EBlasso"} 
			helps to reduce the interference of spurious correlation and noise accumulation.}
  \item{verbose}{This parameter controls the level of message output from EBglment. It takes values from 0 to 5; larger verbose displays more messages. 0 is recommended for CV to avoid excessive outputs. Default value for \code{verbose} is minimum message output.}
}
\details{
	EBglmnet implements three set of hierarchical prior distributions for the regression parameters \eqn{\beta}: \cr\cr
	\bold{lasso prior}:
	\deqn{\beta_j  \sim N(0,\sigma_j^2),}
	\deqn{\sigma_j^2  \sim exp(\lambda), j = 1, \dots, p.}

	\bold{lasso-NEG prior}:
	\deqn{\beta_j  \sim N(0,\sigma_j^2),}
	\deqn{\sigma_j^2  \sim exp(\lambda),}
	\deqn{\lambda  \sim gamma(a,b), j = 1, \dots, p.}

	\bold{elastic net prior}:
	\deqn{\beta_j  \sim N[0,(\lambda_1 + \tilde{\sigma_j}^{-2})^{-2}],}
	\deqn{\tilde{\sigma_j}^{2}  \sim generalized-gamma(\lambda_1, \lambda_2), j = 1, \dots,p}

	The prior distributions are peak zero and flat tail probability distributions that assign a high 
	probability mass to zero and still allow heavy probability on the two tails, which reflect the 
	prior belief of a sparse solution exists: most of the variables will have no effects on the response 
	variable modeled, and only some of the variables will have non-zero effects contributing the variance in \code{y}. \cr

	The three priors all contain hyperparameters that control how heavy the tail probability is, 
	and different values of them will yield different number of non-zero effects retained in the model. 
	Appropriate selection of their values is required to obtain optimal results, and CV is the most 
	oftenly used method. For Gaussian model, CV determines the optimal hyperparameter values that yield 
	the minimum prediction error. In Binomial model, CV calculate the mean logLikelihood in each of 
	the left out fold, and choose the values that yield the maximum mean logLikelihood value of the k-folds.
	See \code{EBglmnet} for the support of hyperparameters in each prior distribution. \cr \cr

}		
		
		
\value{
	\item{CrossValidation}{This is a matrix containing all hyperparameters tested in the Cross Validation (CV). Each row corresponds to one candidate hyperparameter, and the columns are: \cr
	hyperparameter1, \cr
	hyperparameter2(if exists), \cr
	prediction metrics, \cr
	standard error in the n-fold CV. \cr
	Prediction metrics is the mean prediction error for Gaussian model: \cr
	\eqn{1/n||y-\hat{y}||^2}, where \eqn{\hat{y}} is the predicted outcome in the holdout fold using the parameters estimated by the (nfold -1) training set during one fold of CV. \cr 
		Prediction metrics is logLikelihood for the Binomial model: \cr
	\eqn{1/n\sum{[y_ilog(p_i) + (1-y_i)log(1-p_i)}]},  where \eqn{\hat{p_i}} is the predicted outcome in the holdout fold using the parameters estimated by the (nfold -1) training set during one fold of CV.
	}
	\item{optimal hyperparameter}{the optimal hyperparameter as computed}
	\item{fit}{the model fit using the optimal parameters as computed by CV. See \code{EBglmnet} for values in this item. }
	\item{WaldScore}{the Wald Score for the posterior distribution. It is computed as \eqn{\beta^T\Sigma^{-1}\beta}. 
	See (Huang A, 2014b) for using Wald Score to identify significant  effect set.}
	\item{Intercept}{the intercept in the linear regression model. This parameter is not shrunk as the regression coefficients}
	\item{residual variance}{the residual variance if the Gaussian family is assumed in the GLM}
	\item{logLikelihood}{the log Likelihood if the Binomial family is assumed in the GLM}
	\item{hyperparameters}{the optimal hyperparameter as computed}
	\item{family}{the GLM family specified in this function call}
	\item{prior}{the prior used in this function call}
	\item{call}{the call that produced this object}	
	\item{nobs}{number of observations}	
	\item{nfolds}{number of folds in CV}	
}

\references{	Cai, X., Huang, A., and Xu, S. (2011). Fast empirical Bayesian LASSO for multiple quantitative trait locus mapping. BMC Bioinformatics 12, 211.\cr\cr
	Huang A, Xu S, Cai X. (2013). Empirical Bayesian LASSO-logistic regression for multiple binary trait locus mapping. BMC genetics  14(1):5. \cr\cr
	Huang, A., Xu, S., and Cai, X. (2014a). Empirical Bayesian elastic net for multiple quantitative trait locus mapping. Heredity 10.1038/hdy.2014.79 \cr\cr 
uang, A., E. Martin, et al. (2014b). "Detecting genetic interactions in pathway-based genome-wide association studies." Genet Epidemiol 38(4): 300-309.	
	\cr	}
\author{Anhui Huang and Dianting Liu \cr Dept of Electrical and Computer Engineering, Univ of Miami, Coral Gables, FL}



\examples{

rm(list = ls())
library(EBglmnet)
#Use R built-in data set state.x77
y= state.x77[,"Life Exp"]
xNames = c("Population","Income","Illiteracy", "Murder","HS Grad","Frost","Area")
x = state.x77[,xNames]
#
#Gaussian Model
#lassoNEG prior as default
out = cv.EBglmnet(x,y)
out$fit
#lasso prior
out = cv.EBglmnet(x,y,prior= "lasso")
out$fit
#elastic net prior
out = cv.EBglmnet(x,y,prior= "elastic net")
out$fit
#
#Binomial Model
#create a binary response variable
yy = y>mean(y);
out = cv.EBglmnet(x,yy,family="binomial")
out$fit
#with epistatic effects
out = cv.EBglmnet(x,yy,family="binomial",prior= "elastic net",Epis =TRUE)
out$fit
}

