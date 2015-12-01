\name{EBglmnet-package}
\alias{EBglmnet-package}
\docType{package}
\title{Empirical Bayesian Lasso (EBlasso) and Elastic Net (EBEN) Methods for Generalized Linear Models}

\description{ Fast Empirical Bayesian Lasso (EBlasso) and Elastic Net (EBEN) are generalized linear regression methods for variable selections and effect estimations. 
Similar as \code{lasso} and \code{elastic net} implemented in the package \pkg{glmnet}, \pkg{EBglmnet} features 
the capabilities of handling \eqn{p>>n} data, where \code{p} is the number of variables and \code{n} is 
the number of samples in the regression model, and inferring a sparse solution such that irrelevant variables 
will have exactly zero value on their regression coefficients. Additionally, there are several unique features in \pkg{EBglmnet}: \cr\cr
1) Both \code{EBlasso} and \code{EBEN} can select more than \code{n} nonzero effects. \cr
2) EBglmnet also performs hypothesis testing for the significance of nonzero estimates. \cr
3) EBglmnet includes built-in functions for epistasis analysis.\cr\cr
			 There three sets of hierarchical prior distributions implemented in \pkg{EBglmnet}: \cr\cr
1) EBlasso-NE is a two-level prior with (Normal + Exponential) distributions for the regression coefficients.\cr
2) EBlasso-NEG is a three-level hierarchical prior with (Normal + Exponential + Gamma) distributions.\cr
3) EBEN implements a normal and generalized gamma hierarchical prior.\cr

While those sets of priors are all "peak zero and flat tails", \code{EBlasso-NE} assigns more probability mass to the tails, resulting in more nonzero estimates having large \eqn{p-values}. In contrast, \code{EBlasso-NEG} has a third level constraint on the \code{lasso} prior, which results in higher probability mass around zero, thus more sparse results in the final outcome. Meanwhile, \code{EBEN} encourages a grouping effect such that highly correlated variables can be selected as a group.
			 Similar as the relationship between \code{elastic net} and \code{lasso}, while two parameters \eqn{(\alpha, \lambda)} are required for \code{EBEN}, it is reduced to \code{EBlasso-NE} when parameter \eqn{\alpha = 1}. We recommend to use `EBlasso-NEG` when there are a large number of candidate effects such as whole-genome epistasis analysis and GWAS, use `EBlasso-NE` when the number of candidate effects is small such as in multiple QTL mapping, and use `EBEN` when there are groups of highly correlated variables such as linkage disequilibrium and co-regulated gene expressions.
 \cr \cr 
			Two models are available for both methods: General linear regression model and Logistic regression model. Other features in this package includes:\cr 
			* 1 *     Epistasis (two-way interactions) can be included for all models/priors \cr
			* 2 *     model implemented with memory efficient c code. \cr
			* 3 *     LAPACK/BLAS are used for most linear algebra computations. \cr \cr
			
Several simulation and real data analysis in the reference papers demonstrated that \pkg{EBglmnet} enjoys better performance than \code{lasso} and \code{elastic net} methods in terms of Power of Detection, 
			False Discover Rate, as well as Encouraging Grouping Effect when applicable. The computational speed is comparable with that of 
			glmnet package implemented with fast Fortran code.\cr
			
			Key Algorithms are described in the following paper: \cr
			1. EBlasso-NEG: (Cai, Huang et al. 2011), (Huang, Xu et al. 2013) \cr
			2. EBlasso-NE: (Huang, Xu et al. 2013) \cr
			3. group EBlasso: (Huang, Martin et al. 2014) \cr
			4. EBEN: (Huang, Xu et al. 2015)\cr
			5. Whole-genome QTL mapping: (Huang, Xu et al, 2014) \cr
}
\details{
\tabular{ll}{
Package: \tab EBglmnet\cr
Type: \tab Package\cr
Version: \tab 4.0\cr
Date: \tab 2015-11-30\cr
License: \tab gpl\cr
}

}
\author{
Anhui Huang, Dianting Liu \cr
Maintainer: Anhui Huang <a.huang1@umiami.edu>
}
\references{
	Huang, A., Xu, S., and Cai, X. (2015). Empirical Bayesian elastic net for multiple quantitative trait locus mapping. Heredity 114(1): 107-115.\cr\cr	
	Huang, A., E. Martin, et al. (2014). "Detecting genetic interactions in pathway-based genome-wide association studies." Genet Epidemiol 38(4): 300-309.\cr\cr
	Huang, A., S. Xu, et al. (2014). "Whole-genome quantitative trait locus mapping reveals major role of epistasis on yield of rice." PLoS ONE 9(1): e87330.\cr\cr
	Huang, A. (2014). "Sparse model learning for inferring genotype and phenotype associations." Ph.D Dissertation. University of Miami(1186).\cr \cr	
	Huang A, Xu S, Cai X. (2013). Empirical Bayesian LASSO-logistic regression for multiple binary trait locus mapping. BMC genetics  14(1):5. \cr\cr
	Cai, X., Huang, A., and Xu, S. (2011). Fast empirical Bayesian LASSO for multiple quantitative trait locus mapping. BMC Bioinformatics 12, 211.\cr\cr

 
}
\keyword{ package }
