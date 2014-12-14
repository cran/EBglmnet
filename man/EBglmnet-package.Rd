\name{EBglmnet-package}
\alias{EBglmnet-package}
\alias{EBglmnet}
\docType{package}
\title{Empirical Bayesian Lasso (EBlasso) and Elastic Net (EBEN) Methods for Generalized Linear Models}

\description{ Fast EBlasso and EBEN algorithms are generalized linear regression method for variable selections and effect estimations. 
			The algorithms are able to handle p>>n data with better performance than Lasso or elastic net in terms of Power of Detection, 
			False Discover Rate, as well as Encouraging Grouping Effect when applicable. The computational speed is comparable with that of 
			glmnet package implemented with fast Fortran code. \cr\cr
			EBlasso has implemented a three-level hierarchical prior: EBlasso-NEG (Normal + Exponential + Gamma) and a two-level prior EBlasso-NE (Normal + Exponential).
			EBEN implements a normal and generalized gamma hierarchical prior. Parameters controlling degree of shrinkage (i.e., degree of sparseness) need to be determined
			by cross validation (CV). For example, two parameters (alpha, lambda) are required for EBEN, which is equivalent with that of elastic net method. 
					When parameter alpha = 1, it is equivalent with EBlasso-NE (normal + exponential). \cr \cr 
			Two models are available for both methods: General linear regression model and Logistic regression model. Other features in this package includes:\cr 
			* 1 *     Epistasis (two-way interactions) can be included for all models/priors \cr
			* 2 *     model implemented with memory efficient c code. \cr
			* 3 *     LAPACK/BLAS are used for most linear algebra computations. \cr \cr
			Note: The EBglmnet will substitute previous R package of EBEN.\cr
			Key Algorithms are described in the following paper: \cr
			1. EBlasso-NEG: (Cai, Huang et al. 2011), (Huang, Xu et al. 2013) \cr
			2. EBlasso-NE: (Huang, Xu et al. 2013) \cr
			3. group EBlasso: (Huang, Martin et al. 2014) \cr
			4. EBEN: (Huang, Xu et al. 2014)
}
\details{
\tabular{ll}{
Package: \tab EBglmnet\cr
Type: \tab Package\cr
Version: \tab 3.5\cr
Date: \tab 2014-12-12\cr
License: \tab gpl\cr
}

}
\author{
Anhui Huang
Maintainer: Anhui Huang <a.huang1@umiami.edu>
}
\references{
key algorithms: \cr
	Cai, X., Huang, A., and Xu, S. (2011). Fast empirical Bayesian LASSO for multiple quantitative trait locus mapping. BMC Bioinformatics 12, 211.\cr
	Huang A, Xu S, Cai X. (2013). Empirical Bayesian LASSO-logistic regression for multiple binary trait locus mapping. BMC genetics  14(1):5. \cr
	Huang, A., Xu, S., and Cai, X. (2014). Empirical Bayesian elastic net for multiple quantitative trait locus mapping. Heredity 10.1038/hdy.2014.79 \cr \cr
Other publications: \cr
	Huang, A., E. Martin, et al. (2014). "Detecting genetic interactions in pathway-based genome-wide association studies." Genet Epidemiol 38(4): 300-309.\cr
	Huang, A., S. Xu, et al. (2014). "Whole-genome quantitative trait locus mapping reveals major role of epistasis on yield of rice." PLoS ONE 9(1): e87330.\cr
	Huang, A. (2014). "Sparse model learning for inferring genotype and phenotype associations." Ph.D Dissertation. University of Miami(1186).
}
\keyword{ package }
