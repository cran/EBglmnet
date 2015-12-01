## ---- eval=FALSE---------------------------------------------------------
#  install.packages("EBglmnet", repos = "http://cran.us.r-project.org")

## ------------------------------------------------------------------------
rm(list = ls())
set.seed(1)
library(EBglmnet)

## ------------------------------------------------------------------------
varNames = colnames(state.x77);
varNames
y= state.x77[,"Life Exp"]
xNames = c("Population","Income","Illiteracy", "Murder","HS Grad","Frost","Area")
x = state.x77[,xNames]

## ------------------------------------------------------------------------
output = EBglmnet(x,y,hyperparameters = c(0.1, 0.1))

## ------------------------------------------------------------------------
glmfit = output$fit
variables = xNames[glmfit[,1,drop=FALSE]]
cbind(variables,as.data.frame(round(glmfit[,3:6,drop=FALSE],4)))

## ------------------------------------------------------------------------
cvfit = cv.EBglmnet(x, y)

## ------------------------------------------------------------------------
cvfit$CrossValidation

## ------------------------------------------------------------------------
cvfit$hyperparameters
cvfit$fit

## ------------------------------------------------------------------------
output$Intercept
output$residual

## ------------------------------------------------------------------------
yy = y>mean(y);
output = EBglmnet(x,yy,family="binomial", hyperparameters = c(0.1, 0.1))

## ------------------------------------------------------------------------
output = EBglmnet(x,yy,family="binomial", prior = "elastic net", hyperparameters = c(0.1, 0.1))

## ------------------------------------------------------------------------
output = EBglmnet(x,yy,family="binomial", prior = "elastic net", hyperparameters = c(0.1, 0.1),Epis = TRUE)
output$fit

