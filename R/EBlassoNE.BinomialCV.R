EBlassoNE.BinomialCV <- function(BASIS,Target,nFolds,foldId,Epis=FALSE, verbose = 0)
{
nStep= 19
	cat("EBLASSO Logistic Model, NE prior,Epis: ",Epis, ";", nFolds, "fold cross-validation\n");
	N 					= nrow(BASIS);
	K 					= ncol(BASIS);
	#set.seed(proc.time())
	if (missing(foldId)) 
	{
		if(N%%nFolds!=0){
			foldId 			= sample(c(rep(1:nFolds,floor(N/nFolds)),1:(N%%nFolds)),N);
		}else{
			foldId 			= sample(rep(1:nFolds,floor(N/nFolds)),N);
		}
	}
	lambda_Max			= log(1.1);
	ynorm 				= sqrt(sum(Target*Target));
	ynormal 			= Target/ynorm;
	
	for(i_b in 1:K){
		basis 			= BASIS[,i_b];
		basis 			= basis/sqrt(sum(basis*basis));
		corBy 			= basis%*%ynormal;
		if(corBy>lambda_Max) lambda_Max = corBy;
	}	
	if(Epis){
		for(i_b in 1:(K-1)){
			for(i_bj in (i_b + 1):K){
				basis 	= BASIS[,i_b]*BASIS[,i_bj];
				basis 	= basis/sqrt(sum(basis*basis));
				corBy 	= basis%*%ynormal;
				if(corBy>lambda_Max) lambda_Max = corBy;
			}
		}		
	}
	lambda_Max 			= lambda_Max;
	lambda_Min 			= log(0.0001*lambda_Max);
	step 				= (log(lambda_Max) - lambda_Min)/nStep;
	Lambda 				= exp(seq(from = lambda_Min,to=log(lambda_Max),by=step))
	N_step 				= length(Lambda);

	step 				= 1;
	alpha 				= 1;
	nAlpha 				= 1;
		
	Likelihood 			= mat.or.vec((N_step*nAlpha),4);
	logL 				= mat.or.vec(nFolds,1);

	for (i_s in 1:N_step){			
			lambda 		= Lambda[i_s];

	if(verbose >=0) cat("\tTesting step", step, "\t\tlambda: ",lambda,"\t")
			for(i in 1:nFolds){
				index  			= which(foldId!=i);
				Basis.Train 	= BASIS[index,];
				Target.Train 	= Target[index];
				index  			= which(foldId == i);
				Basis.Test  	= BASIS[index,];
				Target.Test 	= Target[index];
				SimF2fEB 		<-EBelasticNet.Binomial(Basis.Train,Target.Train,lambda,alpha,Epis, verbose);
				M				= length(SimF2fEB$fit)/6;
				Betas 			<- matrix(SimF2fEB$fit,nrow= M,ncol =6, byrow= FALSE);
				Mu  			= Betas[,3];
				Mu0 			= SimF2fEB$Intercept[1];
				
				rm(list="SimF2fEB");
				ntest 			= nrow(Basis.Test);
				#M 		= nrow(Betas);
				basisTest 	= matrix(rep(0,ntest*M),ntest,M);
				for(i_basis in 1:M){
					loc1 = Betas[i_basis,1];
					loc2 = Betas[i_basis,2];
					if(loc1==loc2){ 	basisTest[,i_basis] =  Basis.Test[,loc1];}
					else{			basisTest[,i_basis] =  Basis.Test[,loc1]* Basis.Test[,loc2];}
				}
				temp 		= exp(Mu0 + basisTest%*%Mu);
				logL[i] 	= mean(Target.Test*log(temp/(1+temp)) + (1-Target.Test)*log(1/(1+temp)));
			}
			Likelihood[step,] = c(alpha, lambda,mean(logL),sd(logL));
			if(verbose >=0) cat("log Likelihood",mean(logL),"\n");
			step 			= step + 1;
	}
		colnames(Likelihood) = c("alpha","lambda","logLikelihood","standard error");	

	index 				= which.max(Likelihood[,3]);
	Res.lambda			= Likelihood[index,2];
	Res.alpha 			= Likelihood[index,1];
	result 				<- list(Likelihood,Res.lambda);
	names(result)		<-c("CrossValidation","optimal hyperparameter");
	return(result);
}
