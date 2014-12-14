EBlassoNEG.Gaussian <-
function(BASIS,Target,a_gamma,b_gamma,Epis = "no",verbose = 0,group = 1){
	N 					= nrow(BASIS);
	K 					= ncol(BASIS);
	if (verbose>0) cat("EBLASSO Gaussian Model, NEG prior, N: ",N,",K: ",K,", Epis: ",Epis,"\n");
	if(Epis == "yes"){
		#N_effect 		= (K+1)*K/2;
		N_effect 		= 2*K;
		Beta 			= rep(0,N_effect *4);

		#dyn.load("fEBLinearFullFloat.so")

		output<-.C("fEBLinearEpisEff",
			BASIS 		= as.double(BASIS),
			Target 		= as.double(Target),
			a_gamma 	= as.double(a_gamma),
			b_gamma 	= as.double(b_gamma),
			Beta 		= as.double(Beta),
			WaldScore 	= as.double(0),
			Intercept 	= as.double(0),
			N 			= as.integer(N),
			K 			= as.integer(K),
			ver 		= as.integer(verbose),
			bMax 		= as.integer(N_effect),
			residual 	= as.double(0),			
			group 		= as.integer(group),
			PACKAGE 	="EBglmnet");
		#dyn.unload("fEBLinearFullFloat.so")
	}else {
		N_effect 		= K;
		Beta 			= rep(0,N_effect *4);
		#dyn.load("fEBLinearMainEff.so")

		output<-.C("fEBLinearMainEff",
			BASIS 		= as.double(BASIS),
			Target 		= as.double(Target),
			a_gamma 	= as.double(a_gamma),
			b_gamma 	= as.double(b_gamma),
			Beta 		= as.double(Beta),
			WaldScore 	= as.double(0),
			Intercept 	= as.double(0),
			N 			= as.integer(N),
			K 			= as.integer(K),
			ver 		= as.integer(verbose),
			residual 	= as.double(0),
			PACKAGE		="EBglmnet");
#		dyn.unload("fEBLinearMainEff.so")
	}	
	
	result 				= matrix(output$Beta,N_effect,4);
	ToKeep 				= which(result[,3]!=0);
	if(length(ToKeep)==0) { Blup = matrix(0,1,4)
	}else
	{
		nEff 	= length(ToKeep);
		Blup 		= matrix(result[ToKeep,],nEff,4);
	}
	if(Epis == "yes"){
		blupMain 		= Blup[Blup[,1] ==Blup[,2],];
		nMain 			= length(blupMain)/4;
		blupMain 		= matrix(blupMain,nMain,4);
		#
		blupEpis 		= Blup[Blup[,1] !=Blup[,2],];
		nEpis 			= length(blupEpis)/4;
		blupEpis 		= matrix(blupEpis,nEpis,4);
		
		order1 			= order(blupMain[,1]);
		order2 			= order(blupEpis[,1]);
		Blup 			= rbind(blupMain[order1,],blupEpis[order2,]);	
	}
	t 				= abs(Blup[,3])/(sqrt(Blup[,4])+ 1e-20);
	pvalue 			= 1- pt(t,df=(N-1));
	Blup 			= cbind(Blup,t,pvalue); 			#M x 6
	#col1: index1
	#col2: index2
	#col3: beta
	#col4: variance
	#col5: t-value
	#col6: p-value
	fEBresult 			<- list(Blup,output$WaldScore,output$Intercept,output$residual,a_gamma,b_gamma);
	rm(list= "output")	
	names(fEBresult)	<-c("weight","WaldScore","Intercept","residVar","a","b")
	return(fEBresult)
	
}
