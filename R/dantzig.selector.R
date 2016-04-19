dantzig.selector <- function(lambdalist, BETA0, lambda){
	for (i in 1:length(lambdalist)) {
		if (lambdalist[i] < lambda) {
			break;
		}
	}
	beta0<-BETA0[,i]
	return(beta0)
}