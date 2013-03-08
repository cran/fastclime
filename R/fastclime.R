#-------------------------------------------------------------------------------#
# Package: fastclime                                                            #
# fastclime(): Main Function                                                    #
# Authors: Haotian Pang, Han Liu and Robert Vanderbei                           #
# Emails: <hpang@princeton.edu>, <hanliu@princeton.edu> and <rvdb@princetonedu> #
# Date: Mar 13th 2013                                                           #
# Version: 0.9						                        #
#-------------------------------------------------------------------------------#
fastclime <- function(x, lambda.min.ratio = NULL)
{
  
  gcinfo(FALSE)
  if(is.null(lambda.min.ratio)) lambda.min.ratio=0.1
  
  cov.input<-1
  SigmaInput<-x

  if(!isSymmetric(x))
  {
     SigmaInput<-cor(x)
     cov.input<-0
  }

  d<-dim(SigmaInput)[2]
 

  if(d<=50){
    nlambda=d
   }
   else if(d<=1000){
   nlambda=50
   }
   else{
   nlambda=15
   }

  cat("Allocating memory \n")
  maxnlambda=0
  mu_input<-matrix(0,nlambda,d)
  ipath<-matrix(0,nlambda,d*d)
  iicov<-matrix(0,nlambda,d*d)  
  loglik<-rep(0,nlambda)
  ratio<-lambda.min.ratio


 cat("start recovering \n")  
     str=.C("parametric", as.double(SigmaInput), as.integer(d), as.double(mu_input), as.double(ratio),as.integer(nlambda), as.integer(ipath), as.integer(maxnlambda), as.double(iicov), PACKAGE="fastclime")

 cat("preparing precision and path matrix list \n")  
  sigmahat<-matrix(unlist(str[1]),d)   
  mu<-matrix(unlist(str[3]),nlambda,d)
 
  ipath<-matrix(unlist(str[6]),nlambda,d*d)
  maxnlambda<-unlist(str[7])+1
  iicov<-matrix(unlist(str[8]),nlambda,d*d)
  
  sparsity<-rep(0,maxnlambda)
  df<-matrix(0,d,maxnlambda)
  loglik<-rep(0,maxnlambda)
   mu<-mu[1:maxnlambda,]
  path<-list()
  icov<-list()
  for (i in 1:maxnlambda)
  {
     tmppath<-matrix(ipath[i,],d,d)
     tmppath<-Matrix(ceiling((tmppath+t(tmppath))/2),sparse=TRUE)     
     path[i]<-list(tmppath)

     tmpicov<-matrix(iicov[i,],d,d)
     tmpicov<-(tmpicov+t(tmpicov))/2
     icov[i]<-list(tmpicov)

     sparsity[i]<-sum(tmppath)/(d^2-d)
     df[,i]=rowSums(tmppath)
     
  }


  result<-list("data" = x, "cov.input" = cov.input, "sigmahat" = sigmahat, "nlambda" = maxnlambda, "lambda" = mu,"path" = path, "sparsity" = sparsity, "icov" = icov, "df" = df)

  rm(x,cov.input,sigmahat,maxnlambda,mu,path,sparsity,icov,df,tmppath,iicov,ipath,
    nlambda,ratio,lambda.min.ratio,mu_input,SigmaInput,d)
  gc()
  class(result) = "fastclime"
  cat("Done! \n")

  return(result)

}

print.fastclime = function(x, ...)
{	

	if(x$cov.input) cat("Input: The Covariance Matrix\n")
	if(!x$cov.input) cat("Input: The Data Matrix\n")
	
	cat("Path length:",x$nlambda,"\n")
	cat("Graph dimension:",ncol(x$data),"\n")
	cat("Sparsity level:",min(x$sparsity),"----->",max(x$sparsity),"\n")
}


plot.fastclime = function(x, ...){
        gcinfo(FALSE)
        s<-x$lambda[,1]
        poslambda<-s[s>0]
        npos<-length(poslambda)
	
	plot(x$lambda[1:npos,1], x$sparsity[1:npos], log = "x", xlab = "Regularization Parameter", ylab = "Sparsity Level", type = "l", main = "Sparsity vs. Regularization")
	

}
