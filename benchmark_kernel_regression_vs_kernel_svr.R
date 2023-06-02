### functions
classical_kernel_regression <- function(x,y,bandwidth){
  model <- kreg(x,y,grid=x,bandwidth=bandwidth)
  return(list(x=x,y=y,bandwidth=bandwidth,fitted=model$y))
}

idempotent_kernel_regression <- function(x,y,bandwidth,eps=10^-11){
   kernel_matrix <- KRLS::gausskernel(X=x,sigma=bandwidth)
   model <- list(Q=t(kernel_matrix)%*%kernel_matrix + diag(rep(eps,length(x))),obj=-2*t(y)%*%kernel_matrix)
   model$A <- array(1,c(1,length(x)))
   model$sense <- "<="
   model$rhs <- sum(y)+1000000
   model$lb <- rep(-Inf,length(x))
   model$ub <- rep(Inf,length(x))
   solution <- try(solution <- gurobi(model,list(outputflag=0)), silent=TRUE)
   while(is.null(solution$x) | solution$status=="SUBOPTIMAL"){
        eps <- eps * 9/8
	# print(eps)
	model$Q <- t(kernel_matrix)%*%kernel_matrix +diag(rep(eps,length(x)))	
        solution <- try(gurobi(model,list(outputflag=0)), silent=TRUE)
   }
   # solution <<- solution
   # print(solution$status)
   return(list(x=x,y=y,fitted=kernel_matrix%*%solution$x,bandwidth=bandwidth,eps=eps))
   }

	  
bandwidth_selection <- function(x,y_true,sd,learner,bandwidths,n_rep, ...){
	l <- 1
	mses <- rep(0,length(bandwidths))
	for(bandwidth in bandwidths){
	   
	   temp <- 0
	   for(k in (1:n_rep)){
		y <- y_true + rnorm(length(y),sd=sd)  
		 model <- learner(x=x,y=y,bandwidth=bandwidth,...)
		 temp <- temp + mean((y_true-model$fitted)^2)
	   }
	       mses[l] <- temp/n_rep
		
	       l <- l+1
	       
	}
	       
	       return(list(optimal_bandwidth=bandwidths[which.min(mses)],bandwidths=bandwidths,mese=mses))
   	
}
	  
	  

### simulation scenarios
# scenario 1

n_sample <- 200
sd <- 2
x <- seq(1,n_sample,1)
y_true <- 5*sin(x/10)+0.01*x
y <- y_true+ sd * rnorm(n_sample)
plot(x,y)
lines(x,y_true)
# idempotent model
optimal_bandwidth <- bandwidth_selection(x=x,y_true=y_true,sd=sd,learner=idempotent_kernel_regression,bandwidths=seq(2000,5000,length.out=40),n_rep=10,eps=10^-8)$optimal_bandwidth
idempotent_model <- idempotent_kernel_regression(x=x,y=y,bandwidth=optimal_bandwidth,eps=10^-8)
lines(x,idempotent_model$fitted,col="blue")
mean((y_true-idempotent_model$fitted)^2)	 
# classical kernelregression
	 
optimal_bandwidth <- bandwidth_selection(x=x,y_true=y_true,sd=sd,learner=classical_kernel_regression,bandwidths=seq(1,20,length.out=200),n_rep=10)$optimal_bandwidth
classical_model <- classical_kernel_regression(x=x,y=y,bandwidth=optimal_bandwidth)
lines(x,classical_model$fitted,col="red")
mean((y_true-classical_model$fitted)^2)





x <- seq(1,200,1)
y_true <- 5*sin(x^4/5000000)+0.01*x
y <- y_true+ 2 * rnorm(100)
