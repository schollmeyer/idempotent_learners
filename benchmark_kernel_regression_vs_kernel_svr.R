idempotent_kernel_regression <- function(x,y,bandwidth,eps=10^-9){
   K <- gausskernel(X=x,sigma=bandwidth)
   model <- list(Q=t(K)%*%K + diag(rep(eps,length(x))),obj=-2*t(y)%*%K)
   model$A <- array(1,c(1,length(x)))
   model$sense <-"<="
   model$rhs <- sum(y)+1000000
   model$lb <- rep(-Inf,length(x))
   model$ub <- rep(Inf,length(x))
   solution <- try(solution <- gurobi(model,list(outputflag=0)), silent=TRUE)
   while(is.null(solution$x) | solution$status=="SUBOPTIMAL"){
        eps <- eps * 9/8
	# print(eps)
	model$Q <- t(K)%*%K +diag(rep(eps,length(x)))	
        solution <- try(gurobi(model,list(outputflag=0)), silent=TRUE)
   }
   # solution <<- solution
   # print(solution$status)
   return(K%*%solution$x)
   }
