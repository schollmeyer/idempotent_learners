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
	model$Q <- t(kernel_matrix)%*%K +diag(rep(eps,length(x)))	
        solution <- try(gurobi(model,list(outputflag=0)), silent=TRUE)
   }
   # solution <<- solution
   # print(solution$status)
   return(kernel_matrix%*%solution$x)
   }
