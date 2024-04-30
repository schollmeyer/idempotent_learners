### functions

optimize_twodimensional <- function(f, x0=(lb+ub)/2, lb=c(0,0),ub=c(10^8,10^8),tol=10^-10){

inner_function <- function(x){f(c(x,x0[2]))}

outer_function <- function(x){
  x0[2] <<- x
  result <- optimize(f=inner_function, interval =c(lb[1],ub[1]),tol=tol)
 # x0[1] <- result$objective
  return(result$objective)
 }
 result_2 <-optimize(f=outer_function,interval = c(lb[2],ub[2]),tol=tol)
 x0 <- c(0,result_2$minimum)
 
 result_1 <- optimize(f=inner_function, interval =c(lb[1],ub[1]),tol=tol)
return(list(solution=c(result_1$minimum,x0[2]),objective=result_1$objective))

}


classical_kernel_regression <- function(x,y,bandwidth,lambda){
  model <- gplm::kreg(x,y,grid=x,bandwidth=bandwidth)
  return(list(x=x,y=y,bandwidth=bandwidth,fitted=model$y))
}

classical_support_vector_regression <- function(x,y,bandwidth){
  model <- KRLS::krls(x,y,sigma=bandwidth,lambda=lambda,vcov=FALSE,derivative=FALSE,print.level=0)
  return(list(x=x,y=y,bandwidth=bandwidth,fitted=model$fitted))
}

idempotent_kernel_regression <- function(x,y,bandwidth,method="dirty_inversion",eps=10^-10,lambda=0,spiky=FALSE,rho=0,factor=10){
   kernel_matrix <- KRLS::gausskernel(X=x,sigma=bandwidth)
   if(spiky==TRUE){
	   kernel_matrix <- (1-rho)*kernel_matrix + rho*KRLS::gausskernel(X=x,sigma=bandwidth/factor)
	   }
   dirty_kernel_matrix <- kernel_matrix+diag(rep(lambda,length(x)))
   if(method=="quadratic_programming"){	
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
	solution <- solution$x   
   }
  if(method=="dirty_inversion"){
  q <- t(dirty_kernel_matrix)%*%dirty_kernel_matrix
  solution <- try(solution <- solve ( q, dirty_kernel_matrix%*%y),silent=TRUE)
  while(is.character(solution)){
	  q <- q + diag(rep(eps,length(x)))
	  eps <- eps*1.5
  #dirty_kernel_matrix <- dirty_kernel_matrix +diag(rep(eps,length(x)))
  #q <- t(dirty_kernel_matrix)%*%dirty_kernel_matrix
  solution <- try(solution <- solve ( q, dirty_kernel_matrix%*%y),silent=TRUE)

  }
	  hat_matrix <- kernel_matrix%*%solve(q,dirty_kernel_matrix)
  }
   # solution <<- solution
   # print(solution$status)
   return(list(x=x,y=y,fitted=kernel_matrix%*%solution,bandwidth=bandwidth,eps=eps,optimal_solution=solution,kernel_matrix=kernel_matrix,model=model,hat_matrix=hat_matrix))
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
	       
	       return(list(optimal_bandwidth=bandwidths[which.min(mses)],bandwidths=bandwidths,mses=mses))
   	
}
	  
bandwidth_optimization <- function(x,y_true,sd,learner=idempotent_kernel_regression,n_rep,method,opts=NULL,lb=c(0,0),ub=c(10000,100),bandwidth_only=FALSE){
	lambda <<- NULL
	bandwidth <<- NULL
	loss <<- NULL

	if(bandwidth_only){
	   f <- function(param){
	  result <-0
	  for(k in (1:n_rep)){
		  #y <- y_true+rnorm(length(y_true),sd=sd)
		result <- result+ mean((y_true-learner(x,y,bandwidth=param[1],lambda=0)$fitted)^2)
		  
	}
	 # lambda <<- c(lambda,param[2])
	 # bandwidth <<- c(bandwidth,(param[1]))
	 # loss <<- c(loss,result/n_rep)
  return(result/n_rep)}
		}
else{	
  f <- function(param){
	  result <-0
	  for(k in (1:n_rep)){
		  #y <- y_true+rnorm(length(y_true),sd=sd)
		result <- result+ mean((y_true-learner(x,y,bandwidth=param[1],lambda=param[2])$fitted)^2)
		  
	}
	# lambda <<- c(lambda,param[2])
	 # bandwidth <<- c(bandwidth,2^(param[1]))
	 # loss <<- c(loss,result/n_rep)
  return(result/n_rep)}
	}
	if(is.null(opts)){opts <- list(algorithm="NLOPT_GN_ESCH",maxeval=10^9)} #NLOPT_GN_DIRECT
        if(bandwidth_only){
	return(optimize(f=f,interval=c(0,50000),tol=tol))
		#nloptr(x0=c(100),eval_f=f,opts=opts,lb=0,ub=50000))
		}
	else{
	return(optimize_twodimensional(f=f,lb=c(0,0),ub=c(50000,10000),tol=tol))#nloptr(x0=c(100,10^-8),eval_f=f,opts=opts,lb=lb,ub=ub))
		}
	#return(optim(fn=f,par=c(10,lambdastart),method=method,control=control))
}

bandwidth_heuristic <- function(x,y,n_iter=100){
	bw <<- NULL
	rel_diff <<- NULL
	y_hat_1 <- mean(y)
	y_hat_2 <- y
	eps_1 <- y_hat_1 - y
	eps_2 <-y_hat_2 - y
	sd_1 <- sd(eps_1)
	sd_2 <- sd(eps_2)
	y_hat_start <- (y_hat_1 + y_hat_2)/2
	sd_start <- (sd_1+sd_2)/2
	sd <- sd_start
	y_start <- y_hat_start +rnorm(length(y),sd=sd)
	bandwidth <- bw_start
	y_true <- y#y_start
	y_new <- y_start +rnorm(length(x),sd=sd)
	for(k in (1:n_iter)){
		
		
		y_true <- idempotent_kernel_regression(x,y,bandwidth=bandwidth,lambda=0)$fitted
		sd <- sd(y-y_true)
		y_new <- y_true +rnorm(length(x),sd=sd)
		true_total_variation <- total_variation(x,y_true)
		estimated_model <- idempotent_kernel_regression(x,y_new,bandwidth=bandwidth,lambda=0)
		estimated_total_variation <- total_variation(x,estimated_model$fitted)
		if(true_total_variation < estimated_total_variation){bandwidth <- bandwidth * relaxation}
		else{bandwidth <- bandwidth / relaxation_2}
		
		print(true_total_variation)
		print(estimated_total_variation)
		bw <<- c(bw, bandwidth)
		rel_diff <<- c(rel_diff,abs(true_total_variation-estimated_total_variation)/true_total_variation)
		plot(rel_diff)
	}
	
	return(bandwidth)
}	

bwh <- function(x,y, bw_start){
	bw <<- NULL
	bandwidth <- bw_start
	while(TRUE){
		model <- idempotent_kernel_regression(x,y,bandwidth=bandwidth,lambda=0)
		sd <- sd(model$fitted-y)
		temp <- 0
		for(k in (1:10)){
		model_2 <- idempotent_kernel_regression(x,model$fitted+rnorm(length(y),sd=sd),bandwidth=bandwidth,lambda=0)
		temp <- temp + total_variation(x,model_2$fitted)
	        }
	        temp <- temp/10
		
		delta <- total_variation(x,model$fitted) - temp
		
		if(delta >0){
			bandwidth=bandwidth/relaxation
		}
		   else{bandwidth=bandwidth*relaxation}
		 bw <<- c(bw, delta)  
		plot(bw)	
		print(bandwidth)
	}}	

total_variation <- function(x,y){
	o <- order(x)
	x_ordered <- x[o]
	y_ordered <- y[o]
	result <- 0
	for(k in (1:(length(x)-1))){
		result <- result + abs(y_ordered[k+1]-y_ordered[k])
	}
	    
return(result)}
		
		
	
	
	
	
	
	
	


library(nloptr)

bandwidth_optimization(x,y,sd=400,n_rep=10)

### simulation scenarios
# scenario 1

n_sample <- 200
sd <- 5
x <- seq(1,n_sample,1)
y_true <- 5*sin(x/10)+.1*x
y <- y_true+ sd * rnorm(n_sample)
plot(x,y)
lines(x,y_true)
# idempotent model
ans <- bandwidth_optimization(x,y_true,sd=sd,n_rep=1,method="Nelder-Mead",control=list(maxit=100000000)

idempotent_tuning <- bandwidth_selection(x=x,y_true=y_true,sd=sd,learner=idempotent_kernel_regression,bandwidths=seq(100,5000,length.out=100),n_rep=10,eps=10^-8)
idempotent_model <- idempotent_kernel_regression(x=x,y=y,bandwidth=idempotent_tuning$optimal_bandwidth,eps=10^-8)
print(idempotent_tuning$mses)
lines(x,idempotent_model$fitted,col="blue")
mean((y_true-idempotent_model$fitted)^2)


# classical kernelregression
	 
classical_tuning <- bandwidth_selection(x=x,y_true=y_true,sd=sd,learner=classical_kernel_regression,bandwidths=seq(1,60,length.out=100),n_rep=10)
classical_model <- classical_kernel_regression(x=x,y=y,bandwidth=classical_tuning$optimal_bandwidth)
lines(x,classical_model$fitted,col="red")
mean((y_true-classical_model$fitted)^2)

# classical supportvectorregression
	 
classical_svr_tuning <- bandwidth_selection(x=x,y_true=y_true,sd=sd,learner=classical_support_vector_regression,bandwidths=seq(.1,5,length.out=100),n_rep=10)
classical_svr_model <- classical_support_vector_regression(x=x,y=y,bandwidth=classical_svr_tuning$optimal_bandwidth)
lines(x,classical_svr_model$fitted,col="purple")
mean((y_true-classical_svr_model$fitted)^2)


# scenario 2:

			      −2x+x sin(x)+ε

n_sample <- 200
sd <- 1
x <- seq(0,10,length.out=n_sample)
y_true <- -2*x+x*sin(2*x)
y <- y_true+  rnorm(n_sample,sd=sd)
plot(x,y)
lines(x,y_true)
