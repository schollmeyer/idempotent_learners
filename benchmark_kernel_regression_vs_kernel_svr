library(gplm)
library(krls)

n_rep <- 1000
b<- rep(0,20)
a <- rep(0,n_rep)
for(bw1 in (1:20)){

for(k in (1:n_rep)){

x <- seq(1,200,1)
y_true <- 5*sin(x/10)
y <- y_true+ 4 * rnorm(100)

M <- kreg(x,y,bandwidth=bw1,grid=x)

a[k] <- mean((M$y-y_true)^2)
}
b[bw1] <- mean(a)
print(bw1)}

which.min(b)

# 14
#####




bw1 <- 14
bw2 <- 500

x <- seq(1,200,1)
y_true <- 5*sin(x^4/5000000)+0.01*x
y <- y_true+ 2 * rnorm(100)

plot(x,y)
lines(x,y_true)

M <- kreg(x,y,bandwidth=bw1,grid=x)

lines(M$x,M$y,col="red")

y_id <- test(x,y,bandwidth=bw2,eps=eps)
lines(x,y_id,col="cyan")

mean((M$y-y_true)^2)
mean((y_id-y_true)^2)

for(k in (1:1000)){
M <- kreg(x,M$y,bandwidth=10,grid=x)
}
points(M$x,M$y,col="blue")




idempotent_kernel_regression <- function(x,y,bandwidth){
  l <- function(y_in){
    model <- kreg(x,y_in,grid=x,bandwidth=bandwidth)
	y_out <- model$y
  return(mean(abs(y_out-y)^2)  + 0.00*mean(abs(y_in-y)))}
  
  result <- optim(par=y,fn=l,control=list(maxit=10^8))
  return(kreg(x=x,y=result$par,grid=x,bandwidth=bandwidth)$y)
  }
  
  
test <- function(x,y,bandwidth,eps=10^-9){
   K <- gausskernel(X=x,sigma=bandwidth)
   model <- list(Q=t(K)%*%K +diag(rep(eps,length(x))),obj=-2*t(y)%*%K)
   model$A <- array(1,c(1,length(x)))
   model$sense <-"<="
   model$rhs <- sum(y)+1000000
   model$lb <- rep(-Inf,length(x))
   model$ub <- rep(Inf,length(x))
  # diag(rep(1,length(x)))
  # model$sense=rep("<=",length(x))
   #model$rhs <- rep(10000,length(x))
   solution <- try(solution <- gurobi(model,list(outputflag=0)), silent=TRUE)
   while(is.null(solution$x) | solution$status=="SUBOPTIMAL"){
    eps <- eps * 9/8
	print(eps)
	model$Q <- t(K)%*%K +diag(rep(eps,length(x)))
	
    solution <- try(gurobi(model,list(outputflag=0)), silent=TRUE)
   }
     solution <<- solution
   print(solution$status)
   return(K%*%solution$x)
   }
