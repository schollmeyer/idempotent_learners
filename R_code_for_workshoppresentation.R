library(ggplot2)
set.seed(1234567)
n_sample <- 75
sd <- 5
x <- seq(1,200,length.out=75)
y_true <- 5*sin(x/10)+.1*x
y <- y_true+ sd * rnorm(n_sample)
pplot <- ggplot(data=data.frame(x=x,y=y), aes(x=x,y=y))
pplot + layer( mapping = NULL,   position = "identity",   stat="identity",   geom = "point") +labs(x = "x",y="y")




bandwidths_classical <- c(seq(1,100,length.out=1000),seq(101,10000,length.out=200))
classical_tuning <- bandwidth_selection(x=x,y_true=y_true,sd=sd,learner=classical_kernel_regression,bandwidths=bandwidths_classical,n_rep=1000)
bw_classical <- bandwidths_classical[which.min(classical_tuning$mses)]
#[1] 13.54955

mse_classical <- classical_tuning$mses[which.min(classical_tuning$mses)]
#[1] 1.612067

pplot <- ggplot(data=data.frame(x=bandwidths,y=idempotent_tuning$mses), aes(x=x,y=y))
pplot + layer( mapping = NULL,   position = "identity",   stat="identity",   geom = "point") +labs(x = "bandwidth",y="MSE")

classical_result <- classical_kernel_regression(x=x,y=y,bandwidth=bw_classical)

pplot <- ggplot(data=data.frame(x=x,y=y), aes(x=x,y=y))
pplot + layer( mapping = NULL,   position = "identity",   stat="identity",   geom = "point") +labs(x = "x",y="y") + layer( mapping = NULL,  data=data.frame(x=x,y=classical_result$fitted), position = "identity",   stat="identity",   geom = "line") + layer( mapping = NULL,   position = "identity",   stat="identity",   geom = "point") +labs(x = "x",y="y") + layer( mapping = NULL,  data=data.frame(x=x,y=idempotent_result$fitted), position = "identity",   stat="identity",   geom = "line")



## idempotent kernel regression

bandwidths_idempotent  <- c(seq(1,1000,length.out=500),seq(1001,10000,length.out=100))
idempotent_tuning <- bandwidth_selection(x=x,y_true=y_true,sd=sd,learner=idempotent_kernel_regression,bandwidths=bandwidths_idempotent,n_rep=1000)
bw_idempotent <- bandwidths_idempotent[which.min(idempotent_tuning$mses)]
idempotent_result <- idempotent_kernel_regression(x=x,y=y,bandwidth=bw_idempotent)
mse_idempotent <- idempotent_tuning$mses[which.min(idempotent_tuning$mses)]
mse_idempotent

saveRDS(idempotent_tuning,"idempotent_tuning.RDS")
saveRDS(classical_tuning,"classical_tuning.RDS")
##Plots


pplot <- ggplot(data=data.frame(x=x,y=y), aes(x=x,y=y))
pplot  + layer( mapping = NULL,   position = "identity",   stat="identity",geom = "point") +labs(x = "x",y="y") + layer( mapping = NULL,  data=data.frame(x=x,y=y_true), aes(color = "darkgreen",lwd=.8),position = "identity",   stat="identity",   geom = "line")


pplot <- ggplot(data=data.frame(x=x,y=y), aes(x=x,y=y))
pplot  + layer( mapping = NULL,   position = "identity",   stat="identity",geom = "point") +labs(x = "x",y="y") + layer( mapping = NULL,  data=data.frame(x=x,y=y_true), aes(color = "darkgreen",lwd=.8),position = "identity",   stat="identity",   geom = "line") + layer( mapping = NULL,  data=data.frame(x=x,y=classical_result$fitted), aes(color = "black",lwd=.8),position = "identity",   stat="identity",   geom = "line")
#


# Spiky kernel
x_new <- c(x,x+10^-6,x-10^-6)
y_new <- c(y,classical_result$fitted,classical_result$fitted)
o <- order(x_new)
pplot <- ggplot(data=data.frame(x=x,y=y), aes(x=x,y=y))
pplot  + layer( mapping = NULL,   position = "identity",   stat="identity",geom = "point") +labs(x = "x",y="y") + layer( mapping = NULL,  data=data.frame(x=x_new[o],y=y_new[o]), aes(color = "black",lwd=.4),position = "identity",   stat="identity",   geom = "line")

#  layer( mapping = NULL,   position = "identity",   stat="identity",geom = "point") +labs(x = "x",y="y")
#+ layer( mapping = NULL,  data=data.frame(x=x,y=y_true), aes(color = "darkgreen"",lwd=.8),position = "identity",   stat="identity",   geom = "line")

#+ layer( mapping = NULL,   position = "identity",   stat="identity",   geom = "point") +labs(x = "x",y="y")+layer( mapping = NULL,  data=data.frame(x=x,y=idempotent_result$fitted), aes(color = "darkblue",lwd=.8),position = "identity",   stat="identity",   geom = "line") #+ geom_line(color='darkblue')


pplot <- ggplot(data=data.frame(x=x,y=y), aes(x=x,y=y))
pplot + layer( mapping = NULL,   position = "identity",   stat="identity",geom = "point") +labs(x = "x",y="y") + layer( mapping = NULL,  data=data.frame(x=x,y=classical_result$fitted), aes(color = "black",lwd=.8),position = "identity",   stat="identity",   geom = "line") + layer( mapping = NULL,   position = "identity",   stat="identity",   geom = "point") +labs(x = "x",y="y")+layer( mapping = NULL,  data=data.frame(x=x,y=idempotent_result$fitted), aes(color = "darkblue",lwd=.8),position = "identity",   stat="identity",   geom = "line") #+ geom_line(color='darkblue')


pplot <- ggplot(data=data.frame(x=x,y=y), aes(x=x,y=y))
pplot + layer( aes(color = "black",lwd=.8),mapping = NULL,   data=data.frame(x=log2(bandwidths_classical),y=classical_tuning$mses), position = "identity",   stat="identity",geom = "line") +labs(x = "log2(bandwidth)",y="MSE") + layer( mapping = NULL,  data=data.frame(x=log2(bandwidths_idempotent),y=idempotent_tuning$mses), aes(color = "darkblue",lwd=.8),position = "identity",   stat="identity",   geom = "line")


#+ layer( mapping = NULL,   position = "identity",   stat="identity",   geom = "point") +labs(x = "log2(bandwidth",y="MSE")+layer( mapping = NULL,  data=data.frame(x=x,y=idempotent_tuning$mses), position = "identity",   stat="identity",   geom = "line") #+ geom_line(color='darkblue')

# Notizen:
# wahre Kurve einzeichnen?
# Plots: kernelregression mit GIB, mit non-greedy ib, MSE Kurven



# spiky indicator kernel

library(ggplot2)
set.seed(1234567)
n_sample <- 75
sd <- 5
x <- seq(1,200,length.out=500)
y_true <- (sin(x/10)+1)/4 +.2
y <- runif(500) <= y_true
y <- y*1


classical_result <- classical_kernel_regression(x=x,y=y,bandwidth=bw_classical)
y_fitted <- (classical_result$fitted) >=0.5
y_fitted=1*y_fitted
lines(x,classical_result$fitted>=0.5,col="red")
#y <- y_true+ sd * rnorm(n_sample)
pplot <- ggplot(data=data.frame(x=x,y=y), aes(x=x,y=y))
pplot  + layer( mapping = NULL,   position = "identity",   stat="identity",geom = "point") +labs(x = "x",y="y") + layer( mapping = NULL,  data=data.frame(x=x,y=y_fitted), aes(color = "darkblue",lwd=.8),position = "identity",   stat="identity",   geom = "line") + layer( mapping = NULL,   position = "identity",   stat="identity",geom = "point",data=data.frame(x=x,y=y_true)) +labs(x = "x",y="y")
