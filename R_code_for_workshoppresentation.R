n_sample <- 200
sd <- 5
x <- seq(1,n_sample,1)
y_true <- 5*sin(x/10)+.1*x
y <- y_true+ sd * rnorm(n_sample)
pplot <- ggplot(data=data.frame(x=x,y=y), aes(x=x,y=y))
pplot + layer( mapping = NULL,   position = "identity",   stat="identity",   geom = "point") +labs(x = "x",y="y")




bandwidths <- seq(1,200,length.out=1000)
classical_tuning <- bandwidth_selection(x=x,y_true=y_true,sd=sd,learner=classical_kernel_regression,bandwidths=bandwidths,n_rep=500)
bw_classical <- bandwidths[which.min(classical_tuning$mses)]
#[1] 13.54955

mse_classical <- classical_tuning$mses[which.min(classical_tuning$mses)]
#[1] 1.612067

pplot <- ggplot(data=data.frame(x=bandwidths,y=idempotent_tuning$mses), aes(x=x,y=y))
pplot + layer( mapping = NULL,   position = "identity",   stat="identity",   geom = "point") +labs(x = "bandwidth",y="MSE")

classical_result <- classical_kernel_regression(x=x,y=y,bandwidth=bw_classical)

pplot <- ggplot(data=data.frame(x=x,y=y), aes(x=x,y=y))
pplot + layer( mapping = NULL,   position = "identity",   stat="identity",   geom = "point") +labs(x = "x",y="y") + layer( mapping = NULL,  data=data.frame(x=x,y=classical_result$fitted), position = "identity",   stat="identity",   geom = "line")



## idempotent kernel regression

bandwidths <- seq(1,8000,length.out=1000)
idempotent_tuning <- bandwidth_selection(x=x,y_true=y_true,sd=sd,learner=idempotent_kernel_regression,bandwidths=bandwidths,n_rep=100)
bw_idempotent <- bandwidths[which.min(idempotent_tuning$mses)]
idempotent_result <- idempotent_kernel_regression(x=x,y=y,bandwidth=bw_idempotent)
mse_idempotent <- idempotent_tuning$mses[which.min(idempotent_tuning$mses)]
mse_idempotent

saveRDS(idempotent_tuning,"idempotent_tuning.RDS")
saveRDS(classical_tuning,"classical_tuning.RDS")
##Plots


pplot <- ggplot(data=data.frame(x=x,y=y), aes(x=x,y=y))
pplot + layer( mapping = NULL,   position = "identity",   stat="identity",geom = "point") +labs(x = "x",y="y") + layer( mapping = NULL,  data=data.frame(x=x,y=classical_result$fitted), position = "identity",   stat="identity",   geom = "line") + layer( mapping = NULL,   position = "identity",   stat="identity",   geom = "point") +labs(x = "x",y="y")+layer( mapping = NULL,  data=data.frame(x=x,y=idempotent_result$fitted), position = "identity",   stat="identity",   geom = "line") #+ geom_line(color='darkblue')

# Notizen:
# wahre Kurve einzeichnen?
# Plots: kernelregression mit GIB, mit non-greedy ib, MSE Kurven



