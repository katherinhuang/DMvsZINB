load("C:/Users/kathe/OneDrive/00 Universität/Masterarbeit/R_Masterarbeit/Data/DM_results/DM.liver.4mcmc.RData")
load("C:/Users/kathe/OneDrive/00 Universität/Masterarbeit/R_Masterarbeit/Data/DM_results/DM.mela.results.4MCMC.RData")

beta <- list()
for( i in 1:4){
  beta[[i]]<-  mylist[[i]][[1]][["beta"]]
}



selected <- which(selected.DM)
fortraces = selected[sample(length(selected), 10)]
plot.ts(beta[[1]][,fortraces], main = "Some selected beta traceplots",
        xlab = "iteration", ylab = "")



inx <- 163
in2 <- 300
in4 <- 45

df <- NULL
for (i in 1:4){
  df <- cbind(df, beta[[i]][,inx])
}

colnames(df) <- c("MCMC1", "MCMC2", "MCMC3", "MCMC4")

df <- as.data.frame(df)


df1 <- NULL
for (i in 1:4){
  df1 <- cbind(df1, beta[[i]][,in2])
}

colnames(df1) <- c("MCMC1", "MCMC2", "MCMC3", "MCMC4")

df1 <- as.data.frame(df1)

df2 <- NULL
for (i in 1:4){
  df2 <- cbind(df2, beta[[i]][,in4])
}

colnames(df2) <- c("MCMC1", "MCMC2", "MCMC3", "MCMC4")

df2 <- as.data.frame(df2)

par(mfrow=c(3,1))

ggplot(data = df, aes(x = 1:2600)) + 
  geom_line(aes( y = MCMC1), color = "darkgreen") + 
  geom_line(aes( y = MCMC2), color="steelblue") + 
  geom_line(aes( y = MCMC3), color="darkred") + 
  geom_line(aes( y = MCMC4)) + 
  ggtitle("Traceplot of beta values for Roseburia_intestinalis")  

ggplot(data = df1, aes(x = 1:2600)) + 
  geom_line(aes( y = MCMC1), color = "darkgreen") + 
  geom_line(aes( y = MCMC2), color="steelblue") + 
  geom_line(aes( y = MCMC3), color="darkred") + 
  geom_line(aes( y = MCMC4)) + 
  ggtitle("Traceplot of beta values for s__Parabacteroides_unclassified")  

ggplot(data = df2, aes(x = 1:2600)) + 
  geom_line(aes( y = MCMC1), color = "darkgreen") + 
  geom_line(aes( y = MCMC2), color="steelblue") + 
  geom_line(aes( y = MCMC3), color="darkred") + 
  geom_line(aes( y = MCMC4)) + 
  ggtitle("Traceplot of beta values for f__Coriobacteriaceae")  


colnames(Y.mela)[c(inx, in2, in4)]


----------------
  
  x <-  list( mela.2[[1]][["mu0_store"]])
for( i in 1:4){
  x[[i]]<-  mela.2[[1]][["mu0_store"]]
}

mu[[1]]
  
 