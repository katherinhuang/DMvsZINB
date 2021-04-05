library(IntegrativeBayes)
library(dirmult)
library(dplyr)
library(ggplot2)
library(xtable)
source('C:/Users/kathe/OneDrive/00 Universität/Masterarbeit/R_Masterarbeit/dmbvs/dmbvs-master/code/helper_functions.R')
source('C:/Users/kathe/OneDrive/00 Universität/Masterarbeit/R_Masterarbeit/Jiang helper.R')

### load the data from 4 MCMCs 
load("C:/Users/kathe/OneDrive/00 Universität/Masterarbeit/R_Masterarbeit/Data/DM_results/DM.liver.4mcmc.RData")
load("C:/Users/kathe/OneDrive/00 Universität/Masterarbeit/Data/DM_results/DM.mela.results.4MCMC.RData")

# transform DM into gamma_PPI ---------------------------------------------



beta1 <-  mylist[[1]][[1]][["beta"]]
beta2 <-  mylist[[2]][[1]][["beta"]]
beta3 <-  mylist[[3]][[1]][["beta"]]
beta4 <-  mylist[[4]][[1]][["beta"]]

mppi1 = colMeans((beta1 != 0) + 0)
mppi2 = colMeans((beta2 != 0) + 0)
mppi3 = colMeans((beta3 != 0) + 0)
mppi4 = colMeans((beta4 != 0) + 0)

mppi1 <- matrix(mppi1, nrow = 7, ncol = 528, byrow = F )
mppi2 <- matrix(mppi2, nrow = 7, ncol = 528, byrow = F )
mppi3 <- matrix(mppi3, nrow = 7, ncol = 528, byrow = F )
mppi4 <- matrix(mppi4, nrow = 7, ncol = 528, byrow = F )

gammappi1 <- colMeans(mppi1)
gammappi2 <- colMeans(mppi2)
gammappi3 <- colMeans(mppi3)
gammappi4 <- colMeans(mppi4)

gamma.df <- cbind(gammappi1, gammappi2, gammappi3, gammappi4)

gamma_ppi <- rowMeans(gamma.df)




# Manhattan plot ----------------------------------------------------------



sig.level = 0.05
gammaPPI = gamma_ppi
th.gamma = bfdr(gammaPPI, sig.level )[[2]]
selected.DM = bfdr(gammaPPI, sig.level )[[1]]

col.code <- ifelse(selected,'red','grey')
p.code <- ifelse(selected,19,20)
cex.code <- ifelse(selected,1,0.25)


#### plot ####
layout(rbind(1,2), heights  = c(7,1))  # put legend on bottom 1/8th of the chart
plot(round(gammaPPI,2), type='h', ylim = c(0,1),
     ylab = "Posterior Probability of Inclusion",xlab = "Feature Index",
     main = expression(paste("Liver cirrhosis: Evaluation of the Feature Selection for ", gamma )),
     sub = "4 different MCMC "
)
abline(h = th.gamma,lty = 2,col = 'darkred')

# plots:


points(gammaPPI,
       col = col.code,
       pch = p.code, 
       cex = cex.code)


text(x = 450, y = 0.95 ,                # Add text element
     paste("threshold = " , 0.9))

# setup for no margins on the legend
par(mar=c(0, 0, 0, 0))
c(bottom, left, top, right)
plot.new()
legend(x = 'center',c("Discriminating Feature (DF)","Non-discriminating Feature"),
       pch = c( 19,20),
       col = c('red','grey'),
       ncol = 2, 
       bty ="n" ,
       cex = 1.5)



# find out the taxa ---------------------------------------------------
load("C:/Users/kathe/OneDrive/00 Universität/Masterarbeit/Data/Real_data/cirrhosis_datasets.Rdata")
Y.liver <- Y.filter(Y_cirrhosis, zvec = z_cirrhosis, min.number = 2)[[2]]

liver.df <- cbind(colnames(Y.liver)[selected], round(gammaPPI[selected], digits = 4))

print(xtable(liver.df, type = "latex"), file = "DM.liver.tex")

### no idea how to get a good dentro gram. 

liver.gamma.dm.df <- as.data.frame( cbind(colnames(Y.liver), round(gammaPPI, digits = 4)))

hc <- hclust(dist(colnames(Y.liver)[selected]), method = "single")

plot(hc)




#find out how many of the DM selected are in the ZINB selected

selected.DM <- selected

sig.level = 0.05
th.gamma = bfdr(liver.fin.gamma_PPI, sig.level )[[2]]
selected.ZINB = bfdr(liver.fin.gamma_PPI, sig.level )[[1]]



sum((selected.DM  == 1 )  &( selected.ZINB == 1))
sum((selected.DM  == 0 )  &( selected.ZINB == 0))


colnames(Y.liver)[ind]



ind <- which((selected.DM  == 1 )  &( selected.ZINB == 1))

# Pearson Correlation coefficient -----------------------------------------

beta1 <-  mylist[[1]][[1]][["beta"]]
beta2 <-  mylist[[2]][[1]][["beta"]]
beta3 <-  mylist[[3]][[1]][["beta"]]
beta4 <-  mylist[[4]][[1]][["beta"]]

mppi1 = colMeans((beta1 != 0) + 0)
mppi2 = colMeans((beta2 != 0) + 0)
mppi3 = colMeans((beta3 != 0) + 0)
mppi4 = colMeans((beta4 != 0) + 0)

mppi1 <- matrix(mppi1, nrow = 7, ncol = 528, byrow = F )
mppi2 <- matrix(mppi2, nrow = 7, ncol = 528, byrow = F )
mppi3 <- matrix(mppi3, nrow = 7, ncol = 528, byrow = F )
mppi4 <- matrix(mppi4, nrow = 7, ncol = 528, byrow = F )

gammappi1 <- colMeans(mppi1)
gammappi2 <- colMeans(mppi2)
gammappi3 <- colMeans(mppi3)
gammappi4 <- colMeans(mppi4)



liver.DM.gammappi <- cbind(gammappi1 ,gammappi2,gammappi3,gammappi4)
liver.DM.beta <- cbind(mppi1,mppi2,mppi3,mppi4 )



getPairwisePearson(liver.DM.beta)

## for mela


is.numeric(35)

---------------------------------
        
        # setup -------------------------------------------------------------------
library(IntegrativeBayes)
source('C:/Users/kathe/OneDrive/00 Universität/Masterarbeit/R_Masterarbeit/dmbvs/dmbvs-master/code/helper_functions.R', echo=TRUE)


#load data

load("C:/Users/kathe/OneDrive/00 Universität/Masterarbeit/Data/DM_results/DM.liver.4mcmc.RData")
load("C:/Users/kathe/OneDrive/00 Universität/Masterarbeit/Data/ZINB_results/ZINB_liver_4MCMC_5000.RData")
load("C:/Users/kathe/OneDrive/00 Universität/Masterarbeit/Data/ZINB_results/ZINB_liver2_4MCMC_5000.RData")
load("C:/Users/kathe/OneDrive/00 Universität/Masterarbeit/Data/Real_data/cirrhosis_datasets.Rdata")
Y.liver <- Y.filter(Y_cirrhosis, zvec = z_cirrhosis, min.number = 2)[[2]]

# data preparation --------------------------------------------------------
## we only need to look at 14 taxa 
##

sig.level = 0.05
gammaPPI = gamma_ppi
th.gamma = bfdr(gammaPPI, sig.level )[[2]]
selected.DM = bfdr(gammaPPI, sig.level )[[1]]


sig.level = 0.05
th.gamma.ZINB = bfdr(liver.fin.gamma_PPI, sig.level )[[2]]
selected.ZINB = bfdr(liver.fin.gamma_PPI, sig.level )[[1]]


save( gamma_ppi, liver.fin.gamma_PPI, selected.DM, selected.ZINB, overlap.selected, ZINB.beta.mat, file = "compareB.RData" )


overlap.selected <- selected.DM & selected.ZINB   ####### this is what we need!

## ZINB Beta subset
ZINB.beta.mat <- liver.beta.mat[,overlap.selected] 


## DM Beta subset


hist(DM.beta.mat)
DM.beta.mat <- (beta1+beta2+beta3+beta4)/4

D <- DM.beta.mat[selected.DM]
DM.beta.mat <- DM.beta.mat[,overlap.selected]
DM.beta.mat <- DM.beta.mat[,ind]

summary(as.vector(DM.beta.mat))
summary(D)
####

MPPI.ZINB = data.frame(expand.grid(covariates = colnames(X_cirrhosis_processed),
                                   taxa = colnames(Y.liver)[ind]),
                       mppi = as.vector(liver.delta.mat[,ind]),
                       beta =as.vector(liver.beta.mat[,ind]))
thz <- BayFDR(liver.delta.mat, 0.05)
MPPI.ZINB <- subset(MPPI.ZINB, mppi > thz )

summary(MPPI.ZINB$beta)    #Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-0.86561 -0.27333 -0.16745  0.07807  0.41908  1.34541 

association_plot(MPPI.ZINB[,-5], graph_layout = "circular", main = "Sample Results")



########################################################################################################################



    
      Y.mela <- Y.filter(Y_melanoma, zvec = z_melanoma, min.number = 2)[[2]]  

      MPPI.DM = data.frame(expand.grid(covariates = colnames(X_metabolite_processed),
                                         taxa = colnames(Y.mela)),
                             mppi = gamma,
                             beta =DM.beta.mat)
      
      association_plot(MPPI.ZINB[,-5], graph_layout = "bipartite", main = "Sample Results")
      

mela.MPPI.pos <- subset(mela.MPPI , beta > 0 )
mela.MPPI.pos <- arrange(mela.MPPI.pos, desc(beta))
mela.MPPI.neg <- subset(mela.MPPI , beta < 0 )
mela.MPPI.neg <- arrange(mela.MPPI.neg, beta)

th <- bfdr(mela.fin.gamma_PPI, 0.05)[[2]]
association_plot(mela.MPPI.pos, mppi_threshold = th, graph_layout = "bipartite")
association_plot(mela.MPPI.neg, mppi_threshold = th, graph_layout = "bipartite")


# new  --------------------------------------------------------------------

## get beta: 

beta <- list()
for( i in 1:4){
  beta[[i]] <-  mylist[[i]][[1]][["beta"]]  
}


beta_posterior <- lapply(beta, colMeans) ############# use this one 

fun <- function(x){x != 0}


delta <- lapply(beta, fun)
delta <- lapply(delta, colMeans) ################ use this one 



B <- (beta_posterior[[1]] +  beta_posterior[[2]]  +  beta_posterior[[3]] + beta_posterior[[4]]) /4
summary(B) # -3.475827 -0.000229  0.000000  0.075827  0.023989  6.919350 
B <- matrix(B,  nrow = 9, ncol = 338, byrow = F )



B[,ind] <- bfdr(gamma, 0.05)[[2]]


D <- (delta[[1]] + delta[[2]] + delta[[3]] + delta[[4]]) /4
summary(D)
D <- matrix(D,  nrow = 7, ncol = 528, byrow = F )
th <- BayFDR(D, 0.05)


MPPI.DM= data.frame(expand.grid(covariates = colnames(X_cirrhosis_processed),
                                taxa = colnames(Y.liver)[ind]),
                    mppi = as.vector(D[,ind]),
                    beta =as.vector(B[,ind]))

MPPI.DM <- subset(MPPI.DM, mppi > th)
summary(MPPI.DM$beta) # -2.9957 -0.2663  0.1001  0.1763  0.6615  3.2899 

association_plot(MPPI.DM[,-5], graph_layout = "circular", main = "Sample Results")

###################
### now for mela: 

beta <- list()
for( i in 1:4){
  beta[[i]] <-  mylist[[i]][[1]][["beta"]]  
}


fun <- function(x){x != 0}

beta_posterior <- lapply(beta, colMeans) ############# use this one 
delta <- lapply(beta, fun)
delta <- lapply(delta, colMeans) ################ use this one 


B <- (beta_posterior[[1]] +  beta_posterior[[2]]  +  beta_posterior[[3]] + beta_posterior[[4]]) /4
summary(B) # -14.439821  -0.122713   0.000966  -0.180659   0.370432   5.131753 
B <- matrix(B,  nrow = 9, ncol = 338, byrow = F )


th.gamma <-bfdr(gamma, 0.05)[[2]]
ind <- which(gamma > th.gamma)
 



D <- (delta[[1]] + delta[[2]] + delta[[3]] + delta[[4]]) /4
summary(D)
D <- matrix(D,  nrow = 9, ncol = 338, byrow = F )
th <- BayFDR(D, 0.05)





ind

mela.MPPI= data.frame(expand.grid(covariates = colnames(X_metabolite_processed),
                                  taxa = colnames(Y.mela)[ind]),
                      mppi = as.vector(D[,ind]),
                      beta =as.vector(B[,ind]))



th <- BayFDR( D, 0.05)
mela.MPPI <- subset(mela.MPPI, mppi > th)

mela.pos <- subset(mela.MPPI, beta > 0 )
mela.neg <- subset(mela.MPPI, beta < 0 )

#par(mfrow = c(2, 1))

association_plot(mela.pos[,-5], graph_layout = "bipartite", main = "Sample Results")
association_plot(mela.neg[,-5], graph_layout = "bipartite", main = "Sample Results")
summary(mela.MPPI$beta) # -4.1913 -1.1202  0.6089 -0.3247  0.9475  1.3428

B[,ind]

### 


#convergence 

beta[[1]][, 78]






