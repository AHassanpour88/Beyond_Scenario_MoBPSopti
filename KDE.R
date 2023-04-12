
############# Comparing estimate of number of needed simulations ###############
library(data.table)
library(tidyverse)
library(MASS)

load("Bootstrapping_100rep.Rdata")
head(genetic_gain_list_samples)

genetic_gain_list_samples$dtest <- (10000000-3000*genetic_gain_list_samples$n_bull)/4000

colnames(genetic_gain_list_samples) <- c("Nr.Sim", "repetition", "bull", "sel", "target", "test")


### only keep the samples you want to see ###
genetic_gain_list_samples <-
  genetic_gain_list_samples[
    genetic_gain_list_samples$Nr.Sim != 100 &
      genetic_gain_list_samples$Nr.Sim != 200 &
      genetic_gain_list_samples$Nr.Sim != 300 &
      genetic_gain_list_samples$Nr.Sim != 400 &
      genetic_gain_list_samples$Nr.Sim != 500 &
      genetic_gain_list_samples$Nr.Sim != 600 &
      genetic_gain_list_samples$Nr.Sim != 700 &
      genetic_gain_list_samples$Nr.Sim != 800 &
      genetic_gain_list_samples$Nr.Sim != 900 &
      genetic_gain_list_samples$Nr.Sim != 2000 &
      genetic_gain_list_samples$Nr.Sim != 3000 &
      genetic_gain_list_samples$Nr.Sim != 4000 &
      genetic_gain_list_samples$Nr.Sim != 6000 &
      genetic_gain_list_samples$Nr.Sim != 7000 &
      genetic_gain_list_samples$Nr.Sim != 8000 &
      genetic_gain_list_samples$Nr.Sim != 9000 &
      genetic_gain_list_samples$Nr.Sim != 30000, ]


genetic_gain_list_samples = as.matrix(genetic_gain_list_samples)
storage.mode(genetic_gain_list_samples) = 'numeric'

gain_1000 <- genetic_gain_list_samples[genetic_gain_list_samples[,1]=="1000",]
gain_10000 <- genetic_gain_list_samples[genetic_gain_list_samples[,1]=="10000",]
gain_20000 <- genetic_gain_list_samples[genetic_gain_list_samples[,1]=="20000",]
gain_60000 <- genetic_gain_list_samples[genetic_gain_list_samples[,1]=="60000",]

gain = list(gain_1000, gain_10000, gain_20000, gain_60000)
#
# bla <- which( 150<gain_1000[,3] & gain_1000[,3]<200 & 10<gain_1000[,4] & gain_1000[,4]<15)
# gain_1000[bla,]


par(mar=c(5,6,4,3) + 0.1)
for(index3 in 1:4){

  approx_f <- function(x1, x2, x3, sd = 15, r=5, results = NULL){

    dist <- sqrt((results[,6] - x3)^2+ (results[,3] - x1)^2 + (30*(results[,4] -x2))^2 )
    weigth <- dnorm(dist, sd=sd)

    return ( sum(weigth) / length(weigth)/ sd )
  }

  n_cat <- seq(100,700, by =3)
  n_sel <- seq(3, 30, by=1)

  rel <- rel1 <- rel2 <- ir_sim <- matrix(0, nrow=length(n_cat), ncol=length(n_sel))

  for(index in 1:length(n_cat)){
    print(index)
    for(index2 in 1:length(n_sel)){
      ir_sim[index,index2] <- approx_f(n_cat[index], n_sel[index2], (10000000-3000*n_cat[index])/4000, results = gain[[index3]])
      rel1[index,index2] <- n_cat[index]
      rel2[index,index2] <- n_sel[index2]
    }
  }

  require(akima)

  fld <- interp(x=rel1,y=rel2,z=ir_sim,
                linear = T, extrap=T, duplicate = "mean", dupfun = NULL,
                nx = 400, ny = 400,
                jitter = 10^-15, jitter.iter = 10, jitter.random = FALSE)


  max1 = which.max(ir_sim)
  rel1[max1]
  rel2[max1]

  a = filled.contour(fld,
                     plot.title = {par(cex.main=2,
                                       adj=0.5,
                                       cex.lab=3,
                                       cex.axis =2)},
                     zlim = c(0, 0.0003), # this is not calculated but just me looking at what maximum values across plots where

                     xlim = c(100, 400),
                     ylim = c(4, 30))



}
