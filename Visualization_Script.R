################################################################################
##################################Visualization#################################
results_full <- NULL

#first iteration results
load("first_iteration.RData")

set.seed(2091367)

##using all data and using smoothing
###Visualizing each objective separately

n_cat <- seq(100,700, by =3)
n_sel <- seq(3, 30, by=1)

rel1 <- rel2 <- ir_sim <- matrix(0, nrow=length(n_cat), ncol=length(n_sel))


#x1 = testbulls, x2 = proven sires, x3 = test_daughters

###genetic gain
approx_f <- function(x1, x2, x3, sd = 30, r=5){
  dist <- sqrt((results[,2] - x3)^2+ (results[,3] - x1)^2 + (30*(results[,4] -x2))^2 )
  weigth <- dnorm(dist, sd=sd)
  
  return (sum(results[,r] * weigth / sum(weigth)))
}

###Inbreeding rate
approx_f <- function(x1, x2,x3, sd = 10, r=6){
  dist <- sqrt((results[,2] - x3)^2+ (results[,3] - x1)^2 + (30*(results[,4] -x2))^2 )
  weigth <- dnorm(dist, sd=sd)
  kin <- -(sum(results[,r] * weigth / sum(weigth)))
  return (-kin)
}


###Composite function
approx_f <- function(x1, x2, x3){
  dist <- sqrt((results[,2] - x3)^2+ (results[,3] - x1)^2 + (30*(results[,4] -x2))^2 )
  weigth <- dnorm(dist, sd=30)  # reduce it by half at each round
  weigth2 <- dnorm(dist, sd=10) # reduce it by half at each round
  gain <- sum(results[,5] * weigth / sum(weigth))
  kin <- sum(results[,6] * weigth2 / sum(weigth2))
  return ((gain - 50*kin)) } # higher weight on inbreeding rate

for(index in 1:length(n_cat)){
  print(index)
  for(index2 in 1:length(n_sel)){
    ir_sim[index,index2] <- approx_f(n_cat[index], n_sel[index2], (10000000-3000*n_cat[index])/4000)
    rel1[index,index2] <- n_cat[index]
    rel2[index,index2] <- n_sel[index2]
  }
}


par(mar=c(5,6,4,3) + 0.1)
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
                   xlim = c(100, 700), 
                   ylim = c(3.09, 30))

#################################################################################
################################################################################
############################# Zoom-In ##########################################
par(mfrow=c(1,1))
results_full <- NULL

levels = seq(106, 108.3, length.out=20)

load("first_iteration.RData")

approx_f <- function(x1, x2, x3, bw=1){
  dist <- sqrt((results[,2] - x3)^2+ (results[,3] - x1)^2 + (30*(results[,4] -x2))^2 )
  weigth <- dnorm(dist, sd=30*bw)
  weigth2 <- dnorm(dist, sd=10*bw)
  gain <- sum(results[,5] * weigth / sum(weigth))
  kin <- sum(results[,6] * weigth2 / sum(weigth2))
  return ((gain - 50*kin)) }

n_cat <- seq(100,700, by =3)
n_sel <- seq(3, 30, by=1)

ir <- rel <- rel1 <- rel2 <-  ir_sim <-
  r_tab <-  r_pheno <- r_bve <- matrix(0, nrow=length(n_cat), ncol=length(n_sel)) 


bw =1
for(index in 1:length(n_cat)){
  print(index)
  for(index2 in 1:length(n_sel)){
    rel[index,index2] <- paste0(n_cat[index],",",n_sel[index2])
    rel1[index,index2] <- n_cat[index]
    rel2[index,index2] <- n_sel[index2]
    
    ir_sim[index,index2] <- approx_f(n_cat[index], n_sel[index2], (10000000-3000*n_cat[index])/4000, bw=bw)
  }
}

library(akima)
fld <- interp(x=rel1,y=rel2,z=ir_sim)


max1 = which.max(ir_sim)
rel1[max1]
rel2[max1]


par(mar=c(5,6,4,3) + 0.1)
a = filled.contour(fld, 
                   plot.title = {par(cex.main=2,
                                     adj=0.5,
                                     cex.lab=3, 
                                     cex.axis =2)},
                   levels = levels,
                   xlim = c(100, 400), 
                   ylim = c(8, 25), xaxs="i", yaxs="i",
                   plot.axes = {points(c(rel1[max1],189,172), c(rel2[max1],12,14), lwd=3, pch=19); axis(1);axis(2)})


bx = c(120, 260)
by = c(10,20)

load("second_iteration.RData")

results = rbind(results, LOAD_NEW)
bw =0.5
for(index in 1:length(n_cat)){
  print(index)
  for(index2 in 1:length(n_sel)){
    rel[index,index2] <- paste0(n_cat[index],",",n_sel[index2])
    rel1[index,index2] <- n_cat[index]
    rel2[index,index2] <- n_sel[index2]
    
    ir_sim[index,index2] <- approx_f(n_cat[index], n_sel[index2], (10000000-3000*n_cat[index])/4000, bw=bw)
  }
}

fld <- interp(x=rel1,y=rel2,z=ir_sim)
max1 = which.max(ir_sim)
rel1[max1]
rel2[max1]


par(mar=c(5,6,4,3) + 0.1)
a = filled.contour(fld, 
                   plot.title = {par(cex.main=2,
                                     adj=0.5,
                                     cex.lab=3, 
                                     cex.axis =2)},
                   levels = levels,
                   xlim = c(100, 400), 
                   ylim = c(8, 25), xaxs="i", yaxs="i",
                   plot.axes = {points(c(rel1[max1]),
                                       c(rel2[max1]), lwd=3, pch=19);
                     polygon(c(bx[1], bx[2], bx[2], bx[1]),
                             c(by[1], by[1], by[2],by[2]), lwd=2);
                     axis(1);axis(2)})

bx = c(150, 200)
by = c(10,15)

load("third_iteration.RData")

results = rbind(results, third_sim)
bw =0.25
for(index in 1:length(n_cat)){
  print(index)
  for(index2 in 1:length(n_sel)){
    rel[index,index2] <- paste0(n_cat[index],",",n_sel[index2])
    rel1[index,index2] <- n_cat[index]
    rel2[index,index2] <- n_sel[index2]
    
    ir_sim[index,index2] <- approx_f(n_cat[index], n_sel[index2], (10000000-3000*n_cat[index])/4000, bw=bw)
  }
}


fld <- interp(x=rel1,y=rel2,z=ir_sim)

par(mar=c(5,6,4,3) + 0.1)
a = filled.contour(fld, 
                   plot.title = {par(cex.main=2,
                                     adj=0.5,
                                     cex.lab=3, 
                                     cex.axis =2)},
                   levels = levels,
                   xlim = c(100, 400), 
                   ylim = c(8, 25), xaxs="i", yaxs="i",
                   plot.axes = {points(c(179),
                                       c(13), lwd=3, pch=19);
                     polygon(c(bx[1], bx[2], bx[2], bx[1]),
                             c(by[1], by[1], by[2],by[2]), lwd=2);
                     axis(1);axis(2)})


############# Comparing estimate of number of needed simulations ###############
library(data.table)
library(tidyverse)

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




gain_1000 <- genetic_gain_list_samples %>% filter(Nr.Sim == "1000")
gain_5000 <- genetic_gain_list_samples %>% filter(Nr.Sim == "5000")
gain_10000 <- genetic_gain_list_samples %>% filter(Nr.Sim == "10000")
gain_20000 <- genetic_gain_list_samples %>% filter(Nr.Sim == "20000")
gain_40000 <- genetic_gain_list_samples %>% filter(Nr.Sim == "40000")
gain_60000 <- genetic_gain_list_samples %>% filter(Nr.Sim == "60000")


library(MASS)
# Calculate the density using kde2d
col_palette <- hcl.colors(5, "YlOrRd", rev = TRUE)
a <- kde2d(gain_60000$bull, gain_60000$sel, n=200)
image(a, breaks = seq(0, max(a$z), length.out=11), plot.title = {par(cex.main=2,
                                                                     adj=0.5,
                                                                     cex.lab=3,
                                                                     cex.axis =2)},
      col = hcl.colors(10, "YlOrRd", rev = TRUE))

my_seq <- (seq(0, max(a$z), length.out=5))

legend("topleft", fill = col_palette,
       legend = format(my_seq, digits = 1),
       cex = 1.5,
       title = "Density")

