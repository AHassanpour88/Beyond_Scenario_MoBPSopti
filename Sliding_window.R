library(data.table)
library(ggplot2)

  load(paste0("data/firs_iteration.RData")) #10,15 -121,219

n_window <- 10 # number of windows
share_window <- 5*2 # relative size of window --> 5*2 --> 1/5 of range
y <- seq(min(results[,3]),max(results[,3]),length.out=n_window)
z <- seq(min(results[,4]),max(results[,4]),length.out=n_window)

vars <- list()
n <- 1
#for(i_x in x){
for(i_y in y){
  for(i_z in z){
    vars[[n]] <- results[  results[,3] >= i_y - abs(y[1]-y[length(y)])/share_window & results[,3] <= i_y + abs(y[1]-y[length(y)])/share_window &
                             results[,4] >= i_z - abs(z[1]-z[length(z)])/share_window & results[,4] <= i_z + abs(z[1]-z[length(z)])/share_window ,]
    vars[[n]] <- data.table(
      y = i_y,
      z = i_z,
      sd = sd(vars[[n]][,5]  - 50*vars[[n]][,6])
    )
    n <- n+1
  }
  
  # }
  
}
vars <- rbindlist(vars)
vars$scenario <- rep("Not Smooth",100)
colnames(vars) <- c("bull", "sel","sd", "Scenario")
ggplot(vars,
       aes(
         x = y,
         y = z,
         fill = sd
       ))+
  geom_tile()

ggplot(vars,
       aes(
         x = sd
       ))+
  geom_histogram(bins=10)

################################################################################
n_cat <- seq(100,700, by =1)
n_sel <- seq(3, 30, by=1)

rel1 <- rel2 <- ir_sim <- matrix(0, nrow=length(n_cat), ncol=length(n_sel))

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

n_bull <- matrix(rel1, dimnames=list(t(outer(colnames(rel1), rownames(rel1), FUN=paste)), NULL))
n_sel <- matrix(rel2, dimnames=list(t(outer(colnames(rel2), rownames(rel2), FUN=paste)), NULL))
ir_sim <- matrix(ir_sim, dimnames=list(t(outer(colnames(ir_sim), rownames(ir_sim), FUN=paste)), NULL))
soomth_full <- cbind(n_bull,n_sel,ir_sim)
colnames(soomth_full) <- c("bull", "sel","full_smooth")


n_window <- 10 # number of windows
share_window <- 5*2 # relative size of window --> 5*2 --> 1/5 of range
yy <- seq(min(n_bull),max(n_bull),length.out=n_window)
zz <- seq(min(n_sel),max(n_sel),length.out=n_window)


vars1 <- list()
n <- 1
#for(i_x in x){
for(i_yy in yy){
  for(i_zz in zz){
    vars1[[n]] <- soomth_full[  soomth_full[,1] >= i_yy - abs(yy[1]-yy[length(yy)])/share_window & soomth_full[,1] <= i_yy + abs(yy[1]-yy[length(yy)])/share_window &
                                  soomth_full[,2] >= i_zz - abs(zz[1]-zz[length(zz)])/share_window & soomth_full[,2] <= i_zz + abs(zz[1]-zz[length(zz)])/share_window ,]
    vars1[[n]] <- data.table(
      yy = i_yy,
      zz = i_zz,
      sd = sd(vars1[[n]][,3])
    )
    n <- n+1
  }
  
  # }
  
}
vars1 <- rbindlist(vars1)
vars1$scenario <- rep("Smooth",100)
colnames(vars1) <- c("bull", "sel","sd", "Scenario")


both <- rbind(vars1,vars)


##Statistics of the standard deviation of the target score distributions across all parameter combinations
ggplot(both, aes(x = sd, fill = Scenario)) +                       # Draw overlaying histogram
  geom_histogram(position = "identity", alpha = 0.5, bins = 5) +
  theme_classic()+ 
  labs(x = "",
       y = "") +
  # coord_cartesian(ylim = c(111.2,112.7))+
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20,face = "bold"),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1.2),
        legend.position = "top")  +
  # scale_fill_brewer(palette = "YlOrRd") +
  
  scale_x_continuous(breaks=seq(0,2,by=0.5))+
  # scale_fill_manual(values = c("Not Smooth" = "yellow", "Smooth" = "red"))
  guides(col = guide_legend(title="Scenario"))

##Statistics of the standard deviation of the target score distributions across all parameter combinations
ggplot(both, aes(x = sd, fill = Scenario)) +                       # Draw overlaying histogram
  # geom_histogram(position = "identity", alpha = 0.5, bins = 5) +
  geom_density(alpha = 0.3)+
  theme_classic()+ 
  labs(x = "Standard deviation of the target score",
       y = "Density") +
  # coord_cartesian(ylim = c(111.2,112.7))+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20,face = "bold"),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1.2),
        legend.position = "top")  +
  # scale_fill_brewer(palette = "YlOrRd") +
  
  scale_x_continuous(breaks=seq(0,2,by=0.5))+
  guides(col = guide_legend(title="Scenario"))
###############################################################################

