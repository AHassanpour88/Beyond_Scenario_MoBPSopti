library(data.table)
library(tidyverse)

load("results_full.RData")

# sensitivity -------------------------------------------------------------

sensivity_analysis <- function(data
                       , sample_size 
                       , rep 
                       , n_bull = seq(100, 700, by = 1)
                       , n_sel = seq(6, 30, by = 1)
                       ) {

sim_n <- function(data, sample_size, rep) {
  # loop with the name of the samples
  1:length(sample_size) %>% 
    # adding the name of the list
    purrr::set_names(sample_size) %>% 
    # map is sorta equal to loop, it apply a function to each element of a vector or list
    # map always returns a list
    map(function(x) {
      
      # access the size of the samples
      size <- sample_size[x] %>% as.numeric()
      
      ## do the sampling for number of repetitions
      1:rep %>% 
        map( function(y) {
          
          # sampling from initial 60000 simulations depends on the sample size
          data[sample(nrow(data), size = size, replace=FALSE),] 
          
        }) 
      
    })
  
} 

# use the function for smoothing procedure, with the same name to don't redo it, it can have different length 1000,20000
mysim <- sim_n(data
               , sample_size 
               , rep
               )

allsim <- 1:length(mysim) %>% 
  # add the names again
  purrr::set_names(sample_size) %>% 
  map(function(x) {
    
    ## extracting the list by name
    results <- mysim[[names(mysim)[[x]]]]
    
    1:rep %>% map(function(p) {
      
      #smoothing function
      approx_f <- function(x1, x2, x3, bw=1){
        dist <- sqrt((results[[p]][,2] - x3)^2+ (results[[p]][,3] - x1)^2 + (30*(results[[p]][,4] -x2))^2 )
        weigth <- dnorm(dist, sd=30*bw)
        weigth2 <- dnorm(dist, sd=10*bw)
        gain <- sum(results[[p]][,5] * weigth / sum(weigth))
        kin <- sum(results[[p]][,6] * weigth2 / sum(weigth2))
        both <- (gain - 50*kin)
        return (both)
      }
      
      rel1 <- rel2 <-  ir_sim <- matrix(0, nrow=length(n_bull), ncol=length(n_sel)) 
      
      #apply smoothing for n_bull and n_sel in search space
      for(index in 1:length(n_bull)){
        print(index)
        for(index2 in 1:length(n_sel)){
          rel[index,index2] <- paste0(n_bull[index],",",n_sel[index2])
          rel1[index,index2] <- n_bull[index]
          rel2[index,index2] <- n_sel[index2]
          ir_sim[index,index2] <- approx_f(n_bull[index], 
                                           n_sel[index2], 
                                           (10000000-3000*n_bull[index])/4000,
                                           bw=1)
        }
      }
      
      max1 = which.max(ir_sim)
      
      data.table(
        scenario = names(mysim)[[x]],
        repetition = p,
        n_bull = rel1[max1],
        n_sel = rel2[max1],
        target = max(ir_sim)
      )
      
    }) %>% 
      #join the repetitions by scenario
      bind_rows()
    
  }) 

}

rsl <- sensivity_analysis(data = results
                  , rep = 50
                  , sample_size = c(1000, 10000,20000,30000,40000,50000,60000)
                  , n_bull = 
                  )

save(rsl, file = "Results_50rep.Rdata")

rsl %>% 
  #make it as a table
  bind_rows() %>% 
  mutate(across(scenario, as.factor)) %>% 
  ggplot(aes(x = scenario, y = target, group = scenario, fill = scenario)) +
  geom_boxplot() +
  theme_classic()

rsl %>% 
  #make it as a table
  bind_rows() %>% 
  mutate(across(scenario, as.factor)) %>% 
ggplot(aes(x = ordered(scenario),
           y = abs(target-107.8),
           color=scenario,
           group = scenario))+
  geom_violin(scale = "width")+
  geom_boxplot(width = 0.2)+
  geom_point(position = position_jitter(width = 0.2,height = 0))+ 
  labs(x = "Scenarios",
       y = "Target") +
  theme_classic()



rsl %>% 
  #make it as a table
  bind_rows() %>% 
  mutate(across(scenario, as.factor)) %>% 
  ggplot(aes(x = ordered(scenario),
             y = abs(n_bull),
             color=scenario,
             group = scenario))+
  geom_violin(scale = "width")+
  geom_boxplot(width = 0.2)+
  geom_hline(yintercept = 224)+
  geom_point(position = position_jitter(width = 0.2,height = 0))+
  theme_classic()+ 
  labs(x = "Scenarios",
       y = "Number of Bulls") +
  theme(axis.text = element_text(size = 12),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1.2))  +
  scale_y_continuous(breaks=seq(100,700,by=100))
