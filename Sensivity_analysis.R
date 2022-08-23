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
    # map is equal to loop
    map(function(x) {
      
      # access the size of the samples
      size <- sample_size[x] %>% as.numeric()
      
      ## do the sampling for number of repetitions
      1:rep %>% 
        # map creates a list itself, doesn't need to recreate an empty list
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
      
      approx_f <- function(x1, x2, x3, bw=1){
        dist <- sqrt((results[[p]][,2] - x3)^2+ (results[[p]][,3] - x1)^2 + (30*(results[[p]][,4] -x2))^2 )
        weigth <- dnorm(dist, sd=30*bw)
        weigth2 <- dnorm(dist, sd=10*bw)
        gain <- sum(results[[p]][,5] * weigth / sum(weigth))
        kin <- sum(results[[p]][,6] * weigth2 / sum(weigth2))
        both <- (gain - 50*kin)
        return (both)
      }
      
      rel <- rel1 <- rel2 <-  ir_sim <- matrix(0, nrow=length(n_bull), ncol=length(n_sel)) 
      
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
      #we want to join the repetitions by sceanrio
      bind_rows()
    
  }) 

}

rsl <- sensivity_analysis(data = results
                  , rep = 5
                  , sample_size = c(1000, 2000)
                  , n_bull = 
                  )

rsl %>% 
  #make it as a table
  bind_rows() %>% 
  mutate(across(scenario, as.factor)) %>% 
  ggplot(aes(x = scenario, y = target, group = scenario, fill = scenario)) +
  geom_boxplot() +
  theme_classic()



