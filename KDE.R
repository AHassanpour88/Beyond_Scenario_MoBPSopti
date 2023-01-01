
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




gain_1000 <- genetic_gain_list_samples %>% filter(Nr.Sim == "1000")
gain_5000 <- genetic_gain_list_samples %>% filter(Nr.Sim == "5000")
gain_10000 <- genetic_gain_list_samples %>% filter(Nr.Sim == "10000")
gain_20000 <- genetic_gain_list_samples %>% filter(Nr.Sim == "20000")
gain_40000 <- genetic_gain_list_samples %>% filter(Nr.Sim == "40000")
gain_60000 <- genetic_gain_list_samples %>% filter(Nr.Sim == "60000")



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
