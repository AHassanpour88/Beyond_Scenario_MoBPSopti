############## local variance using Kernel regression ###############

load("data/first_iteration.RData")

approx_f <- function(x1, x2, x3, sd1=30, sd2 = 10){
  dist <- sqrt((results[,2] - x3)^2+ (results[,3] - x1)^2 + (30*(results[,4] -x2))^2 )
  weigth <- dnorm(dist, sd=sd1)  
  weigth2 <- dnorm(dist, sd=sd2)
  gain <- sum(results[,5] * weigth / sum(weigth))
  kin <- sum(results[,6] * weigth2 / sum(weigth2))
  return (c(gain, kin)) } # higher weight on inbreeding rate


approx_var <- function(x1, x2, x3, sd1=30, sd2 = 10){
  dist <- sqrt((results[,2] - x3)^2+ (results[,3] - x1)^2 + (30*(results[,4] -x2))^2 )
  weigth <- dnorm(dist, sd=sd1)  
  weigth2 <- dnorm(dist, sd=sd2) 
  exp = approx_f(x1, x2, x3, sd1=sd1, sd2 = sd2)
  gain <- sum((results[,5]^2 - exp[1]^2) * weigth / sum(weigth))
  kin <- sum((results[,6]^2 - exp[2]^2) * weigth2 / sum(weigth2))
  return (c(gain, kin)) } # higher weight on inbreeding rate



sqrt(approx_var(100,3, (10000000-3000*200)/4000, sd1=30, sd2=30))
# [1] 1.32326726 0.02514313
sqrt(approx_var(100,3, (10000000-3000*200)/4000, sd1=10, sd2=10))

100 - (1.30021684 * 100) / 1.3552929   
100 - (0.01514772 * 100) / 0.0239601



sqrt(approx_var(103,30, (10000000-3000*200)/4000, sd1=30, sd2=30))
# [1] 1.32326726 0.02514313
sqrt(approx_var(103,30, (10000000-3000*200)/4000, sd1=10, sd2=10))


100 - (0.54115061  * 100) / 0.580868073     
100 - (0.00266228 * 100) / 0.002849617
