args <- commandArgs(TRUE)
run <- as.numeric(args[1])

results <- NULL

library(MoBPS)
library(RandomFieldsUtils)
library(miraculix)

RandomFieldsUtils::RFoptions(cores = 1)

print(RFoptions()$basic$cores)

for(setting in 1:500){

  h2 <- 0.3 # heritability
  n_cow <- 3000 # number of new cows generated per year
  # n_test <- 1000 # number of cows used for testing
  # n_bull <- 100 # number of new bulls generated (This should be variable!)
  # n_bull_sel <- 10 # number of sires used for the next generation (This should be variable!)
  # n_test <- sample(1000:1500,1)
  # n_bull <- floor((10000000-4000*n_test)/3000)
  # 
  n_bull <-sample(150:200,1)
  n_test <- floor((10000000-3000*n_bull)/4000)
  n_bull_sel <- sample(10:15, 1)
  
  if(4000*n_test + 3000 * n_bull <10000000){
  print(c(n_test, n_bull, n_bull_sel))
  randomSeed <- sample(1:100, 1)


  {


    population <- creating.diploid(map=MoBPSmaps::map_cattle3,
                                   nindi= 100,
                                   n.additive = 1000,
                                   mean.target =100,
                                   var.target = 10, miraculix=TRUE)

    population <- breeding.diploid(population,
                                   heritability = h2,
                                   display.progress = FALSE)

    population <- breeding.diploid(population,
                                   breeding.size = c(n_bull, 0), ##hier
                                   name.cohort = "Testbullen_0",
                                   display.progress = FALSE)

    population <- breeding.diploid(population,
                                   breeding.size = c(n_bull, 0), ##hier
                                   name.cohort = "Testbullen_-1",
                                   display.progress = FALSE)

    population <- breeding.diploid(population,
                                   breeding.size = c(n_bull_sel, 0), ##hier
                                   name.cohort = "Zuchtbullen_0",
                                   display.progress = FALSE)

    population <- breeding.diploid(population,
                                   breeding.size = c(0,n_cow),
                                   name.cohort = "Kuehe_0",
                                   display.progress = FALSE)

    population <- breeding.diploid(population,
                                   breeding.size = c(0,n_cow),
                                   name.cohort = "Kuehe_-1",
                                   display.progress = FALSE)

    population <- breeding.diploid(population,
                                   breeding.size = c(0,n_cow),
                                   name.cohort = "Kuehe_-2",
                                   display.progress = FALSE)

  }

  # 5 Generations of Burn-in to initialize the breeding program
  for(index in 1:15){


    population <- breeding.diploid(population,
                                   breeding.size = c(n_bull, 0), #hier
                                   selection.m.cohorts = paste0("Zuchtbullen_", index-1),
                                   selection.f.cohorts = paste0("Kuehe_", index-1),
                                   max.offspring = c(Inf,1),
                                   name.cohort = paste0("Testbullen_", index),
                                   display.progress = FALSE)

    population <- breeding.diploid(population,
                                   breeding.size = c(0,n_test),
                                   selection.m.cohorts = paste0("Testbullen_", index-1),
                                   selection.f.cohorts = paste0("Kuehe_", index-1),
                                   # max.offspring = c(ceiling(n_test/n_bull),1),
                                   max.offspring = c(Inf,1),
                                   name.cohort = paste0("TestToechter_", index),
                                   display.progress = FALSE)


    if(index>1){ #phenoype the TestToechter if they are above one year old
      population <- breeding.diploid(population,
                                     phenotyping.cohorts = paste0("TestToechter_", index-1))
    }



    if(index>2){
      population <- breeding.diploid(population,
                                     offspring.bve.parents.cohorts = paste0("Testbullen_", index-2), #Cohorts to consider to derive phenotype from offspring phenotypes
                                     bve = TRUE,
                                     bve.cohorts = paste0("Testbullen_", index-2),
                                     input.phenotype = "off") #Select what to use in BVE, offspring phenotype ("off")
    }

    population <- breeding.diploid(population, breeding.size = c(n_bull_sel,0), #hier
                                   selection.size = c(n_bull_sel,0), #hier
                                   copy.individual.m = TRUE,
                                   selection.criteria = "bve",
                                   selection.m.cohorts =  paste0("Testbullen_", index-2),
                                   name.cohort = paste0("Zuchtbullen_",index),
                                   display.progress = FALSE)

    if(index>1){
      population <- breeding.diploid(population,
                                     bve=TRUE,
                                     relationship.matrix = "kinship",#Method to calculate relationship matrix for the breeding value estimation
                                     bve.cohorts = paste0("Kuehe_", c(index-1, index-2, index-3)))

      #1 "Kuehe_0"  "Kuehe_-1" "Kuehe_-2"
      #2 "Kuehe_1"  "Kuehe_0"  "Kuehe_-1"
      #3 "Kuehe_2" "Kuehe_1" "Kuehe_0"
      #4 "Kuehe_3" "Kuehe_2" "Kuehe_1"
      # ...
      ###
    }


    population <- breeding.diploid(population, breeding.size = c(0, n_cow),
                                   selection.size = c(n_bull_sel, n_cow*2),
                                   selection.criteria = c("bve"),
                                   selection.m.cohorts = paste0("Zuchtbullen_", index-1),
                                   selection.f.cohorts = paste0("Kuehe_", c(index-1, index-2, index-3)),
                                   max.offspring = c(Inf,1),
                                   name.cohort=paste0("Kuehe_", index),
                                   display.progress = FALSE)

    population <- breeding.diploid(population,
                                   phenotyping.cohorts = paste0("Kuehe_", index-1))

    if(index==5){
      population <- bv.standardization(population,  #Function to get mean and genetic variance of a trait to a fixed value
                                       cohorts=paste0("Kuehe_", index),
                                       adapt.bve = TRUE, #Modify previous breeding value estimations by scaling
                                       adapt.pheno = TRUE, #Modify previous phenotypes by scaling
                                       mean.target = 100,
                                       var.target = 10)
    }


  }

  acc1 <- acc2 <- 0
  for(index3 in 4:13){
    acc1 <- acc1 + cor(get.bv(population, cohorts=paste0("Testbullen_",index3))[1,], get.pheno.off(population, cohorts=paste0("Testbullen_",index3))[1,], use="complete.obs")
    acc2 <- acc2 + cor(get.bv(population, cohorts=paste0("Testbullen_",index3))[1,], get.bve(population, cohorts=paste0("Testbullen_",index3))[1,], use="complete.obs")
  }
  bv <- mean(get.bv(population, cohorts=paste0("Kuehe_", 15)))
  kin <- rbind(kinship.emp.fast(population=population, cohorts=paste0("Kuehe_", c(15, 14, 13))))
  results <- rbind(results, c(randomSeed, n_test, n_bull, n_bull_sel, bv, kin, acc1, acc2))
  save(file=paste0("results_", run, ".RData"), list=c("results"))

}


}

