
# Resource optimization in breeding program using Kernel regression
## Table of contents
[0 Introduction](https://github.com/AHassanpour88/Beyond_Scenario_MoBPSopti/edit/main/README.md#0-introduction) <br />
[1 Description of the underlying simulation case](https://github.com/AHassanpour88/Beyond_Scenario_MoBPSopti/edit/main/README.md#1-Simulation_Script) <br />
[2 Nadaraya-Watson kernel function estimator](https://github.com/AHassanpour88/Beyond_Scenario_MoBPSopti/edit/main/README.md#2-Visualization_Script) <br />
[3 Estimation of the number of needed simulations](https://github.com/AHassanpour88/Beyond_Scenario_MoBPSopti/edit/main/README.md#4-Kernel_for_local_variance) <br />
[4 ](https://github.com/AHassanpour88/Beyond_Scenario_MoBPSopti/edit/main/README.md#3-Sensivitiy_Script) <br />

## 0 Introduction
This repository contains scripts for resource optimization in breeding program using Kernel regression. Breeding programs nowadays have become increasingly larger and more structurally complex with various interdependent parameters and contrasting objectives. Therefore, it is practically impossible to derive a strategy intuitively to optimize a breeding program for available resources. As a result, it is a common practice to narrow down the optimization problem to a set of scenarios and focus on a smaller subset of possibilities that, in turn, are deeply analyzed in detail. We aim to provide guidance on constructing a multi-objective optimization problem using stochastic simulations and integrating its use into breeding research programs seeking to design, implement and evaluate the best resource utilization or near-optimal resource combination beyond just analyzing scenario differences. 


## 1 Description of the underlying simulation case
We will propose a general pipeline for optimizing breeding programs in the following. We will provide an exemplary breeding scenario of a simplified classical dairy cattle breeding scheme. In dairy cattle breeding program, the performance traits of a bull cannot be determined phenotypically and can be expressed only by cows, so pre-selected bulls must be mated to cows to produce test daughters. The offspring performance of the test bulls will be used as criteria for selection decisions, and selected sires and cows with the desired trait will then be used as parents for the next cycle. In our example, a genomic breeding value estimation was used to select bulls, whereas cows were selected based on pedigree breeding value estimation. With the two older generations of cows having phenotypic information available, the last three generations have been considered for selection. Using 1,000 purely additive QTLs, a mean genomic value of 100, and a genomic standard deviation of 10, we simulated a single trait with a heritability of 0.3. [Simulation_Script](https://github.com/AHassanpour88/Beyond_Scenario_MoBPSopti/edit/main/README.md#1-Simulation_Script)


$4000x_1 + 3000x_2 - 10000000 &\leq 0\\$
$x_1 + x_2 + x_3 &\geq 0\\$
$100 &\leq x_2 \leq 700\\$
$3 &\leq x_3 \leq 30$


## 2 Nadaraya-Watson kernel function estimator


## 3 Estimation of the number of needed simulations


## 4 
