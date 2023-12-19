# Threshold-awareness in adaptive cancer therapy
This repository contains the source code used to generate all examples presented in "Threshold-awareness in adaptive cancer therapy" manuscript (both the main text and the Supplementary Materials) by MingYi Wang, Jacob G. Scott, and Alexander Vladimirsky.

# License #
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.

## Abstract of our manuscript ##
While adaptive cancer therapy is beginning to prove a promising approach of building evolutionary dynamics into therapeutic scheduling, the stochastic nature of cancer evolution has rarely been incorporated. Various sources of random perturbations can impact the evolution of heterogeneous tumors. In this paper, we propose a method that can effectively select optimal adaptive treatment policies under randomly evolving tumor dynamics based on Stochastic Optimal Control theory. 

We first construct a stochastic model of cancer dynamics under drug therapy based on Evolutionary Game theory. 
That model is then used to improve the cumulative "cost", a combination of the total amount of drugs used and the time to recovery. 
As this cost becomes random in a stochastic setting, we maximize the probability of recovery under a pre-specified cost threshold (or a "budget"). 
We can achieve our goal for a range of threshold values simultaneously using the tools of dynamic programming.

We then compare our threshold-aware policies with the policies previously shown to be optimal in the deterministic setting.
We show that this threshold-awareness yields a significant improvement in the probability of under-the-budget recovery, which is correlated with a lower general drug usage.

The particular model underlying our discussion has originated in [Kaznatcheev2017 et al.](https://www.nature.com/articles/bjc20175), but the presented approach is far more general and provides a new tool for optimizing adaptive therapies based on a broad range of stochastic cancer models.

# Manuscript #
The "Stochastic optimal control to guide adaptive cancer therapy" manuscript can be found [here](https://www.biorxiv.org/content/10.1101/2022.06.17.496649v2).

# Contributions & Acknowledgements # 
  * The problem statement, the specific stochastic cancer dynamics, and the numerical scheme were developed by MingYi Wang, Jacob Scott, and Alexander Vladimirsky.
  * The implementation was carried out by MingYi Wang.
  * The manuscript was written by MingYi Wang, Jacob Scott, and Alexander Vladimirsky.
  * The authors acknowledge Mark Gluzman, Artem Kaznatcheev, Robert Vander Velde, and David Basanta for inspiring our work. The authors are also grateful to Roberto Ferretti and Lars Grüne for their advice on some aspects of numerical methods used in this project.

# Instructions #
  
## Requirements: ## 
* The C++ code requires the users to install the "[Boost](https://www.boost.org/)" library (external). 
    * We install it in the same directory where `include` and `src` are located. Please feel free to install it anywhere you want on your end but you would have to modify its path in the ***Makefile*** described below.

* The CDFs and the deterministic-optimal policy are generated with Matlab code.

## Running the C++ Code: ##
The following instructions explain how to run the Solver for our threshold-aware optimal value functions/policies using the ***Makefile***. 

To change compilers, edit `CC=` by putting your desired compiler after the `=` sign. The default compiler is set to `g++`. 

To update the path/version of the "Boost" library you've installed, edit `BOOSTPATH` by putting your own path/version after the `=` sign. The current path/version is `./boost_1_79_0`.

To compile the code, type `make main` in the terminal in this folder. 
* Our C++ code is capable of handling both Examples (models) presented in the main text. To switch between models, modify the global variable `gCase` in the `main.cpp` to either **"Example1"** or **"Example2"**, as appropriate.

To delete the executable and the data files, type `make clean`.

## Running the Matlab Code: ##
To generate a CDF with threshold-aware optimal policies (Please go to the respective folder (`Example1` / `Example2`) to run the code for the model you would like to test):
  * `CDF_threshold_aware_policy_withDetOpt.m` OR `CDF_threshold_aware_policy_EG2.m`
      * Produces a CDF $y(s) = {\mathbb{P}}(J \le s)$ that measures the probability of keeping the accumulative cost $J$ under any threshold value $s$ but maximized at a given initial threshold/budget $s_0$ using threshold-aware policies. It will produce a plot of the CDF and a plot of a sample path at the end of execution. 
      * Note: both of them require a data matrix of the deterministic-optimal policy and a data file path of threshold-aware policies (with the same spatial dimensions) as inputs. 
          * E.g., a 1601x1601 data matrix for the deterministic-optimal policy and a data file path containing a 1601x1601x1201 multi-dimensional array for the threshold-aware policies.
      * We also provide `CDF_threshold_aware_policy_withMTD.m` to generate CDFs with the MTD-based therapy used after the budget runs out for Example 1. This version does not require a data matrix of the deterministic-optimal policy as an input. 

To generate a CDF with the deterministic-optimal policy:
   * `CDF_stationary_policy.m` OR `CDF_stationary_policy_EG2.m`
      * Produces a CDF $y(s) = {\mathbb{P}}(J \le s)$ that measures the probability of success where the cumulative cost $J$ is within any positive threshold value $s$ using the stationary policy. It will produce a plot of the CDF and a superposed plot of a sample path together with the optimal policy in the background at the end of execution.
      * Note: it requires a data matrix of the deterministic-optimal policy as an input.

To generate the deterministic-optimal policy for Model 1 (based on [(Gluzman et al. 2020)](https://royalsocietypublishing.org/doi/10.1098/rspb.2019.2454)):
  * `Deterministic_Cancer_ValuePolicy_Ite_EG1.m` (under the folder `Example1`)
      * Produces the deterministic-optimal policy matrix named `Dmat` and its corresponding value function matrix named `U`.
      * Note: it requires the dimension of the desired output matrix `N` as an input.
 To generate the deterministic-optimal policy for Model 2 (based on [(Carrère 2017)](https://doi.org/10.1016/j.jtbi.2016.11.009)):
  * `Deterministic_Cancer_ValuePolicy_Ite_EG2.m` (under the folder `Example2`)
      * Same input/output format as above
    

Demonstration:
  * `demo_EG1.m` under the folder `Example1`
      * This script is solely for demonstrating Model 1.
      * The demonstration is set to use a coarse spatial grid (401x401). Hence, it will not produce a desirable accuracy. To reproduce the results shown in our paper, please refer to section 4S of our Supplementary Materials for implementation details.
  * `demo_EG2.m` under the folder `Example2`
      * This script is solely for demonstrating Model 2.
 
