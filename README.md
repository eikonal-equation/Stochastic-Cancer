# Stochastic-Cancer 
A (soon to be) public repository with codes &amp; movies from our paper on environmental stochasticity in cancer models.


This repository contains the source code used to generate all examples presented in "Stochastic optimal control to guide adaptive cancer therapy treatment" manuscript (both the main text and the Supplementary Materials) by MingYi Wang, Jacob G. Scott and Alexander Vladimirsky.

# License #
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.

## Abstract of our manuscript ##
While adaptive cancer therapy is beginning to prove a promising approach of building evolutionary dynamics into therapeutic scheduling, the stochastic nature of cancer evolution has rarely been incorporated. Various sources of random perturbations can impact the evolution of heterogeneous tumors.

We propose a method that can effectively select optimal adaptive treatment policies under randomly evolving tumor dynamics based on Stochastic Optimal Control theory. 
We first construct a stochastic model of cancer dynamics under drug therapy based on Evolutionary Game theory. That model is then used to improve the cumulative "cost", a combination of the total amount of drugs used and the time to recovery. As this cost becomes random in a stochastic setting, we maximize the probability of recovery under a pre-specified cost threshold (or a "budget"). 

We can achieve our goal for a range of threshold values simultaneously using the tools of dynamic programming. We then compare our threshold-aware policies with the policies previously shown to be optimal in the deterministic setting. We show that this threshold-awareness yields a significant improvement in the probability of under-the-budget recovery, which is correlated with a lower general drug usage.

# Manuscript #
The Stochastic Optimization of Adaptive Cancer Therapy manuscript can be found [here](https://github.com/eikonal-equation/Stochastic-Cancer).

# Contributions & Acknowledgements # 
  * The problem statement, the specific stochastic cancer dynamics, and the numerical scheme were developed by MingYi Wang, Jacob Scott, and Alexander Vladimirsky.
  * The implementation was carried out by MingYi Wang.
  * The manuscript was written by MingYi Wang, Jacob Scott, and Alexander Vladimirsky.
  * The authors also acknowledge Mark Gluzman, Artem Kaznatcheev, Robert Vander Velde, and David Basanta for their profound previous work.

# Instructions #
  
## Requirements: ## 
* The C++ code requires the users to install the "[Boost](https://www.boost.org/)" library (external). We install it in the same directory where `include` and `src` locate (feel free to install it anywhere you want on your end but you need to modify its path in the ***Makefile*** described below).

* The CDFs and the deterministic-optimal policy are generated with Matlab code.

## Running the C++ Code: ##
The following instructions explain how to run the Solver for threshold-aware optimal policy using the ***Makefile***. 

To change compilers, edit `CC=` by putting your desired compiler after the `=` sign. The default compiler is set to `g++`. 

To update the path/version of the "Boost" library you've installed, edit `BOOSTPATH` by putting your own path/version after the `=` sign. The current path/version is `./boost_1_75_0`.

To compile the code, type `make main` in the terminal at this folder. 

To delete the executable and the date files, type `make clean`.

## Running the Matlab Code: ##
To generate a CDF with threshold-aware optimal policies:
  * `CDF_threshold_aware_policy.m`
      * Produces a CDF y = Pr(J <= s) that measures the probability of keeping the accumulative cost J under a given initial threshold/budget value s using threshold-aware policies. It will produce a plot of the CDF and a plot of a sample path at the end of execution. 
      * Note: it requires a data matrix of the deterministic-optimal policy and a multi-dimensional data array of the threshold-aware policies (with the same same spatial dimensions) as inputs. 
          * E.g., a 1601x1601 data matrix for the deterministic-optimal policy and a 1601x1601x1201 multidimensional array for the threshold-aware policies.

To generate a CDF with the deterministic-optimal policy:
   * `CDF_stationary_policy.m`
      * Produces a CDF y = Pr(J <= s) that measures the probability of success where the cumulative cost J is within any positive threshold value s using the stationary policy. It will produce a plot of the CDF and a superposed plot of a sample path together with the optimal policy in the background at the end of execution.
      * Note: it requires a data matrix of deterministic-optimal policy as an input.

To generate the deterministic-optimal policy from [Gluzman et al.](https://royalsocietypublishing.org/doi/10.1098/rspb.2019.2454) (but with our own sets of parameter values):
  * `Deterministic_Cancer_ValuePolicy_Ite.m`
      * Produces the determinstic-optimal policy matrix named "Dmat" and its corresponding value function matrix named "U".
      * Note: it requies the dimension of the desired output matrix as an input.
