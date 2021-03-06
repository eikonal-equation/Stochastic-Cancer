/*=============================================================================
 * Copyright (C) 2022 MingYi Wang
 *
 * This program is free software: you can redistribute it and/or modify it under the terms
 * of the GNU General Public License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program.
 * If not, see http://www.gnu.org/licenses/.
 *============================================================================*/
 
 
/*==============================================================================
 * File: SemiLagrangian_Cancer.h
 *
 * Author: MingYi Wang
 *
 * Description: This file contains the declarations of functions that compute the 
 * drift (deterministic) portion of the cancer dynamics; the coefficients and 
 * feet of the 1st-order weak approximation of the cancer dynamics; and the main 
 * solver for the value function using semi-Lagrangian scheme.
 *
 * Details of all of these functions are found in SemiLagrangian_Cancer.h.cpp.
 *
 *============================================================================*/
#pragma once
#ifndef SEMILAGRANGIAN_CANCER_H
#define SEMILAGRANGIAN_CANCER_H

//----------------------Project specific header files---------------------------
#include "ENO3_2d.h"

//-------------------------Libraires-------------------------------------------
#include<algorithm>
#include<numeric>

//---------------------------Main Class----------------------------------------
class CancerSL
{
public:

	CancerSL() = default;
	CancerSL(int a_Factor, double a_Budget, double a_Treatment_const, double a_Diff_const);

    // This function find the first index s.t. xloc <= aArray
	int find_index(const ublas::vector<double>& aArray, const double xloc); 

	double drift_p(const double q, const double p, double aDose); // Drift function for q

	double drift_q(const double q, const double p); // Drift function for p

	// This function pre-compute all the coefficients and feet of the 1st-order weak approximation of the dynamics
	tuple<ublas::matrix<vector<double>>, ublas::matrix<vector<double>>, ublas::matrix<vector<double>>, ublas::matrix<vector<double>>,
		ublas::matrix<vector<int>>, ublas::matrix<vector<int>>, ublas::matrix<vector<int>>, ublas::matrix<vector<int>>>
		precompute_coeff_2D(const ublas::vector<double>& qArray, const ublas::vector<double>& pArray);


	// This is the main solver that uses time marching and semi-Lagrangian schemes. 
	// I.e., given a budget s, this function will compute the value function 
	// and return policies in feedbackform for each slice of cost up to the budget.
	ublas::matrix<double> MainSolver_by_SL();

private:
    // An input factor that would be multiplied by 100 used as number of points on each side of spatial grid (Base case: 100*2)
	int fMulti_factor; 
	
	double fBudget; // The cost we are integrating up to.
	int fN; // number of points on each side of spatial grid
	int fM; // number of slices
	double fx0; // starting point at x-axis
	double fy0; // starting point at y-axis
	double fDx; // Delta x(or q)
	double fDy; // Delta y(or p)
	double fDs; // Delta s (discretization of cost)

	//parameters from Gluzman et al. https://royalsocietypublishing.org/doi/10.1098/rspb.2019.2454
	double ba;
	double bv;
	double c;
	double n;
    // Diffusion constants (volatilities)
	double fDiff_const_1; 
	double fDiff_const_2;
	double fDiff_const_3;
	
	double fSigma; //treatment constant
	double fDmax; // MTD (here dmax = 3)
	int index_rec; // index of the recovery barrier on y(p)-axis
	int index_death; // index of death barrier on y(p)-axis
	double tau_mtd; // time step for MTD
	double tau_0; // time step for no drug
};
#endif // ! SEMILAGRANGIAN_CANCER_H

