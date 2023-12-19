/*=============================================================================
 * Copyright (C) 2023 MingYi Wang
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
  * File: Solver_SL.h
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
#ifndef SOLVER_SL_H
#define SOLVER_SL_H

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
	CancerSL(int a_Factor, double a_Slice_factor, double a_Budget, double a_Treatment_const, double a_Diff_const, int a_Num_diffpts_interior, int a_Num_diffpts_boundary,
		double a_Thres_recovery, double a_Thres_death, double a_MaxDose, string a_Case);
	CancerSL(int a_Factor, double a_Slice_factor, double a_Budget, double a_Treatment_const, double a_Diff_const, int a_Num_diffpts_interior, int a_Num_diffpts_boundary,
		double a_Thres_recovery, double a_Thres_death, double a_MaxDose, double a_Alpha, double a_Beta, double a_Carrying_capacity, double a_Time_rescale, string a_Case);

	// This function find the first index s.t. xloc <= aArray
	int find_index(const ublas::vector<double>& aArray, const double xloc);

	double drift_p(const double q, const double p, double aDose, string aCase); // Drift function for q

	double drift_q(const double q, const double p, double aDose, string aCase); // Drift function for p

	array<double,3> diffusion_q_3d(const double q, const double p); // Diffusion function for q with a 3D BM

	array<double, 3> diffusion_p_3d(const double q, const double p); // Diffusion function for p with a 3D BM

	array<double, 2> diffusion_q_2d(const double q, const double p); // Diffusion function for q with a 2D BM

	array<double, 2> diffusion_p_2d(const double q, const double p); // Diffusion function for p with a 2D BM

	double diffusion_q_1d(const double q, const double p); // Diffusion function for q with a 1D BM

	double diffusion_p_1d(const double q, const double p); // Diffusion function for p with a 1D BM

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
	double fSlice_factor;
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
	//parameters from Melbinger et al. https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.105.178101
	double gs;
	double gr;
	double m_size;
	double Kc;
	double alpha;
	double beta;
	double beta_hat;
	// Diffusion constants (volatilities)
	double fDiff_const_1;
	double fDiff_const_2;
	double fDiff_const_3;

	double fDelta; //treatment constant
	double fDmax; // MTD (here dmax = 3)
	double fThres_rec; // threshold of recovery
	double fThres_death; // threshold of failure
	int index_rec; // index of the recovery barrier on y(p)-axis
	int index_start_q; // staring index on the q-axis
	int index_death; // index of death barrier on y(p)-axis
	double tau_mtd; // time step for MTD
	double tau_0; // time step for no drug
	int fNum_diffpts_interior; // number of sample points for the diffusion process in the interior of the domain
	int fNum_diffpts_boundary; // number of sample points for the diffusion process on the boundary of the domain
	string fCase; //Indication of which example to run
	double fTol;
	int fStorage_factor;
};
#endif // ! SOLVER_SL_H

