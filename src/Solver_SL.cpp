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
  * File: Solver_SL.cpp
  *
  * Author: MingYi Wang
  *
  * Description: This file contains the implementation of functions that compute the
  * drift (deterministic) portion of the cancer dynamics; the coefficients and
  * feet of the 1st-order weak approximation of the cancer dynamics; and the main
  * solver for the value function using semi-Lagrangian scheme.
  *
  *============================================================================*/
#include "Solver_SL.h"

  //----------------------Project specific header files---------------------------
#include "WriteToFile.h"

//------------------------------Libraries------------------------------------
#include<chrono>
#include<omp.h>
#include<iomanip>

// Define the constructor for Example 1
CancerSL::CancerSL(int a_Factor, double a_Slice_factor, double a_Budget, double a_Treatment_const, double a_Diff_const, int a_Num_diffpts_interior, int a_Num_diffpts_boundary,
	double a_Thres_recovery, double a_Thres_death, double a_MaxDose, string a_Case)
{
	// initialization of the 2D grid
	fMulti_factor = a_Factor;
	fSlice_factor = a_Slice_factor;
	fBudget = a_Budget;
	fN = 100 * fMulti_factor;

	fM = int(fSlice_factor * fMulti_factor * fBudget); //for octed ds


	fDs = fBudget / fM;
	fDx = 1 / double(fN);
	fDy = fDx;
	fx0 = 0;
	fy0 = 0;

	// parameter values from Gluzman et al. https://royalsocietypublishing.org/doi/10.1098/rspb.2019.2454
	ba = 2.5;
	bv = 2;
	c = 1;
	n = 4;
	fDmax = a_MaxDose; //(dmax = 3)

	// Diffusion constants (volatilities) and Treatment constants
	//we have only considered cases that all three diffusion constants being the same
	fDiff_const_1 = a_Diff_const;
	fDiff_const_2 = fDiff_const_1;
	fDiff_const_3 = fDiff_const_1;

	fDelta = a_Treatment_const;
	fThres_rec = a_Thres_recovery; // threshold of stabilization
	fThres_death = a_Thres_death; // threshold of failure
	index_rec = int(fThres_rec / fDx); // index of the stabilization barrier
	index_death = fN - index_rec; // index of the death barrier

	tau_mtd = fDs / (fDmax + fDelta); // time-step for MTD-based therapy
	tau_0 = fDs / fDelta; // time-step for no therapy

	fNum_diffpts_interior = a_Num_diffpts_interior; // number of diffusion samples in the interior
	fNum_diffpts_boundary = a_Num_diffpts_boundary; // number of diffusion samples on the boundary
	fCase = a_Case;
    fTol = 1e-14;
    fStorage_factor = 5; // store every 5 slices
	index_start_q = 0; // start with q = 0
}


// Define the overloaded constructor for Example 2
CancerSL::CancerSL(int a_Factor, double a_Slice_factor, double a_Budget, double a_Treatment_const, double a_Diff_const, int a_Num_diffpts_interior, int a_Num_diffpts_boundary,
	double a_Thres_recovery, double a_Thres_death, double a_MaxDose, double a_Alpha, double a_Beta, double a_Carrying_capacity, double a_Time_rescale, string a_Case)
{
	// initialization of the 2D grid
	fMulti_factor = a_Factor;
	fSlice_factor = a_Slice_factor;
	fBudget = a_Budget;
	fN = 100 * fMulti_factor;

	fM = int(fSlice_factor * fMulti_factor * fBudget);


	fDs = fBudget / fM;
	fDx = 1 / double(fN);
	fDy = fDx;
	fx0 = 0;
	fy0 = 0;

	//parameters from Carrère 2017 https://doi.org/10.1016/j.jtbi.2016.11.009
	gs = 0.031*a_Time_rescale;
	gr = 0.026*a_Time_rescale;
	m_size = 30;
	Kc = a_Carrying_capacity;
	alpha = a_Alpha;
	beta = a_Beta;
	beta_hat = beta * Kc;
	fDmax = a_MaxDose; //(dmax = 3)

	// Diffusion constants (volatilities) and Treatment constants
	//we have only considered cases that all three diffusion constants being the same
	fDiff_const_1 = a_Diff_const*sqrt(a_Time_rescale);
	fDiff_const_2 = fDiff_const_1;
	fDiff_const_3 = fDiff_const_1;

	fDelta = a_Treatment_const;
	fThres_rec = a_Thres_recovery; // threshold of remission
	fThres_death = a_Thres_death; // threshold of failure
	index_rec = int(fThres_rec / fDx); // index of the remission barrier
	index_death = fN - index_rec; // index of the death barrier

	tau_mtd = fDs / (fDmax + fDelta); // time-step for MTD-based therapy
	tau_0 = fDs / fDelta; // time-step for no therapy

	fNum_diffpts_interior = a_Num_diffpts_interior; // number of diffusion samples in the interior
	fNum_diffpts_boundary = a_Num_diffpts_boundary; // number of diffusion samples on the boundary
	fCase = a_Case;
    fTol = 1e-14;
    fStorage_factor = 10; // store every 10 slices
	index_start_q = 0; // start with q = 0
}

// Define find_index function
// aArray (input): a (N+1)x1 array of sample points
// xloc (input): coordinate of the query point
//
// index (output): the first index s.t. xloc <= aArray
int CancerSL::find_index(const ublas::vector<double>& aArray, const double xloc)
{
	char flag = 'F';
	int index = 0;
	for (int i = 0; i < fN + 1; i++)
	{
		if (xloc <= aArray[i])
		{
			flag = 'T';
			index = i;
			break;
		}
	}
	return index;
}

// Define drift functions
// q (input): x-coordinate of the query point
// p (input): y-coordinate of the query point
// aDose(input): dosage of drugs (either 0 or dmax=3 in our case)
inline double CancerSL::drift_p(const double q, const double p, const double aDose, const string aCase)
{
    double returnVal = -1;
	if (aCase == "Example1") {
		returnVal = p * (1 - p) * (ba / (n + 1) - q * (bv - c) - aDose) -
			p * (pow(fDiff_const_1, 2) * p - (pow(fDiff_const_1, 2) * pow(p, 2)
				+ pow(fDiff_const_2, 2) * pow(1 - p, 2) * pow(1 - q, 2) + pow(fDiff_const_3, 2) * pow(1 - p, 2) * pow(q, 2)));
	}
	else if (aCase == "Example2") {
		returnVal = p * (1 - p) * (gs * q + gr * (1 - q)) - alpha * q * p * aDose - beta_hat * pow(p, 2) * q * (1 - q);
	}
	return returnVal;
}

inline double CancerSL::drift_q(const double q, const double p, const double aDose, const string aCase)
{
    double returnVal = -1;
	if (aCase == "Example1") {
		returnVal = q * (1 - q) * (bv / (n + 1) * (1 + p + pow(p, 2) + pow(p, 3) + pow(p, 4)) - c) +
			q * (1 - q) * ((1 - q) * pow(fDiff_const_2, 2) - q * pow(fDiff_const_3, 2));
	}
	else if (aCase == "Example2") {
		//The case of 1d BM
		returnVal = (1 - p) * q * (1 - q) * (gs - gr) - alpha * q * (1 - q) * aDose + beta_hat * p * pow(q, 2) * (1 - q)
			+ pow(1 - p, 2) * q * (1 - q) * (pow(fDiff_const_2, 2) * (1 - q) - pow(fDiff_const_1, 2) * q + fDiff_const_1 * fDiff_const_2);
	}
	return returnVal;
}

// Define diffusion functions
// q (input): x-coordinate of the query point
// p (input): y-coordinate of the query point
inline array<double,3> CancerSL::diffusion_p_3d(const double q, const double p)
{
	array<double, 3> diff_arr;
	diff_arr[0] = p * (1 - p); // for dW1
	diff_arr[1] = p * (1 - p) * (1 - q); // for dW2
	diff_arr[2] = p * (1 - p) * q; // for dW3
	return diff_arr;
}

inline array<double, 3> CancerSL::diffusion_q_3d(const double q, const double p)
{
	array<double, 3> diff_arr;
	diff_arr[0] = 0; // for dW1
	diff_arr[1] = -q * (1 - q); // for dW2
	diff_arr[2] = q * (1 - q); // for dW3
	return diff_arr;
}


inline array<double, 2> CancerSL::diffusion_p_2d(const double q, const double p)
{
	array<double, 2> diff_arr;
	diff_arr[0] = p * (1 - p) * q;// for dW1
	diff_arr[1] = p * (1 - p) * (1 - q); // for dW2
	return diff_arr;
}

inline array<double, 2> CancerSL::diffusion_q_2d(const double q, const double p)
{
	array<double, 2> diff_arr;
	diff_arr[0] = (1 - p) * q * (1 - q); // for dW1
	diff_arr[1] = -(1 - p) * q * (1 - q); // for dW2
	return diff_arr;
}


inline double CancerSL::diffusion_p_1d(const double q, const double p)
{
	return p * (1 - p) * (fDiff_const_1 * q + fDiff_const_2 * (1 - q));
}

inline double CancerSL::diffusion_q_1d(const double q, const double p)
{
	double returnVal;
	if (fDiff_const_1 == fDiff_const_2) {
		returnVal = 0;
	}
	else
	{
		returnVal = q * (1 - q) * (1 - p) * (fDiff_const_1 - fDiff_const_2);
	}
	return returnVal;
}


// This function pre-compute all the coefficients and feet of the 1st-order weak approximation of the dynamics
void CancerSL::precompute_coeff_2D(const ublas::vector<double>& q, const ublas::vector<double>& p, ublas::matrix<vector<double>>& sample_x_nodrug,
	ublas::matrix<vector<double>>& sample_x_mtd, ublas::matrix<vector<double>>& sample_y_nodrug, ublas::matrix<vector<double>>& sample_y_mtd,
	ublas::matrix<vector<int>>& index_x_nodrug, ublas::matrix<vector<int>>& index_x_mtd,
	ublas::matrix<vector<int>>& index_y_nodrug, ublas::matrix<vector<int>>& index_y_mtd)
{
	//ublas::matrix<std::vector<double>> sample_x_nodrug(fN + 1, fN + 1); // store all possible locations in x-direction for no therapy
	//ublas::matrix<std::vector<double>> sample_x_mtd(fN + 1, fN + 1); // store all possible locations in x-direction for mtd
	//ublas::matrix<std::vector<double>> sample_y_nodrug(fN + 1, fN + 1); // store all possible locations in y-direction for no therapy
	//ublas::matrix<std::vector<double>> sample_y_mtd(fN + 1, fN + 1); // tore all possible locations in y-direction for mtd
	//ublas::matrix<std::vector<int>> index_x_nodrug(fN + 1, fN + 1); // store all possible coefficients at x-direction for no therapy
	//ublas::matrix<std::vector<int>> index_x_mtd(fN + 1, fN + 1); // store all possible coefficients at x-direction for mtd
	//ublas::matrix<std::vector<int>> index_y_nodrug(fN + 1, fN + 1); // store all possible coefficients at y-direction for no therapy
	//ublas::matrix<std::vector<int>> index_y_mtd(fN + 1, fN + 1); // store all possible coefficients at y-direction for mtd

	// square root of time steps
	double root_tau_0 = sqrt(tau_0);
	double root_tau_mtd = sqrt(tau_mtd);
	array<double, 3> diff_y_3d{0,0,0};
	array<double, 2> diff_y_2d{0,0};
	array<double, 3> diff_x_3d{0,0,0};
	array<double, 2> diff_x_2d{0,0};
	double diff_x_1d = 0;
	double diff_y_1d = 0;

	for (int i = index_rec; i < index_death + 1; i++)
	{
		for (int j = index_start_q; j < fN + 1; j++)
		{
			if (j == 0) // at q = 0
			{
				int k = 0;
				vector<int> ky_list{ 0,0,0,0 };
				vector<double> hat_y{ 0,0,0,0 };
				if (fCase == "Example1") {
					diff_y_3d = diffusion_p_3d(0, p(i)); // On q = 0, there is only diffusion in p direction
				}
				else if (fCase == "Example2") {
					diff_y_1d = diffusion_p_1d(0, p(i));
				}
				// for No therapy, i.e. Dose = 0
				double y_no_drug = p(i) + tau_0 * drift_p(0, p(i), 0, fCase); // foot of characteristics at (0,p(i)) for no therapy


				for (int k = 0; k < fNum_diffpts_boundary; k++) {
					if (fCase == "Example1") {
						hat_y[k] = y_no_drug + (1 - 2 * (k % 2)) * root_tau_0 * diff_y_3d[0] * fDiff_const_1 + (1 - 2 * (k / 2 % 2)) * root_tau_0 * diff_y_3d[1] * fDiff_const_2;
					}
					else if (fCase == "Example2") {
						hat_y[k] = y_no_drug + (1 - 2 * (k % 2)) * root_tau_0 * diff_y_1d;
					}
				}

				for (int m = 0; m < fNum_diffpts_boundary; m++)
				{
					k = find_index(p, hat_y[m]);
					ky_list[m] = k;
				}
				index_y_nodrug(i, j) = ky_list;
				sample_y_nodrug(i, j) = hat_y;

				// for MTD-based therapy, i.e. Dose = fDmax
				double y_mtd = p(i) + tau_mtd * drift_p(0, p(i), fDmax, fCase); // foot of characteristics at (0,p(i)) for mtd

				hat_y = { 0,0,0,0 };

				for (int k = 0; k < fNum_diffpts_boundary; k++) {
					if (fCase == "Example1") {
						hat_y[k] = y_mtd + (1 - 2 * (k % 2)) * root_tau_mtd * diff_y_3d[0] * fDiff_const_1 + (1 - 2 * (k / 2 % 2)) * root_tau_mtd * diff_y_3d[1] * fDiff_const_2;
					}
					else if (fCase == "Example2") {
						hat_y[k] = y_mtd + (1 - 2 * (k % 2)) * root_tau_mtd * diff_y_1d;
					}
				}

				ky_list = { 0,0,0,0 };

				for (int m = 0; m < fNum_diffpts_boundary; m++)
				{
					k = find_index(p, hat_y[m]);
					ky_list[m] = k;
				}
				index_y_mtd(i, j) = ky_list;
				sample_y_mtd(i, j) = hat_y;

			}
			else if (j == fN) // at q = 1
			{
				vector<int> ky_list{ 0,0,0,0 };
				int k = 0;
				vector<double> hat_y{ 0,0,0,0 };
				if (fCase == "Example1") {
					diff_y_3d = diffusion_p_3d(1, p(i)); // On q = 1, there is only diffusion in p direction
				}
				else if (fCase == "Example2") {
					//diff_y_2d = diffusion_p_2d(1, p(i));
					diff_y_1d = diffusion_p_1d(1, p(i));
				}
				// for No therapy, i.e. Dose = 0
				double y_no_drug = p(i) + tau_0 * drift_p(1, p(i), 0, fCase); // foot of characteristics at (1,p(i)) for no drug

				for (int k = 0; k < fNum_diffpts_boundary; k++) {
					if (fCase == "Example1") {
						hat_y[k] = y_no_drug + (1 - 2 * (k % 2)) * root_tau_0 * diff_y_3d[0] * fDiff_const_1 + (1 - 2 * (k / 2 % 2)) * root_tau_0 * diff_y_3d[2] * fDiff_const_3;
					}
					else if (fCase == "Example2") {
						hat_y[k] = y_no_drug + (1 - 2 * (k % 2)) * root_tau_0 * diff_y_1d;
					}
				}

				for (int m = 0; m < fNum_diffpts_boundary; m++)
				{
					k = find_index(p, hat_y[m]);
					ky_list[m] = k;
				}
				index_y_nodrug(i, j) = ky_list;
				sample_y_nodrug(i, j) = hat_y;

				// for MTD-based therapy, i.e. Dose = fDmax
				double y_mtd = p(i) + tau_mtd * drift_p(1, p(i), fDmax, fCase); // foot of characteristics at (1,p(i)) for mtd

				hat_y = { 0,0,0,0 };

				for (int k = 0; k < fNum_diffpts_boundary; k++) {
					if (fCase == "Example1") {
						hat_y[k] = y_mtd + (1 - 2 * (k % 2)) * root_tau_mtd * diff_y_3d[0] * fDiff_const_1 + (1 - 2 * (k / 2 % 2)) * root_tau_mtd * diff_y_3d[2] * fDiff_const_3;
					}
					else if (fCase == "Example2") {
						hat_y[k] = y_mtd + (1 - 2 * (k % 2)) * root_tau_mtd * diff_y_1d;
					}
				}

				ky_list = { 0,0,0,0 };

				for (int m = 0; m < fNum_diffpts_boundary; m++)
				{
					k = find_index(p, hat_y[m]);
					ky_list[m] = k;
				}
				index_y_mtd(i, j) = ky_list;
				sample_y_mtd(i, j) = hat_y;
			}
			else
			{
				vector<int> kx_list{ 0,0,0,0,0,0,0,0 };
				vector<int> ky_list{ 0,0,0,0,0,0,0,0 };
				int kx = 0;
				int ky = 0;
				vector<double> hat_x{ 0,0,0,0,0,0,0,0 };
				vector<double> hat_y{ 0,0,0,0,0,0,0,0 };
				if (fCase == "Example1") {
					diff_y_3d = diffusion_p_3d(q(j), p(i));
					diff_x_3d = diffusion_q_3d(q(j), p(i));
				}
				else if (fCase == "Example2") {
					//diff_y_2d = diffusion_p_2d(q(j), p(i));
					//diff_x_2d = diffusion_q_2d(q(j), p(i));
					diff_y_1d = diffusion_p_1d(q(j), p(i));
					diff_x_1d = diffusion_q_1d(q(j), p(i));
				}


				// for No therapy, i.e. D = 0
				double x_no_drug = q(j) + tau_0 * drift_q(q(j), p(i), 0, fCase); // x-coordinate of foot of characteristics at (q(j),p(i))
				double y_no_drug = p(i) + tau_0 * drift_p(q(j), p(i), 0, fCase); // y-coordinate of foot of characteristics at (q(j),p(i))

				//k / 4 % 2 is checking the operation status for W_1;
				//k % 2 is checking the operation status for W_2;
				//k / 2 % 2 is checking the operation status for W_3;
				for (int k = 0; k < fNum_diffpts_interior; k++) {
					if (fCase == "Example1") {
						hat_x[k] = x_no_drug + (1 - 2 * (k / 2 % 2)) * root_tau_0 * diff_x_3d[2] * fDiff_const_3 + (1 - 2 * (k % 2)) * root_tau_0 * diff_x_3d[1] * fDiff_const_2;
						hat_y[k] = y_no_drug + (1 - 2 * (k / 4 % 2)) * root_tau_0 * diff_y_3d[0] * fDiff_const_1 + (1 - 2 * (k % 2)) * root_tau_0 * diff_y_3d[1] * fDiff_const_2
							+ (1 - 2 * (k / 2 % 2)) * root_tau_0 * diff_y_3d[2] * fDiff_const_3;
					}
					else if (fCase == "Example2") {
						hat_x[k] = x_no_drug + (1 - 2 * (k % 2)) * root_tau_0 * diff_x_1d;
						hat_y[k] = y_no_drug + (1 - 2 * (k % 2)) * root_tau_0 * diff_y_1d;
					}
				}

				for (int m = 0; m < fNum_diffpts_interior; m++)
				{
					kx = find_index(q, hat_x[m]);
					ky = find_index(p, hat_y[m]);
					kx_list[m] = kx;
					ky_list[m] = ky;
				}
				index_x_nodrug(i, j) = kx_list;
				index_y_nodrug(i, j) = ky_list;
				sample_x_nodrug(i, j) = hat_x;
				sample_y_nodrug(i, j) = hat_y;

				// for MTD-based therapy, i.e. D = fDmax
				hat_x = { 0,0,0,0,0,0,0,0 };
				hat_y = { 0,0,0,0,0,0,0,0 };
				kx_list = { 0,0,0,0,0,0,0,0 };
				ky_list = { 0,0,0,0,0,0,0,0 };
				kx = 0;
				ky = 0;

				double x_mtd = q(j) + tau_mtd * drift_q(q(j), p(i), fDmax, fCase); // x-coordinate of foot of characteristics at (q(j),p(i))
				double y_mtd = p(i) + tau_mtd * drift_p(q(j), p(i), fDmax, fCase); // y-coordinate of characteristics at (q(j),p(i))


				//k / 4 % 2 is checking the operation status for W_1;
				//k % 2 is checking the operation status for W_2;
				//k / 2 % 2 is checking the operation status for W_3;
				for (int k = 0; k < fNum_diffpts_interior; k++) {
					if (fCase == "Example1") {
						hat_x[k] = x_mtd + (1 - 2 * (k / 2 % 2)) * root_tau_mtd * diff_x_3d[2] * fDiff_const_3 + (1 - 2 * (k % 2)) * root_tau_mtd * diff_x_3d[1] * fDiff_const_2;
						hat_y[k] = y_mtd + (1 - 2 * (k / 4 % 2)) * root_tau_mtd * diff_y_3d[0] * fDiff_const_1 + (1 - 2 * (k % 2)) * root_tau_mtd * diff_y_3d[1] * fDiff_const_2
							+ (1 - 2 * (k / 2 % 2)) * root_tau_mtd * diff_y_3d[2] * fDiff_const_3;
					}
					else if (fCase == "Example2") {
						hat_x[k] = x_mtd + (1 - 2 * (k % 2)) * root_tau_mtd * diff_x_1d;
						hat_y[k] = y_mtd + (1 - 2 * (k % 2)) * root_tau_mtd * diff_y_1d;
					}
				}

				for (int m = 0; m < fNum_diffpts_interior; m++)
				{
					kx = find_index(q, hat_x[m]);
					ky = find_index(p, hat_y[m]);
					kx_list[m] = kx;
					ky_list[m] = ky;
				}
				index_x_mtd(i, j) = kx_list;
				index_y_mtd(i, j) = ky_list;
				sample_x_mtd(i, j) = hat_x;
				sample_y_mtd(i, j) = hat_y;
			}
		}
	}
	//return make_tuple(sample_x_nodrug, sample_x_mtd, sample_y_nodrug, sample_y_mtd, index_x_nodrug, index_x_mtd, index_y_nodrug, index_y_mtd);
}


//--------------------------------Main Solver-----------------------------------
ublas::matrix<double> CancerSL::MainSolver_by_SL()
{
	//int downsample_size = 1;
	//-------------------------------------Initialization---------------------------------

	// initialize the matrix of the value function at s=0
	ublas::matrix<double> Vmat_old(fN + 1, fN + 1);
	for (int i = 0; i < fN + 1; i++)
	{
		for (int j = 0; j < fN + 1; j++)
		{
			if ((i >= 0) && (i < index_rec))  // boundary conditions
			{
				Vmat_old(i, j) = 1;
			}
			else
			{
				Vmat_old(i, j) = 0;
			}
		}
	}

	// initialize the matrix of the policy at s=0
	ublas::matrix<bool> policy_mat(fN + 1, fN + 1);
	for (int i = 0; i < fN + 1; i++)
	{
		for (int j = 0; j < fN + 1; j++)
		{
			policy_mat(i, j) = 0;
		}
	}
	// construct a uniform grid on [0,1]x[0,1]
	ENO2D myENOgrid2D(fN, fN, fDx, fDy, fx0, fy0);
	// construct a uniform grid on [0,1]
	ENO1D myENOgrid1D_y(fDy, fy0, fN);

	// save the initial value function and policy matrices
	 //string filename1 = "strict corrected ENO3_v2_2D N=1600 s=6 s1=0.5.dat";
	 //io::writeToFile2D<double>(filename1,Vmat_old);
	 //string filename2 = "strict corrected dval 2D N=1600 s=6 s1=0.5 tol=1e-14.dat";
	 //io::writeToFile2D<bool>(filename2,policy_mat);

	string filename1 = "test_valuefn.dat";
	//io::writeToFile2D_downsample<double>(filename1, Vmat_old, downsample_size);
	io::writeToFile2D<double>(filename1, Vmat_old);
	string filename2 = "test_policy.dat";
	io::writeToFile2D<bool>(filename2, policy_mat);


	// precompute all the coefficients
	// initialization of storage
	ublas::matrix<std::vector<double>> sample_x_nodrug(fN + 1, fN + 1); // store all possible locations in x-direction for no therapy
	ublas::matrix<std::vector<double>> sample_x_mtd(fN + 1, fN + 1); // store all possible locations in x-direction for mtd
	ublas::matrix<std::vector<double>> sample_y_nodrug(fN + 1, fN + 1); // store all possible locations in y-direction for no therapy
	ublas::matrix<std::vector<double>> sample_y_mtd(fN + 1, fN + 1); // tore all possible locations in y-direction for mtd
	ublas::matrix<std::vector<int>> index_x_nodrug(fN + 1, fN + 1); // store all possible coefficients at x-direction for no therapy
	ublas::matrix<std::vector<int>> index_x_mtd(fN + 1, fN + 1); // store all possible coefficients at x-direction for mtd
	ublas::matrix<std::vector<int>> index_y_nodrug(fN + 1, fN + 1); // store all possible coefficients at y-direction for no therapy
	ublas::matrix<std::vector<int>> index_y_mtd(fN + 1, fN + 1); // store all possible coefficients at y-direction for mtd


	// construct sample arrays in the horizontal (q) and vertical (p) direction respectively
	ublas::vector<double> q(fN + 1);
	ublas::vector<double> p(fN + 1);
	for (int i = 0; i < fN + 1; i++)
	{
		q(i) = fx0 + i * fDx;
		p(i) = fy0 + i * fDy;
	}

	// start time of computing all coefficients
	auto start = std::chrono::steady_clock::now();

	//pre-compute the coefficients
	precompute_coeff_2D(q, p, sample_x_nodrug, sample_x_mtd, sample_y_nodrug, sample_y_mtd, index_x_nodrug, index_x_mtd, index_y_nodrug, index_y_mtd);

	// end time of computating all coefficients
	auto end = std::chrono::steady_clock::now();

	// print out the computational time used for coefficients
	std::chrono::duration<double> elapsed_seconds = end - start;
	std::cout << "Elapsed time for precomputing all coefficients: " << elapsed_seconds.count() << "s\n";
	//-------------------------------------End of Initialization---------------------------------


	 //-------------------------------------Main Loop--------------------------------------------
	for (int l = 1; l < fM + 1; l++) //outermost loop for s-slices
	{
		auto start = std::chrono::steady_clock::now();

		// initialize the matrix of value function of the current slice at the beginning of each itration
		ublas::matrix<double> Vmat_current(fN + 1, fN + 1);
		for (int i = 0; i < fN + 1; i++)
		{
			for (int j = 0; j < fN + 1; j++)
			{
				if ((i >= 0) && (i < index_rec)) // boundary conditions
				{
					Vmat_current(i, j) = 1;
				}
				else
				{
					Vmat_current(i, j) = 0;
				}
			}
		}

		// initialize the matrix of policy of the current slice at the beginning of each itration
		ublas::matrix<bool> policy_mat(fN + 1, fN + 1);
		for (int i = 0; i < fN + 1; i++)
		{
			for (int j = 0; j < fN + 1; j++)
			{
				policy_mat(i, j) = 0;
			}
		}

		//-------------------Inner loops over the 2D spatial grid-------------------------------------

		#pragma omp parallel for
		for (int i = index_rec; i < index_death + 1; i++)
		{
			for (int j = index_start_q; j < fN + 1; j++)
			{
				if (j == 0) // on boundary q = 0
				{
					// For No therapy i.e. D = 0
					vector<double> hat_v{ 0,0,0,0 };
					vector<double> hat_y = sample_y_nodrug(i, j);
					vector<int> ky_list = index_y_nodrug(i, j);

					for (int m = 0; m < fNum_diffpts_boundary; m++)
					{
						if (hat_y[m] < fThres_rec)
						{
							hat_v[m] = 1;
						}
						else if (hat_y[m] > fThres_death)
						{
							hat_v[m] = 0;
						}
						else
						{
							array<double, 6> first_col;
							for (int i = 0; i < 6; i++)
							{
								first_col[i] = Vmat_old(i + ky_list[m] - 3, 0);
							}
							// 4th-order ENO cubic interpolation in 1D
							hat_v[m] = myENOgrid1D_y.ENO3_interp_1d(first_col, ky_list[m], hat_y[m]);
						}
					}

					double v1 = std::accumulate(hat_v.begin(), hat_v.end(), 0.0);

					v1 = v1 / double(fNum_diffpts_boundary); //take the average
					v1 = min(v1, 1.0); //ensuring the largest possible value is v=1

					// For MTD-based therapy i.e. D = fDmax
					hat_v = { 0,0,0,0 };
					hat_y = sample_y_mtd(i, j);
					ky_list = index_y_mtd(i, j);

					for (int m = 0; m < fNum_diffpts_boundary; m++)
					{
						if (hat_y[m] < fThres_rec)
						{
							hat_v[m] = 1;
						}
						else if (hat_y[m] > fThres_death)
						{
							hat_v[m] = 0;
						}
						else
						{
							array<double, 6> first_col;
							for (int i = 0; i < 6; i++)
							{
								first_col[i] = Vmat_old(i + ky_list[m] - 3, 0);
							}
							// 4th-order ENO cubic interpolation in 1D
							hat_v[m] = myENOgrid1D_y.ENO3_interp_1d(first_col, ky_list[m], hat_y[m]);
						}
					}

					double v2 = std::accumulate(hat_v.begin(), hat_v.end(), 0.0);

					v2 = v2 / double(fNum_diffpts_boundary); //take the average
					v2 = min(v2, 1.0); //ensuring the largest possible value is v=1

					Vmat_current(i, j) = max(v1, v2); //take the larger value

					if ((v1 >= v2) || ((std::abs(v1 - v2) / min(v1, v2)) < fTol))
					{
						// if the difference is small, we don't use drugs
						policy_mat(i, j) = 0;
					}
					else
					{
						policy_mat(i, j) = 1;
					}

				}
				else if (j == fN) // on boundary q = 1
				{
					// For No therapy i.e. D = 0
					vector<double> hat_v{ 0,0,0,0 };
					vector<double> hat_y = sample_y_nodrug(i, j);
					vector<int> ky_list = index_y_nodrug(i, j);

					for (int m = 0; m < fNum_diffpts_boundary; m++)
					{
						if (hat_y[m] < fThres_rec)
						{
							hat_v[m] = 1;
						}
						else if (hat_y[m] > fThres_death)
						{
							hat_v[m] = 0;
						}
						else
						{
							array<double, 6> last_col;
							for (int i = 0; i < 6; i++)
							{
								last_col[i] = Vmat_old(i + ky_list[m] - 3, fN);
							}
							// 4th-order ENO cubic interpolation in 1D
							hat_v[m] = myENOgrid1D_y.ENO3_interp_1d(last_col, ky_list[m], hat_y[m]);
						}
					}
					double v1 = std::accumulate(hat_v.begin(), hat_v.end(), 0.0);

					v1 = v1 / double(fNum_diffpts_boundary); //take the average
					v1 = min(v1, 1.0); //ensuring the largest possible value is v=1

					// For MTD-based therapy i.e. D = fDmax
					hat_v = { 0,0,0,0 };
					hat_y = sample_y_mtd(i, j);
					ky_list = index_y_mtd(i, j);

					for (int m = 0; m < fNum_diffpts_boundary; m++)
					{
						if (hat_y[m] < fThres_rec)
						{
							hat_v[m] = 1;
						}
						else if (hat_y[m] > fThres_death)
						{
							hat_v[m] = 0;
						}
						else
						{

							array<double, 6> last_col;
							for (int i = 0; i < 6; i++)
							{
								last_col[i] = Vmat_old(i + ky_list[m] - 3, fN);
							}
							// 4th-order ENO cubic interpolation in 1D
							hat_v[m] = myENOgrid1D_y.ENO3_interp_1d(last_col, ky_list[m], hat_y[m]);


						}
					}
					double v2 = std::accumulate(hat_v.begin(), hat_v.end(), 0.0);

					v2 = v2 / double(fNum_diffpts_boundary); // take the average
					v2 = min(v2, 1.0); //ensuring the largest possible value is v=1

					Vmat_current(i, j) = max(v1, v2); //take the larger value

					if ((v1 >= v2) || ((std::abs(v1 - v2) / min(v1, v2)) < fTol))
					{
						// if the difference is small, we don't use drugs
						policy_mat(i, j) = 0;
					}
					else
					{
						policy_mat(i, j) = 1;
					}
				}
				else // in the interior domain
				{
					// For No Therapy i.e. D = 0;
					vector<double> hat_v{ 0,0,0,0,0,0,0,0 };
					vector<double> hat_x = sample_x_nodrug(i, j);
					vector<double> hat_y = sample_y_nodrug(i, j);
					vector<int> kx_list = index_x_nodrug(i, j);
					vector<int> ky_list = index_y_nodrug(i, j);

					for (int m = 0; m < fNum_diffpts_interior; m++)
					{
						if (hat_y[m] < fThres_rec)
						{
							hat_v[m] = 1;
						}
						else if (hat_y[m] > fThres_death)
						{
							hat_v[m] = 0;
						}
						else
						{
							if (hat_x[m] < 0)
							{
								// if it diffuses outside q=0, we just project it onto q=0 (never observed actually)
								hat_x[m] = 0;

								array<double, 6> first_col;
								for (int i = 0; i < 6; i++)
								{
									first_col[i] = Vmat_old(i + ky_list[m] - 3, 0);
								}
								// 4th-order ENO cubic interpolation in 1D
								hat_v[m] = myENOgrid1D_y.ENO3_interp_1d(first_col, ky_list[m], hat_y[m]);
							}
							else if (hat_x[m] > 1)
							{
								// if it diffuses outside q=1, we just project it onto q=1 (never observed actually)
								hat_x[m] = 1;

								array<double, 6> last_col;
								for (int i = 0; i < 6; i++)
								{
									last_col[i] = Vmat_old(i + ky_list[m] - 3, fN);
								}
								// 4th-order ENO cubic interpolation in 1D
								hat_v[m] = myENOgrid1D_y.ENO3_interp_1d(last_col, ky_list[m], hat_y[m]);
							}
							else
							{
								// construct the 6x6 matrix of the value function for ENO interpolation in 2D
								ublas::matrix<double> Interp_Mat = myENOgrid2D.Matrix_for_ENO3_Interp(Vmat_old, kx_list[m], ky_list[m]);
								// 4th-order ENO cubic interpolation in 2D
								hat_v[m] = myENOgrid2D.ENO3_interp_2d(Interp_Mat, kx_list[m], ky_list[m], hat_x[m], hat_y[m]);

							}
						}
					}
					double v1 = std::accumulate(hat_v.begin(), hat_v.end(), 0.0);

					v1 = v1 / double(fNum_diffpts_interior); //take the average
					v1 = min(v1, 1.0); //ensuring the largest possible value is v=1


					// For MTD-based Therapy i.e. D = fDmax;
					hat_v = { 0,0,0,0,0,0,0,0 };
					hat_x = sample_x_mtd(i, j);
					hat_y = sample_y_mtd(i, j);
					kx_list = index_x_mtd(i, j);
					ky_list = index_y_mtd(i, j);

					for (int m = 0; m < fNum_diffpts_interior; m++)
					{
						if (hat_y[m] < fThres_rec)
						{
							hat_v[m] = 1;
						}
						else if (hat_y[m] > fThres_death)
						{
							hat_v[m] = 0;
						}
						else
						{
							if (hat_x[m] < 0)
							{
								// if it diffuses outside q=0, we just project it onto q=0 (never observed actually)
								hat_x[m] = 0;

								array<double, 6> first_col;
								for (int i = 0; i < 6; i++)
								{
									first_col[i] = Vmat_old(i + ky_list[m] - 3, 0);
								}
								// 4th-order ENO cubic interpolation in 1D
								hat_v[m] = myENOgrid1D_y.ENO3_interp_1d(first_col, ky_list[m], hat_y[m]);
							}
							else if (hat_x[m] > 1)
							{
								// if it diffuses outside q=1, we just project it onto q=1 (never observed actually)
								hat_x[m] = 1;

								array<double, 6> last_col;
								for (int i = 0; i < 6; i++)
								{
									last_col[i] = Vmat_old(i + ky_list[m] - 3, fN);
								}
								// 4th-order ENO cubic interpolation in 1D
								hat_v[m] = myENOgrid1D_y.ENO3_interp_1d(last_col, ky_list[m], hat_y[m]);
							}
							else
							{
								// construct the 6x6 matrix of value function for ENO interpolation in 2D
								ublas::matrix<double> Interp_Mat = myENOgrid2D.Matrix_for_ENO3_Interp(Vmat_old, kx_list[m], ky_list[m]);
								// 4th-order ENO cubic interpolation in 2D
								hat_v[m] = myENOgrid2D.ENO3_interp_2d(Interp_Mat, kx_list[m], ky_list[m], hat_x[m], hat_y[m]);

							}
						}
					}

					double v2 = std::accumulate(hat_v.begin(), hat_v.end(), 0.0);

					v2 = v2 / double(fNum_diffpts_interior); //take the average
					v2 = min(v2, 1.0); //ensuring the largest possible value is v=1


					Vmat_current(i, j) = max(v1, v2); //take the larger value

					if ((v1 >= v2) || ((std::abs(v1 - v2) / min(v1, v2)) < fTol))
					{
						// if the difference is small, we don't use drugs
						policy_mat(i, j) = 0;
					}
					else
					{
						policy_mat(i, j) = 1;
					}

				}

			}
		}

		Vmat_old = Vmat_current; //swap old slice and current slice at the end

	   // save the data for every X slices
		if (l % fStorage_factor == 0)
		{
			//io::AppendToFile2D_downsample<double>(filename1, Vmat_old, downsample_size);
			io::AppendToFile2D<double>(filename1, Vmat_old);
			io::AppendToFile2D<bool>(filename2, policy_mat);
		}

		cout << "The current slice is l = " << l << endl;
		auto end = std::chrono::steady_clock::now();
		std::chrono::duration<double> elapsed_seconds = end - start;
		std::cout << "Elapsed time for the current slice: " << elapsed_seconds.count() << "s\n";
	}

	return Vmat_old;
}
