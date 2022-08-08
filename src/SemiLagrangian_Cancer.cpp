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
  * File: SemiLagrangian_Cancer.cpp
  *
  * Author: MingYi Wang
  *
  * Description: This file contains the implementation of functions that compute the
  * drift (deterministic) portion of the cancer dynamics; the coefficients and
  * feet of the 1st-order weak approximation of the cancer dynamics; and the main
  * solver for the value function using semi-Lagrangian scheme.
  *
  *============================================================================*/
#include "SemiLagrangian_Cancer.h"

  //----------------------Project specific header files---------------------------
#include "WriteToFile.h"

//------------------------------Libraries------------------------------------
#include<chrono>
#include<omp.h>

// Define the constructor
CancerSL::CancerSL(int a_Factor, double a_Budget, double a_Treatment_const, double a_Diff_const)
{
	// initialization of the 2D grid
	fMulti_factor = a_Factor;
	fBudget = a_Budget;
	fN = 100 * fMulti_factor;

	fM = int(125.0 / 2 * fMulti_factor * fBudget); //for octed ds


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

	// Diffusion constants (volatilities) and Treatment constants
	//we have only considered cases that all three diffusion constants being the same
	fDiff_const_1 = a_Diff_const;
	fDiff_const_2 = fDiff_const_1;
	fDiff_const_3 = fDiff_const_1;

	fSigma = a_Treatment_const;
	fDmax = 3; //(dmax = 3)

	index_rec = int(0.01 / fDx); // index of the recovery barrier

	index_death = fN - index_rec; // index of the death barrier

	tau_mtd = fDs / (fDmax + fSigma); // time-step for MTD-based therapy
	tau_0 = fDs / fSigma; // time-step for no therapy

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
inline double CancerSL::drift_p(const double q, const double p, const double aDose)
{
	return p * (1 - p) * (ba / (n + 1) - q * (bv - c) - aDose) -
		p * (pow(fDiff_const_1, 2) * p - (pow(fDiff_const_1, 2) * pow(p, 2)
			+ pow(fDiff_const_2, 2) * pow(1 - p, 2) * pow(1 - q, 2) + pow(fDiff_const_3, 2) * pow(1 - p, 2) * pow(q, 2)));
}

inline double CancerSL::drift_q(const double q, const double p)
{
	return q * (1 - q) * (bv / (n + 1) * (1 + p + pow(p, 2) + pow(p, 3) + pow(p, 4)) - c) +
		q * (1 - q) * ((1 - q) * pow(fDiff_const_2, 2) - q * pow(fDiff_const_3, 2));
}

// This function pre-compute all the coefficients and feet of the 1st-order weak approximation of the dynamics
tuple<ublas::matrix<vector<double>>, ublas::matrix<vector<double>>, ublas::matrix<vector<double>>, ublas::matrix<vector<double>>,
	ublas::matrix<vector<int>>, ublas::matrix<vector<int>>, ublas::matrix<vector<int>>, ublas::matrix<vector<int>>>
	CancerSL::precompute_coeff_2D(const ublas::vector<double>& q, const ublas::vector<double>& p)
{
	// initialization of storage
	ublas::matrix<std::vector<double>> sample_x_nodrug(fN + 1, fN + 1); // store all possible locations in x-direction for no therapy
	ublas::matrix<std::vector<double>> sample_x_mtd(fN + 1, fN + 1); // store all possible locations in x-direction for mtd
	ublas::matrix<std::vector<double>> sample_y_nodrug(fN + 1, fN + 1); // store all possible locations in y-direction for no therapy
	ublas::matrix<std::vector<double>> sample_y_mtd(fN + 1, fN + 1); // tore all possible locations in y-direction for mtd
	ublas::matrix<std::vector<int>> index_x_nodrug(fN + 1, fN + 1); // store all possible coefficients at x-direction for no therapy
	ublas::matrix<std::vector<int>> index_x_mtd(fN + 1, fN + 1); // store all possible coefficients at x-direction for mtd
	ublas::matrix<std::vector<int>> index_y_nodrug(fN + 1, fN + 1); // store all possible coefficients at y-direction for no therapy
	ublas::matrix<std::vector<int>> index_y_mtd(fN + 1, fN + 1); // store all possible coefficients at y-direction for mtd

	// square root of time steps
	double root_tau_0 = sqrt(tau_0);
	double root_tau_mtd = sqrt(tau_mtd);
	for (int i = index_rec; i < index_death + 1; i++)
	{
		for (int j = 0; j < fN + 1; j++)
		{
			if (j == 0) // at q = 0
			{
				vector<int> ky_list{ 0,0,0,0 };
				int k = 0;
				vector<double> hat_y{ 0,0,0,0 };

				// for No therapy, i.e. Dose = 0
				double y_no_drug = p(i) + tau_0 * drift_p(0, p(i), 0); // foot of characteristics at (0,p(i)) for no therapy


				hat_y[0] = y_no_drug + root_tau_0 * p(i) * (1 - p(i)) * fDiff_const_1 + root_tau_0 * p(i) * (1 - p(i)) * fDiff_const_2;
				hat_y[1] = y_no_drug - root_tau_0 * p(i) * (1 - p(i)) * fDiff_const_1 + root_tau_0 * p(i) * (1 - p(i)) * fDiff_const_2;
				hat_y[2] = y_no_drug + root_tau_0 * p(i) * (1 - p(i)) * fDiff_const_1 - root_tau_0 * p(i) * (1 - p(i)) * fDiff_const_2;
				hat_y[3] = y_no_drug - root_tau_0 * p(i) * (1 - p(i)) * fDiff_const_1 - root_tau_0 * p(i) * (1 - p(i)) * fDiff_const_2;

				for (int m = 0; m < 4; m++)
				{
					k = find_index(p, hat_y[m]);
					ky_list[m] = k;
				}
				index_y_nodrug(i, j) = ky_list;
				sample_y_nodrug(i, j) = hat_y;

				// for MTD-based therapy, i.e. Dose = fDmax
				double y_mtd = p(i) + tau_mtd * drift_p(0, p(i), fDmax); // foot of characteristics at (0,p(i)) for mtd

				hat_y = { 0,0,0,0 };

				hat_y[0] = y_mtd + root_tau_mtd * p(i) * (1 - p(i)) * fDiff_const_1 + root_tau_mtd * p(i) * (1 - p(i)) * fDiff_const_2;
				hat_y[1] = y_mtd - root_tau_mtd * p(i) * (1 - p(i)) * fDiff_const_1 + root_tau_mtd * p(i) * (1 - p(i)) * fDiff_const_2;
				hat_y[2] = y_mtd + root_tau_mtd * p(i) * (1 - p(i)) * fDiff_const_1 - root_tau_mtd * p(i) * (1 - p(i)) * fDiff_const_2;
				hat_y[3] = y_mtd - root_tau_mtd * p(i) * (1 - p(i)) * fDiff_const_1 - root_tau_mtd * p(i) * (1 - p(i)) * fDiff_const_2;

				ky_list = { 0,0,0,0 };

				for (int m = 0; m < 4; m++)
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

				// for No therapy, i.e. Dose = 0
				double y_no_drug = p(i) + tau_0 * drift_p(1, p(i), 0); // foot of characteristics at (1,p(i)) for no drug

				hat_y[0] = y_no_drug + root_tau_0 * p(i) * (1 - p(i)) * fDiff_const_1 + root_tau_0 * p(i) * (1 - p(i)) * fDiff_const_3;
				hat_y[1] = y_no_drug - root_tau_0 * p(i) * (1 - p(i)) * fDiff_const_1 + root_tau_0 * p(i) * (1 - p(i)) * fDiff_const_3;
				hat_y[2] = y_no_drug + root_tau_0 * p(i) * (1 - p(i)) * fDiff_const_1 - root_tau_0 * p(i) * (1 - p(i)) * fDiff_const_3;
				hat_y[3] = y_no_drug - root_tau_0 * p(i) * (1 - p(i)) * fDiff_const_1 - root_tau_0 * p(i) * (1 - p(i)) * fDiff_const_3;

				for (int m = 0; m < 4; m++)
				{
					k = find_index(p, hat_y[m]);
					ky_list[m] = k;
				}
				index_y_nodrug(i, j) = ky_list;
				sample_y_nodrug(i, j) = hat_y;

				// for MTD-based therapy, i.e. Dose = fDmax
				double y_mtd = p(i) + tau_mtd * drift_p(1, p(i), fDmax); // foot of characteristics at (1,p(i)) for mtd

				hat_y = { 0,0,0,0 };

				hat_y[0] = y_mtd + root_tau_mtd * p(i) * (1 - p(i)) * fDiff_const_1 + root_tau_mtd * p(i) * (1 - p(i)) * fDiff_const_3;
				hat_y[1] = y_mtd - root_tau_mtd * p(i) * (1 - p(i)) * fDiff_const_1 + root_tau_mtd * p(i) * (1 - p(i)) * fDiff_const_3;
				hat_y[2] = y_mtd + root_tau_mtd * p(i) * (1 - p(i)) * fDiff_const_1 - root_tau_mtd * p(i) * (1 - p(i)) * fDiff_const_3;
				hat_y[3] = y_mtd - root_tau_mtd * p(i) * (1 - p(i)) * fDiff_const_1 - root_tau_mtd * p(i) * (1 - p(i)) * fDiff_const_3;

				ky_list = { 0,0,0,0 };

				for (int m = 0; m < 4; m++)
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

				// for No therapy, i.e. D = 0
				double x_no_drug = q(j) + tau_0 * drift_q(q(j), p(i)); // x-coordinate of foot of characteristics at (q(j),p(i))
				double y_no_drug = p(i) + tau_0 * drift_p(q(j), p(i), 0); // y-coordinate of characteristics at (q(j),p(i))

				// 8 samples points of the 1st-order weak approximation of Brownian motions
				hat_x[0] = x_no_drug + root_tau_0 * q(j) * (1 - q(j)) * fDiff_const_3 - root_tau_0 * q(j) * (1 - q(j)) * fDiff_const_2;
				hat_y[0] = y_no_drug + root_tau_0 * p(i) * (1 - p(i)) * fDiff_const_1 + root_tau_0 * p(i) * (1 - p(i)) * (1 - q(j)) * fDiff_const_2
					+ root_tau_0 * p(i) * (1 - p(i)) * q(j) * fDiff_const_3;

				hat_x[1] = x_no_drug + root_tau_0 * q(j) * (1 - q(j)) * fDiff_const_3 + root_tau_0 * q(j) * (1 - q(j)) * fDiff_const_2;
				hat_y[1] = y_no_drug + root_tau_0 * p(i) * (1 - p(i)) * fDiff_const_1 - root_tau_0 * p(i) * (1 - p(i)) * (1 - q(j)) * fDiff_const_2
					+ root_tau_0 * p(i) * (1 - p(i)) * q(j) * fDiff_const_3;

				hat_x[2] = x_no_drug - root_tau_0 * q(j) * (1 - q(j)) * fDiff_const_3 - root_tau_0 * q(j) * (1 - q(j)) * fDiff_const_2;
				hat_y[2] = y_no_drug + root_tau_0 * p(i) * (1 - p(i)) * fDiff_const_1 + root_tau_0 * p(i) * (1 - p(i)) * (1 - q(j)) * fDiff_const_2
					- root_tau_0 * p(i) * (1 - p(i)) * q(j) * fDiff_const_3;

				hat_x[3] = x_no_drug - root_tau_0 * q(j) * (1 - q(j)) * fDiff_const_3 + root_tau_0 * q(j) * (1 - q(j)) * fDiff_const_2;
				hat_y[3] = y_no_drug + root_tau_0 * p(i) * (1 - p(i)) * fDiff_const_1 - root_tau_0 * p(i) * (1 - p(i)) * (1 - q(j)) * fDiff_const_2
					- root_tau_0 * p(i) * (1 - p(i)) * q(j) * fDiff_const_3;

				hat_x[4] = x_no_drug + root_tau_0 * q(j) * (1 - q(j)) * fDiff_const_3 - root_tau_0 * q(j) * (1 - q(j)) * fDiff_const_2;
				hat_y[4] = y_no_drug - root_tau_0 * p(i) * (1 - p(i)) * fDiff_const_1 + root_tau_0 * p(i) * (1 - p(i)) * (1 - q(j)) * fDiff_const_2
					+ root_tau_0 * p(i) * (1 - p(i)) * q(j) * fDiff_const_3;

				hat_x[5] = x_no_drug - root_tau_0 * q(j) * (1 - q(j)) * fDiff_const_3 - root_tau_0 * q(j) * (1 - q(j)) * fDiff_const_2;
				hat_y[5] = y_no_drug - root_tau_0 * p(i) * (1 - p(i)) * fDiff_const_1 + root_tau_0 * p(i) * (1 - p(i)) * (1 - q(j)) * fDiff_const_2
					- root_tau_0 * p(i) * (1 - p(i)) * q(j) * fDiff_const_3;

				hat_x[6] = x_no_drug + root_tau_0 * q(j) * (1 - q(j)) * fDiff_const_3 + root_tau_0 * q(j) * (1 - q(j)) * fDiff_const_2;
				hat_y[6] = y_no_drug - root_tau_0 * p(i) * (1 - p(i)) * fDiff_const_1 - root_tau_0 * p(i) * (1 - p(i)) * (1 - q(j)) * fDiff_const_2
					+ root_tau_0 * p(i) * (1 - p(i)) * q(j) * fDiff_const_3;

				hat_x[7] = x_no_drug - root_tau_0 * q(j) * (1 - q(j)) * fDiff_const_3 + root_tau_0 * q(j) * (1 - q(j)) * fDiff_const_2;
				hat_y[7] = y_no_drug - root_tau_0 * p(i) * (1 - p(i)) * fDiff_const_1 - root_tau_0 * p(i) * (1 - p(i)) * (1 - q(j)) * fDiff_const_2
					- root_tau_0 * p(i) * (1 - p(i)) * q(j) * fDiff_const_3;

				for (int m = 0; m < 8; m++)
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

				double x_mtd = q(j) + tau_mtd * drift_q(q(j), p(i)); // x-coordinate of foot of characteristics at (q(j),p(i))
				double y_mtd = p(i) + tau_mtd * drift_p(q(j), p(i), fDmax); // y-coordinate of characteristics at (q(j),p(i))

				// sample 8 points for Brownian Motions
				hat_x[0] = x_mtd + root_tau_mtd * q(j) * (1 - q(j)) * fDiff_const_3 - root_tau_mtd * q(j) * (1 - q(j)) * fDiff_const_2;
				hat_y[0] = y_mtd + root_tau_mtd * p(i) * (1 - p(i)) * fDiff_const_1 + root_tau_mtd * p(i) * (1 - p(i)) * (1 - q(j)) * fDiff_const_2
					+ root_tau_mtd * p(i) * (1 - p(i)) * q(j) * fDiff_const_3;

				hat_x[1] = x_mtd + root_tau_mtd * q(j) * (1 - q(j)) * fDiff_const_3 + root_tau_mtd * q(j) * (1 - q(j)) * fDiff_const_2;
				hat_y[1] = y_mtd + root_tau_mtd * p(i) * (1 - p(i)) * fDiff_const_1 - root_tau_mtd * p(i) * (1 - p(i)) * (1 - q(j)) * fDiff_const_2
					+ root_tau_mtd * p(i) * (1 - p(i)) * q(j) * fDiff_const_3;

				hat_x[2] = x_mtd - root_tau_mtd * q(j) * (1 - q(j)) * fDiff_const_3 - root_tau_mtd * q(j) * (1 - q(j)) * fDiff_const_2;
				hat_y[2] = y_mtd + root_tau_mtd * p(i) * (1 - p(i)) * fDiff_const_1 + root_tau_mtd * p(i) * (1 - p(i)) * (1 - q(j)) * fDiff_const_2
					- root_tau_mtd * p(i) * (1 - p(i)) * q(j) * fDiff_const_3;

				hat_x[3] = x_mtd - root_tau_mtd * q(j) * (1 - q(j)) * fDiff_const_3 + root_tau_mtd * q(j) * (1 - q(j)) * fDiff_const_2;
				hat_y[3] = y_mtd + root_tau_mtd * p(i) * (1 - p(i)) * fDiff_const_1 - root_tau_mtd * p(i) * (1 - p(i)) * (1 - q(j)) * fDiff_const_2
					- root_tau_mtd * p(i) * (1 - p(i)) * q(j) * fDiff_const_3;

				hat_x[4] = x_mtd + root_tau_mtd * q(j) * (1 - q(j)) * fDiff_const_3 - root_tau_mtd * q(j) * (1 - q(j)) * fDiff_const_2;
				hat_y[4] = y_mtd - root_tau_mtd * p(i) * (1 - p(i)) * fDiff_const_1 + root_tau_mtd * p(i) * (1 - p(i)) * (1 - q(j)) * fDiff_const_2
					+ root_tau_mtd * p(i) * (1 - p(i)) * q(j) * fDiff_const_3;

				hat_x[5] = x_mtd - root_tau_mtd * q(j) * (1 - q(j)) * fDiff_const_3 - root_tau_mtd * q(j) * (1 - q(j)) * fDiff_const_2;
				hat_y[5] = y_mtd - root_tau_mtd * p(i) * (1 - p(i)) * fDiff_const_1 + root_tau_mtd * p(i) * (1 - p(i)) * (1 - q(j)) * fDiff_const_2
					- root_tau_mtd * p(i) * (1 - p(i)) * q(j) * fDiff_const_3;

				hat_x[6] = x_mtd + root_tau_mtd * q(j) * (1 - q(j)) * fDiff_const_3 + root_tau_mtd * q(j) * (1 - q(j)) * fDiff_const_2;
				hat_y[6] = y_mtd - root_tau_mtd * p(i) * (1 - p(i)) * fDiff_const_1 - root_tau_mtd * p(i) * (1 - p(i)) * (1 - q(j)) * fDiff_const_2
					+ root_tau_mtd * p(i) * (1 - p(i)) * q(j) * fDiff_const_3;

				hat_x[7] = x_mtd - root_tau_mtd * q(j) * (1 - q(j)) * fDiff_const_3 + root_tau_mtd * q(j) * (1 - q(j)) * fDiff_const_2;
				hat_y[7] = y_mtd - root_tau_mtd * p(i) * (1 - p(i)) * fDiff_const_1 - root_tau_mtd * p(i) * (1 - p(i)) * (1 - q(j)) * fDiff_const_2
					- root_tau_mtd * p(i) * (1 - p(i)) * q(j) * fDiff_const_3;

				for (int m = 0; m < 8; m++)
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
	return make_tuple(sample_x_nodrug, sample_x_mtd, sample_y_nodrug, sample_y_mtd, index_x_nodrug, index_x_mtd, index_y_nodrug, index_y_mtd);
}


//--------------------------------Main Solver-----------------------------------
ublas::matrix<double> CancerSL::MainSolver_by_SL()
{
	//-------------------------------------Initialization---------------------------------

	double tol = pow(10, -14); //set up a tolerance to choose one of the two values that are close to each other

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

	 //string filename1 = "strict ENO3_v2_2D N=1600 s=6 s1=0.15.dat";
	 //io::writeToFile2D<double>(filename1,Vmat_old);
	 //string filename2 = "strict dval 2D N=1600 s=6 s1=0.15 tol=1e-14.dat";
	 //io::writeToFile2D<bool>(filename2,policy_mat);

	string filename1 = "test_valuefn.dat";
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
	std::tie(sample_x_nodrug, sample_x_mtd, sample_y_nodrug, sample_y_mtd, index_x_nodrug, index_x_mtd, index_y_nodrug, index_y_mtd)
		= precompute_coeff_2D(q, p);

	// end time of computating all coefficients
	auto end = std::chrono::steady_clock::now();

	// print out the computational time used for coefficients
	std::chrono::duration<double> elapsed_seconds = end - start;
	std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
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
		//int max_num_threads = omp_get_max_threads();
		//cout << max_num_threads << endl;

#pragma omp parallel for
		for (int i = index_rec; i < index_death + 1; i++)
		{

			for (int j = 0; j < fN + 1; j++)
			{
				if (j == 0) // on boundary q = 0
				{
					// For No therapy i.e. D = 0
					vector<double> hat_v{ 0,0,0,0 };
					vector<double> hat_y = sample_y_nodrug(i, j);
					vector<int> ky_list = index_y_nodrug(i, j);

					for (int m = 0; m < 4; m++)
					{
						if (hat_y[m] < 0.01)
						{
							hat_v[m] = 1;
						}
						else if (hat_y[m] > 0.99)
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

					v1 = v1 / 4; //take the average
					v1 = min(v1, 1.0); //ensuring the largest possible value is v=1

					// For MTD-based therapy i.e. D = fDmax
					hat_v = { 0,0,0,0 };
					hat_y = sample_y_mtd(i, j);
					ky_list = index_y_mtd(i, j);

					for (int m = 0; m < 4; m++)
					{
						if (hat_y[m] < 0.01)
						{
							hat_v[m] = 1;
						}
						else if (hat_y[m] > 0.99)
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

					v2 = v2 / 4; //take the average
					v2 = min(v2, 1.0); //ensuring the largest possible value is v=1

					Vmat_current(i, j) = max(v1, v2); //take the larger value

					if ((v1 >= v2) || ((std::abs(v1 - v2) / min(v1, v2)) < tol))
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

					for (int m = 0; m < 4; m++)
					{
						if (hat_y[m] < 0.01)
						{
							hat_v[m] = 1;
						}
						else if (hat_y[m] > 0.99)
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

					v1 = v1 / 4; //take the average
					v1 = min(v1, 1.0); //ensuring the largest possible value is v=1

					// For MTD-based therapy i.e. D = fDmax
					hat_v = { 0,0,0,0 };
					hat_y = sample_y_mtd(i, j);
					ky_list = index_y_mtd(i, j);

					for (int m = 0; m < 4; m++)
					{
						if (hat_y[m] < 0.01)
						{
							hat_v[m] = 1;
						}
						else if (hat_y[m] > 0.99)
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

					v2 = v2 / 4; // take the average
					v2 = min(v2, 1.0); //ensuring the largest possible value is v=1

					Vmat_current(i, j) = max(v1, v2); //take the larger value

					if ((v1 >= v2) || ((std::abs(v1 - v2) / min(v1, v2)) < tol))
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

					for (int m = 0; m < 8; m++)
					{
						if (hat_y[m] < 0.01)
						{
							hat_v[m] = 1;
						}
						else if (hat_y[m] > 0.99)
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

					v1 = v1 / 8; //take the average
					v1 = min(v1, 1.0); //ensuring the largest possible value is v=1


					// For MTD-based Therapy i.e. D = fDmax;
					hat_v = { 0,0,0,0,0,0,0,0 };
					hat_x = sample_x_mtd(i, j);
					hat_y = sample_y_mtd(i, j);
					kx_list = index_x_mtd(i, j);
					ky_list = index_y_mtd(i, j);
					for (int m = 0; m < 8; m++)
					{
						if (hat_y[m] < 0.01)
						{
							hat_v[m] = 1;
						}
						else if (hat_y[m] > 0.99)
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

					v2 = v2 / 8; //take the average
					v2 = min(v2, 1.0); //ensuring the largest possible value is v=1


					Vmat_current(i, j) = max(v1, v2); //take the larger value

					if ((v1 >= v2) || ((std::abs(v1 - v2) / min(v1, v2)) < tol))
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

	   // save the data for every 5 slices
		if (l % 5 == 0)
		{
			io::AppendToFile2D<double>(filename1, Vmat_old);
			io::AppendToFile2D<bool>(filename2, policy_mat);
		}

		cout << "The current slice is l = " << l << endl;
		auto end = std::chrono::steady_clock::now();
		std::chrono::duration<double> elapsed_seconds = end - start;
		std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
	}

	return Vmat_old;
}
