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
  * File: main.cpp
  *
  * Author: MingYi Wang
  *
  * Description: This file initializes all the global variables and executes the
  * the corresponding example from the command line.
  *
  *============================================================================*/

  //-----------------------Project specific header files---------------------------
#include "ENO3_1d.h"
#include "ENO3_2d.h"
#include "Solver_SL.h"
#include "WriteToFile.h"

//----------------------------Libraries------------------------------------------
#include <chrono>
#include <boost/numeric/ublas/assignment.hpp>

using namespace std;

int main()
{
	const int gFactor = 16; // Global factor to determine the spatial grid resolution. E.g., gFactor = 4 means a 401 x 401 grid
	const double gTreatment_const = 0.05; // Treatment constant "delta"
	const double gDiff_const = 0.15; // Diffusion constant (here we assumed all volatilities are the same. But feel free to change the setting)
	const double gThres_recovery = 0.01; // Stabilization/remission barrier
	const double gThres_failure = 1.0 - gThres_recovery; // Failure barrier
	const string gCase = "Example2"; // or "Example2" // Choose which Example to run. Current setting: (1) EGT-based model and (2) SR model

	if (gCase == "Example1") {
		const double gSlice_factor = 125.0 / 2; // Global factor to determine how many s-slices to implement
		const int gNum_diffpts_interior = 8; // Number of diffusive sample points for each gridpoint in the interior
		const int gNum_diffpts_boundary = 4; // Number of diffusive sample points for each gridpoint on the boundary
		const double gMaxDose = 3; // MTD rate
		const double gBudget = 6; // Total budget
		//=============================================Start of Our PDE Solver=====================================================
		auto start = chrono::high_resolution_clock::now();
		cout << "Running our PDE Solver for Example 1" << endl;
		CancerSL Example(gFactor, gSlice_factor, gBudget, gTreatment_const, gDiff_const,
			gNum_diffpts_interior, gNum_diffpts_boundary, gThres_recovery, gThres_failure, gMaxDose, gCase);
		ublas::matrix<double> v = Example.MainSolver_by_SL();
		auto end = chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed_seconds = end - start;
		cout << "Total runtime: " << elapsed_seconds.count() << "s\n";
		cout << "Successfully completed!" << endl;
		cout << "========================================" << endl;
		//=============================================END of Our PDE Solver=====================================================
	}
	else if (gCase == "Example2") {
		const double gSlice_factor = 250; // Global factor to determine how many s-slices to implement
		const int gNum_diffpts_interior = 2; // Number of diffusive sample points for each gridpoint in the interior
		const int gNum_diffpts_boundary = 2; // Number of diffusive sample points for each gridpoint on the boundary
		const double gTime_rescale = 15; // Rescaling factor for time
		const double gMaxDose = 3; // MTD rate
		const double gAlpha = 0.06*gTime_rescale; // Rescaled durg efficiency
		const double gBeta = 6.25*gTime_rescale; // Rescaled action of senstivie on resistant (of order 1e-7)
		const double gCarrying_capacity = 0.48; // Petri dish Carrying Capacity (of order 1e7)
		const double gBudget = 6; // Total budget
		//=============================================Start of Our PDE Solver=====================================================
		auto start = chrono::high_resolution_clock::now();
		cout << "Running our PDE Solver for Example 2" << endl;
		CancerSL Example(gFactor, gSlice_factor, gBudget, gTreatment_const, gDiff_const,
			gNum_diffpts_interior, gNum_diffpts_boundary, gThres_recovery, gThres_failure, gMaxDose,
			gAlpha, gBeta, gCarrying_capacity,gTime_rescale, gCase);
		ublas::matrix<double> v = Example.MainSolver_by_SL();
		auto end = chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed_seconds = end - start;
		cout << "Total runtime: " << elapsed_seconds.count() << "s\n";
		cout << "Successfully completed!" << endl;
		cout << "========================================" << endl;
		//=============================================END of Our PDE Solver=====================================================
	}


	return 0;
}


