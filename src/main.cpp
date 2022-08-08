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
#include "SemiLagrangian_Cancer.h"
#include "WriteToFile.h"

//----------------------------Libraries------------------------------------------
#include <chrono>
#include <boost/numeric/ublas/assignment.hpp>

using namespace std;

int main()
{
	const int gFactor = 16;
	const double gBudget = 6;
	const double gTreatment_const = 0.05;
	const double gDiff_const = 0.15;

	CancerSL Example(gFactor, gBudget, gTreatment_const, gDiff_const);
	ublas::matrix<double> v = Example.MainSolver_by_SL();

	return 0;
}


