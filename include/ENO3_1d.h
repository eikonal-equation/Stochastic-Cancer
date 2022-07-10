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
 * File: ENO3_1d.h
 *
 * Author: MingYi Wang
 *
 * Description: This file contains the declarations of functions that compute 
 * Newton 2nd and 3rd divided difference; cubic interpolation in Newton form;
 * and 4th-order ENO cubic interpolation in 1D.

 *
 * Details of all of these functions are found in ENO3_1d.cpp.
 *
 *============================================================================*/

#pragma once
#ifndef ENO3_1D_H
#define ENO3_1D_H

//---------------------------Libraries-----------------------------------------
#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <tuple>
#include <chrono>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

//---------------------------Definitions---------------------------------------
namespace ublas = boost::numeric::ublas;
using namespace std;

//------------------------------Main class------------------------------------
class ENO1D
{
public:
	ENO1D() = default;

	ENO1D(double a_Delta_x, double a_x_initial, int a_Dimension_1D)
	{
		fDx = a_Delta_x;
		fx0 = a_x_initial;
		fN = a_Dimension_1D;
	}

	// This function computes the second divided difference
	tuple<double, double, double> Second_Divided_Difference(const double D2, const double D12, const array<double,4>& aStencilArray, const array<double, 4>& aValueArray);

	// This function computes the third divided difference
	double Third_Divided_Difference(const double D3, const double D23, const double D123, const array<double, 4>& aStencilArray, const array<double, 4>& aValueArray);

	// This function is the cubic interpolation in Newton form given a 4-point stencil
	double NewtonInterp(const array<double, 4>& aStencilArray, const array<double, 4>& aCoefficientArray, const double xloc);

	// Construct an array with length 6 for ENO3 1D interpolation
	array<double, 6> Array_for_ENO3_Interp(const ublas::vector<double>& aValueFunctionArray, const int kIndex);

	// Main ENO3 interploation function in 1D
	double ENO3_interp_1d(const array<double, 6>& aValueArray, const int kIndex, const double xloc);



protected:
	double fDx; // Delta x or h
	double fx0; // starting position
	int fN; // length of the position array in 1D
};



#endif // !ENO3_1D_H

