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
 * File: ENO3_1d.cpp
 *
 * Author: MingYi Wang
 *
 * Description: This file contains the implementation of functions that compute 
 * Newton 2nd and 3rd divided difference; cubic interpolation in Newton form;
 * and 4th-order ENO cubic interpolation in 1D.
 *
 *============================================================================*/

#include "ENO3_1d.h"

//------------------------------Libraries--------------------------------------
#include <boost/math/interpolators/makima.hpp>
using boost::math::interpolators::makima;


// This function computes the second divided difference
// D2(input) : the second 0-th order divided difference
// D12(input) : the first 1st order divided difference
// aStencilArray(input) : 3x1 vector of sample points
// aValueArray(input) : 3x1 vector of corresponding values of sample points
//
// D3(output) : the third 0-th order divided difference
// D23(output) : the second 1st order divided difference
// D123(output) : the first 2nd order divided difference

tuple<double, double, double>
 ENO1D::Second_Divided_Difference(const double D2, const double D12,
	const array<double, 4>& aStencilArray, const array<double, 4>& aValueArray)
{
	// D123 : = f[x1, x2, x3]
	//        = (f[x2, x3] - f[x1, x2]) / (x3 - x1).
	double D3 = aValueArray[2];
	double D23 = (D3 - D2) / (aStencilArray[2] - aStencilArray[1]);
	double D123 = (D23 - D12) / (aStencilArray[2] - aStencilArray[0]);
	return make_tuple(D3, D23, D123);
}

// This function computes the third divided difference
// D3(input) : the third 0-th order divided difference
// D23(input) : the second 1st order divided difference
// D123(input) : the first 2nd order divided difference
// aStencilArray(input) : 4x1 vector of sample points
// aValueArray(input) : 4x1 vector of corresponding values of sample points
//
// D1234(output) : the first 3rd order divided difference
double
 ENO1D::Third_Divided_Difference(const double D3, const double D23, const double D123,
	const array<double, 4>& aStencilArray, const array<double, 4>& aValueArray)
{
	// D1234 : = f[x1, x2, x3, x4]
	//         = (f[x2, x3, x4] - f[x1, x2, x3]) / (x4 - x1).
	double D4 = aValueArray[3];
	double D34 = (D4 - D3) / (aStencilArray[3] - aStencilArray[2]);
	double D234 = (D34 - D23) / (aStencilArray[3] - aStencilArray[1]);
	double D1234 = (D234 - D123) / (aStencilArray[3] - aStencilArray[0]);
	return D1234;
}

// This function computes cubic interpolation in Newton form
// aStencilArray(input) : 4x1 vector of sample points
// aCoefficientArray(input) : 4x1 coefficient vector of x computed from Newton divided difference.
// xloc(input) : coordinate of the query point
// val(output) : interpolated value
double ENO1D::NewtonInterp(const array<double, 4>& aStencilArray, const array<double, 4>& aCoefficientArray, const double xloc)
{
	const int n = aCoefficientArray.size() - 1; //degree of the interpolating polynomial
	double val = aCoefficientArray[n];
	for (int i = 0; i < n; i++)
	{
		val = aCoefficientArray[n - i - 1] + (xloc - aStencilArray[n - i - 1]) * val;
	}
	return val;
}


// This function constructs a value array with length 6 for ENO cubic interpolation in 1D
// aValueFunctionArray(input): a (N+1)x1 vector of values corresponding to a uniform grid in 1D
// kIndex(input): the first index such that xloc <= x, where x is the vector of sample points
//
// ValueArray(output): a 6x1 vector of values corresponding to the 6-point stencil to be used for ENO interpolation
//
array<double,6> inline ENO1D::Array_for_ENO3_Interp(const ublas::vector<double>& aValueFunctionArray, const int kIndex)
{
	array<double, 6> ValueArray{ 0,0,0,0,0,0 };

	if ((kIndex == 1) || (kIndex == 2))
	{
		for (int i = 0; i < 4; i++)
		{
			ValueArray[i] = aValueFunctionArray[i];
		}
	}
	else if ((kIndex == fN - 1) || (kIndex == fN))
	{
		for (int i = 0; i < 4; i++)
		{
			ValueArray[i] = aValueFunctionArray[i + fN - 3];
		}
	}
	else
	{
		for (int i = 0; i < 6; i++)
		{
			ValueArray[i] = aValueFunctionArray[i + kIndex - 3];
		}
	}
	return ValueArray;
}

// This function computes a more efficient version of ENO cubic interpolation on a uniform grid
// aValueArray(input) : a 6x1 vector of values of the 6-point stencil to be used for ENO3 interpolation
// xloc(input) : coordinate of the query point
// kIndex(input) : the first index such that xloc <= x, where x is the vector of sample points
//
// val(output) : interpolated value
//
double  ENO1D::ENO3_interp_1d(const array<double, 6>& aValueArray, const int kIndex, const double xloc)
{

	const double NN = fN; // cast the length of position vector into double for future computation
	const double kk = kIndex; // cast the index into double for future computation
	double val = 0; // initialize the interpolated value to be 0
	char flag = 'r'; //initialize the flag to be 'r', which stands for the right stencil

	// initialize arrays with length 4 as storage
	array<double, 4> left_stencil;
	array<double, 4> right_stencil;
	// values at left stencil and right stencil respectively
	array<double, 4> yleft;
	array<double, 4> yright;



	if ((kIndex == 1) || (kIndex == 2))
	{
		// if kIndex = 1 or kIndex = 2, we simply choose the first 4-point stencil and use Makima interpolation
		vector<double> stencil{ fx0,fx0 + fDx,fx0 + 2 * fDx,fx0 + 3 * fDx };
		vector<double> temp_value_arrary{ aValueArray[0],aValueArray[1],aValueArray[2],aValueArray[3] };
		makima<vector<double>> interp(std::move(stencil), std::move(temp_value_arrary));
		val = interp(xloc);
	}
	else if ((kIndex == fN - 1) || (kIndex == fN))
	{
		// if kIndex = N-1 or kIndex = N, we simply choose the last 4-point stencil AND use Makima interpolation
		vector<double> stencil{ fx0 + (NN - 3) * fDx,fx0 + (NN - 2) * fDx,fx0 + (NN - 1) * fDx,fx0 + NN * fDx };
		vector<double> temp_value_arrary{ aValueArray[0],aValueArray[1],aValueArray[2],aValueArray[3] };
		makima<vector<double>> interp(std::move(stencil), std::move(temp_value_arrary));
		val = interp(xloc);
	}
	else
	{
		// for given kIndex, set the original two-point stencil to be [x(kIndex - 1), x(kIndex)]
		//  since the grid is uniform, and fx0 = x(0), then x(k) = fx0 + kIndex * fDx
		//  Here as we passed the whole 6-point stencil in, the value function corresponding to
		//  the first two-point stencil is fixed to be [aValueArray(2), aValueArray(3)]
		left_stencil[0] = fx0 + (kk - 1) * fDx; left_stencil[1] = fx0 + kk * fDx;
		right_stencil = left_stencil;

		yleft[0] = aValueArray[2]; yleft[1] = aValueArray[3];
		yright = yleft;

		// add one point either from the left or from the right to form left and
		//	right stencil respectively. Note here we just concatenate them because
		//	 the order does not matter for divided differences
		left_stencil[2] = fx0 + (kk - 2) * fDx;
		right_stencil[2] = fx0 + (kk + 1) * fDx;
		yleft[2] = aValueArray[1];
		yright[2] = aValueArray[4];


		// compute first divided difference
		double D1 = yleft[0]; double D2 = yleft[1];
		double D12 = (D2 - D1) / fDx;

		// compute the second divided differences
		double D3l, D23l, D123l, D3r, D23r, D123r, D3, D23, D123;
		std::tie(D3l, D23l, D123l) = Second_Divided_Difference(D2, D12, left_stencil, yleft);
		std::tie(D3r, D23r, D123r) = Second_Divided_Difference(D2, D12, right_stencil, yright);

		if (abs(D123l) <= abs(D123r))
		{
			// if the seond divided difference of the left stencil is less than the one of the right, then choose the left stencil.
			right_stencil = left_stencil; yright = yleft;
			D3 = D3l; D23 = D23l; D123 = D123l;
			flag = 'l'; // change the flag to be 'l', which stands for the left stencil, if we switch to the left

		}
		else // Otherwise, choose the right stencil
		{
			left_stencil = right_stencil; yleft = yright;
			D3 = D3r; D23 = D23r; D123 = D123r;
		}
		// repeat the process to the third-order
		if (flag == 'l')
		{
			left_stencil[3] = fx0 + (kk - 3) * fDx;
			right_stencil[3] = fx0 + (kk + 1) * fDx;
			yleft[3] = aValueArray[0];
			yright[3] = aValueArray[4];

			// compute the third divivded differences
			double D1234l = Third_Divided_Difference(D3, D23, D123, left_stencil, yleft);
			double D1234r = Third_Divided_Difference(D3, D23, D123, right_stencil, yright);

			if (abs(D1234l) <= abs(D1234r))
			{
				array<double, 4> coeff{ D1,D12,D123,D1234l }; // construct the coefficient array
				val = NewtonInterp(left_stencil, coeff, xloc); // compute the interpolated value
			}
			else
			{
				array<double, 4> coeff{ D1,D12,D123,D1234r }; // construct the coefficient array
				val = NewtonInterp(right_stencil, coeff, xloc); // compute the interpolated value
			}
		}
		else
		{
			left_stencil[3] = fx0 + (kk - 2) * fDx;
			right_stencil[3] = fx0 + (kk + 2) * fDx;
			yleft[3] = aValueArray[1];
			yright[3] = aValueArray[5];

			// compute the third divivded differences
			double D1234l = Third_Divided_Difference(D3, D23, D123, left_stencil, yleft);
			double D1234r = Third_Divided_Difference(D3, D23, D123, right_stencil, yright);

			if (abs(D1234l) <= abs(D1234r))
			{
				array<double, 4> coeff{ D1,D12,D123,D1234l }; // construct the coefficient array
				val = NewtonInterp(left_stencil, coeff, xloc); // compute the interpolated value

			}
			else
			{
				array<double, 4> coeff{ D1,D12,D123,D1234r }; // construct the coefficient array
				val = NewtonInterp(right_stencil, coeff, xloc); // compute the interpolated value
			}
		}
	}
	return val;
}
