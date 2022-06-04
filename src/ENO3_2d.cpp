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
 * File: ENO3_2d.cpp
 *
 * Author: MingYi Wang
 *
 * Description: This file contains the implementation of functions that construct
 * the stencil for 3rd-order ENO interpolation in 2D and compute the actual 3rd-order
 * ENO interpolation in 2D
 *
 *============================================================================*/

#include "ENO3_2d.h"

//-----------------------------Libraries--------------------------------------
#include <boost/math/interpolators/makima.hpp>
using boost::math::interpolators::makima;

// This function constructs the 6x6 matrix needed for ENO interpolation in 2D
// aValueFunctionMatrix(input): a (N+1)x(N+1) matrix of the value function in 2D on a uniform grid
// kIndex_x(input): the first index such that xloc <= x, where x is the vector of sample points in the horizontal direciton
// kIndex_y(input): the first index such that yloc <= x, where y is the vector of sample points in the vertical direction
//
// aValueMatrix(output): a 6x6 matrix of values corresponding to the 6x6 stencil to be used for ENO3 interpolation in 2D.
ublas::matrix<double> ENO2D::Matrix_for_ENO3_Interp(const ublas::matrix<double>& aValueFunctionMatrix, const int kIndex_x, const int kIndex_y)
{
	// if kIndex_x = 1 or kIndex_x, we simply choose the first 4-point stencil in the horinzontal direction
	if ((kIndex_x == 1) || (kIndex_x == 2))
	{
		ublas::matrix<double> aValueMatrix(6, 4);
		for (int i = 0; i < 6; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				aValueMatrix(i, j) = aValueFunctionMatrix(i + kIndex_y - 3, j);
			}
		}
		return aValueMatrix;
	}
	else if ((kIndex_x == fN - 1) || (kIndex_x == fN))
	{
		// if kIndex_x = fN-1 or kIndex_x = fN, we simply choose the last 4-point stencil in the horinzontal direction
		ublas::matrix<double> aValueMatrix(6, 4);
		for (int i = 0; i < 6; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				aValueMatrix(i, j) = aValueFunctionMatrix(i + kIndex_y - 3, j + fN - 3);
			}
		}
		return aValueMatrix;
	}
	else
	{
		ublas::matrix<double> aValueMatrix(6, 6);
		for (int i = 0; i < 6; i++)
		{
			for (int j = 0; j < 6; j++)
			{
				aValueMatrix(i, j) = aValueFunctionMatrix(i + kIndex_y - 3, j + kIndex_x - 3);
			}
		}
		return aValueMatrix;
	}


}



// This function computes a more efficient verstion of ENO3 interpolation on a uniform grid for interior points in 2D
// aValueMatrix(input) : 6x6 matrix of values of sample points.
// xloc(input) : x coordinate of the query point
// yloc(input) : y coordinate of the query point
// kIndex_x(input) : the first index such that xloc <= x, where x is the vector of
//					  sample points in horizontal direction
// kIndex_y(input) : the first index such that yloc <= y, where y is the vector of
//			          sample points in vertical direction
//
// val(output): the interpolated value
double ENO2D::ENO3_interp_2d(const ublas::matrix<double>& aValueMatrix, const int kIndex_x, const int kIndex_y, const double xloc, const double yloc)
{
	const double NN = fN; // cast the length of position vector in x-direction into double for future computation
	const double kxx = kIndex_x; // cast the index in x-direction into double for future computation (column # of the 2D grid)
	const double kyy = kIndex_y; // cast the index in y-direction into double for future computation (row # of the 2D grid)
	double val = 0; // initialize the interpolated value to be 0
	char flag = 'r'; //initialize the flag to be 'r', which stands for the right stencil

	//initialize arrays with length 4 as storage
	array<double, 4> left_stencil;
	array<double, 4> right_stencil;
	// values at left stencil and right stencil respectively
	array<double, 4> yleft;
	array<double, 4> yright;
	// create an instance of ENO 1D class object in y-direction
	ENO1D grid_y(fDy, fy0, fNy);

	if ((kIndex_x == 1) || (kIndex_x == 2))
	{
		// if kIndex_x = 1 or kIndex_x = 2, we simply choose the first 4-point stencil in the horinzontal direction
		array<double, 6> col1, col2, col3, col4;

		for (int i = 0; i < 6; i++)
		{
			col1[i] = aValueMatrix(i, 0);
			col2[i] = aValueMatrix(i, 1);
			col3[i] = aValueMatrix(i, 2);
			col4[i] = aValueMatrix(i, 3);
		}

		double y1 = grid_y.ENO3_interp_1d(col1, kIndex_y, yloc);
		double y2 = grid_y.ENO3_interp_1d(col2, kIndex_y, yloc);
		double y3 = grid_y.ENO3_interp_1d(col3, kIndex_y, yloc);
		double y4 = grid_y.ENO3_interp_1d(col4, kIndex_y, yloc);

		vector<double> stencil{ fx0,fx0 + fDx, fx0 + 2 * fDx, fx0 + 3 * fDx };
		vector<double> y{ y1,y2,y3,y4 };

		makima<vector<double>> interp(std::move(stencil), std::move(y));
		val = interp(xloc);
	}
	else if ((kIndex_x == fN - 1) || (kIndex_x == fN))
	{
		// if kIndex_x = fN-1 or kIndex_x = fN, we simply choose the last 4-point stencil in the horinzontal direction
		array<double, 6> col1, col2, col3, col4;

		for (int i = 0; i < 6; i++)
		{
			col1[i] = aValueMatrix(i, 0);
			col2[i] = aValueMatrix(i, 1);
			col3[i] = aValueMatrix(i, 2);
			col4[i] = aValueMatrix(i, 3);
		}

		double y1 = grid_y.ENO3_interp_1d(col1, kIndex_y, yloc);
		double y2 = grid_y.ENO3_interp_1d(col2, kIndex_y, yloc);
		double y3 = grid_y.ENO3_interp_1d(col3, kIndex_y, yloc);
		double y4 = grid_y.ENO3_interp_1d(col4, kIndex_y, yloc);

		vector<double> stencil{ fx0 + (NN - 3) * fDx,fx0 + (NN - 2) * fDx, fx0 + (NN - 1) * fDx, fx0 + NN * fDx };
		vector<double> y{ y1,y2,y3,y4 };

		makima<vector<double>> interp(std::move(stencil), std::move(y));
		val = interp(xloc);

	}
	else
	{
		// for given kIndex_x, set the original two-point stencil in the horizontal direction to be [x(kIndex_x - 1), x(kIndex_x)].
		// Since the grid is uniform, and fx0 = x(0), then x(kx) = fx0 + kIndex_x * fDx
		// we will use 1D ENO to find their corresponding values

		array<double, 6> col1, col2;
		for (int i = 0; i < 6; i++)
		{
			col1[i] = aValueMatrix(i, 2);
			col2[i] = aValueMatrix(i, 3);
		}

		double y1 = grid_y.ENO3_interp_1d(col1, kIndex_y, yloc);
		double y2 = grid_y.ENO3_interp_1d(col2, kIndex_y, yloc);

		left_stencil[0] = fx0 + (kxx - 1) * fDx; left_stencil[1] = fx0 + kxx * fDx;
		right_stencil = left_stencil;

		yleft[0] = y1; yleft[1] = y2;
		yright = yleft;

		// Add one point either from the left or from the right to form left and right stencil respectively.
		// Note here we just concatenate them because the order does not matter for divided differences
		left_stencil[2] = fx0 + (kxx - 2) * fDx;
		right_stencil[2] = fx0 + (kxx + 1) * fDx;

		// copy the column at (kIndex_x - 2) and (kIndex_x + 1)
		array<double, 6> left_col, right_col;
		for (int i = 0; i < 6; i++)
		{
			left_col[i] = aValueMatrix(i, 1);
			right_col[i] = aValueMatrix(i, 4);
		}
		// compute the value for the point added from the left by ENO 1D interpolation
		double yl = grid_y.ENO3_interp_1d(left_col, kIndex_y, yloc);
		// compute the value for the point added from the right by ENO 1D interpolation
		double yr = grid_y.ENO3_interp_1d(right_col, kIndex_y, yloc);

		yleft[2] = yl;
		yright[2] = yr;

		// compute first divided difference in horizontal direction
		double D1 = yleft[0]; double D2 = yleft[1];
		double D12 = (D2 - D1) / fDx;
		// compute second divided differences
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
		else // otherwise, choose the right stencil
		{
			left_stencil = right_stencil; yleft = yright;
			D3 = D3r; D23 = D23r; D123 = D123r;
		}
		// repeat the process to the third-order in the horizontal direction
		if (flag == 'l')
		{
			left_stencil[3] = fx0 + (kxx - 3) * fDx;
			right_stencil[3] = fx0 + (kxx + 1) * fDx;

			// copy the column at (kIndex_x - 3)
			for (int i = 0; i < 6; i++)
			{
				left_col[i] = aValueMatrix(i, 0);
			}

			// compute the value for the point added from the left at x(kIndex_x-3) by ENO 1D interpolation
			yl = grid_y.ENO3_interp_1d(left_col, kIndex_y, yloc);
			// Note in this case, the value for the pointed added from the right is the last entry in the previous yr
			yr = yright[2];

			yleft[3] = yl;
			yright[3] = yr;

			// compute the third divivded differences
			double D1234l = Third_Divided_Difference(D3, D23, D123, left_stencil, yleft);
			double D1234r = Third_Divided_Difference(D3, D23, D123, right_stencil, yright);

			if (abs(D1234l) <= abs(D1234r))
			{
				array<double, 4> coeff{ D1,D12,D123,D1234l };
				val = NewtonInterp(left_stencil, coeff, xloc);
			}
			else
			{
				array<double, 4> coeff{ D1,D12,D123,D1234r };
				val = NewtonInterp(right_stencil, coeff, xloc);
			}
		}
		else
		{
			left_stencil[3] = fx0 + (kxx - 2) * fDx;
			right_stencil[3] = fx0 + (kxx + 2) * fDx;

			// copy the column at (kIndex_x + 2)
			for (int i = 0; i < 6; i++)
			{
				right_col[i] = aValueMatrix(i, 5);
			}

			// compute the value for the point added from the right at x(kIndex_x + 2) by ENO 1D interpolation
			yr = grid_y.ENO3_interp_1d(right_col, kIndex_y, yloc);
			// Note in this case, the value for the pointed added from the left is the last entry in the previous yl
			yl = yleft[2];

			yleft[3] = yl;
			yright[3] = yr;

			// compute the third divivded differences
			double D1234l = Third_Divided_Difference(D3, D23, D123, left_stencil, yleft);
			double D1234r = Third_Divided_Difference(D3, D23, D123, right_stencil, yright);

			if (abs(D1234l) <= abs(D1234r))
			{
				array<double, 4> coeff{ D1,D12,D123,D1234l };
				val = NewtonInterp(left_stencil, coeff, xloc);
			}
			else
			{
				array<double, 4> coeff{ D1,D12,D123,D1234r };
				val = NewtonInterp(right_stencil, coeff, xloc);
			}
		}
	}
	return val;
}
