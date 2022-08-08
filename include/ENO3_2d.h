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
  * File: ENO3_2d.h
  *
  * Author: MingYi Wang
  *
  * Description: This file contains the declarations of functions that construct
  * the stencil for 4th-order ENO cubic interpolation in 2D and compute the actual
  * 4th-order ENO cubic interpolation in 2D
  *
  * Details of all of these functions are found in ENO3_2d.cpp.
  *
  *============================================================================*/

#pragma once
#ifndef ENO3_2D_H
#define ENO3_2D_H

  //----------------------Project specific header files---------------------------
#include "ENO3_1d.h"

//----------------------Libraries----------------------------------------------
#include <functional>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

//--------------------------------Main Class-----------------------------------
class ENO2D : public ENO1D
{
public:
	// default constructor
	ENO2D() = default;
	ENO2D(int a_Dimension_x, int a_Dimension_y, double a_Delta_x, double a_Delta_y, double a_x_initial, double a_y_initial) : ENO1D(a_Delta_x, a_x_initial, a_Dimension_x)
	{
		fNy = a_Dimension_y;
		fDy = a_Delta_y;
		fy0 = a_y_initial;
	}

	// This function constructs the 6x6 matrix needed for ENO interpolation in 2D
	ublas::matrix<double> Matrix_for_ENO3_Interp(const ublas::matrix<double>& aValueFunctionMatrix, const int kIndex_x, const int kIndex_y);

	// Main ENO cubic interpolation function in 2D
	double ENO3_interp_2d(const ublas::matrix<double>& aValueMatrix, const int kIndex_x, const int kIndex_y, const double xloc, const double yloc);


private:
	double fDy; // Delta y or h in y-direction
	double fy0; // starting point on y-axis
	int fNy; // length of the side on y-direction

};


#endif // ! ENO3_2D_H

