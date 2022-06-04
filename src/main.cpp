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
	//ENO1D my1Dgrid(0.005, 0, 200);
	//ublas::vector<double> v(201);
	//for (int i = 0; i < 201; i++)
	//{
	//	v(i) = sqrt(double(i));
	//}
	//const int k = 55;
	//double val = 0;
	////vector<double> valarray = my1Dgrid.Array_for_ENO3_Interp(v, k);
	//array<double, 6> valarray = my1Dgrid.Array_for_ENO3_Interp(v, k);
	//auto start = std::chrono::steady_clock::now();
	//for (int i = 0; i < 1000001; i++)
	//{
	//	val = my1Dgrid.ENO3_interp_1d(valarray, k, 0.272);
	//}
	//auto end = std::chrono::steady_clock::now();
	//std::chrono::duration<double> elapsed_seconds = end - start;
	//std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
	//cout << val << endl;


	//ENO2D my2Dgrid(200, 200, 0.005, 0.005, 0, 0);
	//vector<double> x(201);
	//vector<double> y(201);
	//for (int i = 0; i < 201; i++)
	//{
	//	x[i] = i * 0.005;
	//	y[i] = i * 0.005;
	//}
	//ublas::matrix<double> mat(201, 201);
	//for (int i = 0; i < 201; i++)
	//{
	//	for (int j = 0; j < 201; j++)
	//	{
	//		mat(i, j) = exp(-5 * (pow(x[i], 2) + pow(y[j], 2)));
	//	}
	//}
	//int kx = 23;
	//int ky = 3;

	//auto start = std::chrono::steady_clock::now();
	//for (int m = 0; m < 1001; m++)
	//{
		/*ublas::matrix<double> cf = my2Dgrid.Matrix_for_ENO3_Interp(mat, kx, ky);

		double vv = my2Dgrid.ENO3_interp_2d(cf, kx, ky, 0.111, 0.011);
		cout << vv << endl;*/
	//}
	//auto end = std::chrono::steady_clock::now();
	//std::chrono::duration<double> elapsed_seconds = end - start;
	//std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

	const int factor = 16;
	const double budget = 6;
	const double treatment_const = 0.05;
	const double diff_const = 0.5;

	CancerSL Example(factor, budget, treatment_const, diff_const);
	ublas::matrix<double> v = Example.MainSolver_by_SL();
	//cout << v << endl;
	//cout << int(125.0/2*factor*budget) << endl;

	return 0;

}


