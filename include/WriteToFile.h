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
  * File: WriteToFile.h
  *
  * Author: MingYi Wang (based on the code by Marc Aur¨¨le Gilles)
  *
  * Description: This file contains helper functions for writing multi-dimensional
  * Boost arrays and vectors to file
  *
  *============================================================================*/


#ifndef WRITE_TO_FILE_H
#define WRITE_TO_FILE_H

  /** ----- Libraries ----------------------------------------------------------*/
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
namespace ublas = boost::numeric::ublas;

namespace io {


	//This function writes the 2D Boost::ublas matrix "aMatrix" to a file with name aFilename

	template <class T>
	void writeToFile2D(std::string aFilename, ublas::matrix<T> aMatrix) {
		const int aDim0 = aMatrix.size1();
		const int aDim1 = aMatrix.size2();
		//double element;
		ofstream dataFile;
		dataFile.open(aFilename.c_str(), std::ios::app);
		if (!dataFile) {
			// if can't open the file
			std::cout << "An error occured trying to open the file" << std::endl;
			return;
		}
		else {
			// Success!
			dataFile.close();
		}

		dataFile.open(aFilename.c_str(), std::ios::binary);
		for (int i = 0; i < aDim0; ++i) {
			for (int j = 0; j < aDim1; ++j)
			{
				dataFile.write(reinterpret_cast<char*> (&aMatrix(i, j)), sizeof(T));
			}
		}
		dataFile.close();
	}

	template <class T>
	void AppendToFile2D(std::string aFilename, ublas::matrix<T> aMatrix) {
		const int aDim0 = aMatrix.size1();
		const int aDim1 = aMatrix.size2();
		//double element;
		ofstream dataFile;
		dataFile.open(aFilename.c_str(), std::ios::binary | std::ios::app);
		for (int i = 0; i < aDim0; ++i) {
			for (int j = 0; j < aDim1; ++j)
			{
				dataFile.write(reinterpret_cast<char*> (&aMatrix(i, j)), sizeof(T));
			}
		}
		dataFile.close();
	}


	// This function writes the 1D vector "aVec" to a file with name aFilename

	template <class T>
	void writeVectorToFile(std::string aFilename, std::vector<T> aVec) {
		ofstream dataFile;
		dataFile.open(aFilename.c_str(), std::ios::app);
		if (!dataFile) {
			// if can't open the file
			std::cout << "An error occured trying to open the file" << std::endl;
			return;
		}
		else {
			// Success!
			dataFile.close();
		}
		dataFile.open(aFilename.c_str(), std::ios::binary);
		for (int i = 0; i < aVec.size(); ++i) {
			dataFile.write(reinterpret_cast<char*> (&aVec[i]), sizeof(T));
		}
		dataFile.close();
	}



	//This function writes the 2D Boost::ublas matrix "aMatrix" to a file with name aFilename

	template <class T>
	void writeToFile2D_downsample(std::string aFilename, ublas::matrix<T> aMatrix) {
		const int aDim0 = aMatrix.size1();
		const int aDim1 = aMatrix.size2();
		//double element;
		ofstream dataFile;
		dataFile.open(aFilename.c_str(), std::ios::app);
		if (!dataFile) {
			// if can't open the file
			std::cout << "An error occured trying to open the file" << std::endl;
			return;
		}
		else {
			// Success!
			dataFile.close();
		}

		dataFile.open(aFilename.c_str(), std::ios::binary);
		for (int i = 0; i < aDim0; i += 2) {
			for (int j = 0; j < aDim1; j += 2)
			{
				dataFile.write(reinterpret_cast<char*> (&aMatrix(i, j)), sizeof(T));
			}
		}
		dataFile.close();
	}

	template <class T>
	void AppendToFile2D_downsample(std::string aFilename, ublas::matrix<T> aMatrix) {
		const int aDim0 = aMatrix.size1();
		const int aDim1 = aMatrix.size2();
		//double element;
		ofstream dataFile;
		dataFile.open(aFilename.c_str(), std::ios::binary | std::ios::app);
		for (int i = 0; i < aDim0; i += 2) {
			for (int j = 0; j < aDim1; j += 2)
			{
				dataFile.write(reinterpret_cast<char*> (&aMatrix(i, j)), sizeof(T));
			}
		}
		dataFile.close();
	}



	//This function writes the 2D Boost::ublas matrix "aMatrix" to a file with name aFilename

	template <class T>
	void writeToFile2D_recast(std::string aFilename, ublas::matrix<T> aMatrix) {
		const int aDim0 = aMatrix.size1();
		const int aDim1 = aMatrix.size2();
		//double element;
		ofstream dataFile;
		dataFile.open(aFilename.c_str(), std::ios::app);
		if (!dataFile) {
			// if can't open the file
			std::cout << "An error occured trying to open the file" << std::endl;
			return;
		}
		else {
			// Success!
			dataFile.close();
		}

		dataFile.open(aFilename.c_str(), std::ios::binary);
		for (int i = 0; i < aDim0; i++) {
			for (int j = 0; j < aDim1; j++)
			{
				short p = (short)round(aMatrix(i,j) * 1e4);
				//dataFile.write(reinterpret_cast<char*> (&aMatrix(i, j)), sizeof(T));
				dataFile.write(reinterpret_cast<char*> (&p), sizeof(short));
			}
		}
		dataFile.close();
	}

	template <class T>
	void AppendToFile2D_recast(std::string aFilename, ublas::matrix<T> aMatrix) {
		const int aDim0 = aMatrix.size1();
		const int aDim1 = aMatrix.size2();
		//double element;
		ofstream dataFile;
		dataFile.open(aFilename.c_str(), std::ios::binary | std::ios::app);
		for (int i = 0; i < aDim0; i++) {
			for (int j = 0; j < aDim1; j++)
			{
				dataFile.write(reinterpret_cast<char*> (&aMatrix(i, j)), sizeof(T));
			}
		}
		dataFile.close();
	}

} /* namespace io */

#endif // !WRITE_TO_FILE_H
