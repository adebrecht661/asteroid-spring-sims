/*
 * Matrix.cpp
 *
 *  Created on: Mar 26, 2020
 *      Author: alex
 */

#include "math.h"

// Declare a single identity matrix
double identity[3][3] = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
const static Matrix I(identity);

Matrix::Matrix(double array[3][3]) {
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			this->array[i][j] = array[i][j];
		}
	}
}

