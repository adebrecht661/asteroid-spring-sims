/*
 * test.cpp
 *
 *  Created on: Mar 31, 2020
 *      Author: alex
 */

#include <iostream>
#include "matrix_math.h"

int main() {
	Vector left(zero_vec);
	Vector right(zero_vec);

	std::cout << left << std::endl;
	std::cout << right << std::endl;

	left = 1;
	std::cout << left << std::endl;
	right += left;
	std::cout << right << std::endl;
	left -= 2;
	std::cout << left << std::endl;

	Vector testvec{0, 0, 0};
	std::cout << testvec << std::endl;

	Matrix testmat{{1,2,3},{4,5,6},{7,8,9}};
	std::cout << testmat << std::endl;
}

