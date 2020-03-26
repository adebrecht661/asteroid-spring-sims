/*
 * math.h
 *
 * Implements various mathematical operations, mostly on matrices.
 * 3 dimensions only.
 *
 *  Created on: Mar 26, 2020
 *      Author: alex
 */

#ifndef SRC_SPRING_MATH_H_
#define SRC_SPRING_MATH_H_

class Matrix {
	double array[3][3];
public:
	Matrix(double array[3][3]);
	Matrix operator+(Matrix rhs); // Matrix addition
	Matrix operator-(Matrix rhs); // Matrix subtraction
	Matrix operator*(int scalar); // Multiply matrix by scalar, scalar RHS
	Matrix operator/(int scalar); // Divide matrix by scalar
	friend Matrix operator*(int scalar, Matrix rhs); // Multiply matrix by scalar, scalar LHS
};

class Vector {
	double array[3];
public:
	Vector(double array[3]);
	Vector operator+(Vector rhs); // Vector addition
	Vector operator-(Vector rhs); // Vector subtraction
	Vector operator*(int scalar); // Multiply vector by scalar, scalar RHS
	Vector operator/(int scalar); // Divide vector by scalar
	friend Vector operator*(int scalar, Vector rhs); // Multiply vector by scalar, scalar LHS
};

#endif /* SRC_SPRING_MATH_H_ */
