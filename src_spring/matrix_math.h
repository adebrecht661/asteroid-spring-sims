/*
 * matrix_math.h
 *
 * Implements various mathematical operations, mostly on matrices.
 * 3x3 symmetric matrices only. 3-vectors only.
 *
 *  Created on: Mar 26, 2020
 *      Author: alex
 */

#ifndef SRC_SPRING_MATRIX_MATH_H_
#define SRC_SPRING_MATRIX_MATH_H_

#include <iostream>

class Matrix; // forward declaration

class Vector {
	double array[3];
public:
	// Constructors
	Vector();
	Vector(double scalar);
	Vector(double array[3]);
	Vector(std::initializer_list<double> list);
	Vector(const Vector& vector);

	// Overloaded operators
	Vector operator+(Vector rhs); // Vector addition
	Vector operator+(double scalar); // Add scalar to each term of vector
	Vector operator-(Vector rhs); // Vector subtraction
	Vector operator-(double scalar); // Subtract scalar from each term of vector
	Vector operator*(double scalar); // Multiply vector by scalar, scalar RHS
	friend Vector operator*(Matrix lhs, Vector rhs); // Matrix-vector multiplication
	Vector operator/(double scalar); // Divide vector by scalar

	// Getters
	double getX() const;
	double getY() const;
	double getZ() const;

	// Misc functions
	double len();
	friend Matrix outer(Vector lhs, Vector rhs); // Outer product (much easier with access to private var)
};

class Matrix {
	double array[3][3];
public:
	// Constructors
	Matrix();
	Matrix(double array[3][3]);
	Matrix(std::initializer_list<std::initializer_list<double>> list);
	Matrix(const Matrix& matrix);

	// Overloaded operators
	Matrix operator+(Matrix rhs); // Matrix addition
	Matrix operator-(Matrix rhs); // Matrix subtraction
	friend Vector operator*(Matrix lhs, Vector rhs); // Matrix-vector multiplication
	friend Matrix operator*(Matrix lhs, Matrix rhs);
	Matrix operator*(double scalar); // Multiply matrix by scalar, scalar RHS
	Matrix operator/(double scalar); // Divide matrix by scalar

	// Getters
	double getXX() const;
	double getXY() const;
	double getXZ() const;
	double getYX() const;
	double getYY() const;
	double getYZ() const;
	double getZX() const;
	double getZY() const;
	double getZZ() const;

	// Misc functions
	bool isSym(); // Check if matrix is symmetric
};

/*******************************/
/* Vector non-member functions */
/*******************************/

double dot(Vector lhs, Vector rhs); // Dot product
Vector cross(Vector lhs, Vector rhs); // Cross product
Matrix outer(Vector lhs, Vector rhs); // Outer product

Vector operator*(double scalar, Vector rhs); // Multiply vector by scalar, scalar LHS
void operator+=(Vector &lhs, double scalar); // Casts scalar to vector
void operator+=(Vector &lhs, Vector rhs); // Add vector
void operator-=(Vector &lhs, double scalar); // Casts scalar to vector
void operator-=(Vector &lhs, Vector rhs); // Subtract vector
void operator*=(Vector &lhs, double scalar); // Multiply vector by scalar
void operator/=(Vector &lhs, double scalar); // Divide vector by scalar
std::ostream& operator<<(std::ostream& os, const Vector& vec); // Stream output
std::istream& operator>>(std::istream& is, Vector& vec); // Stream input

/*******************************/
/* Matrix non-member functions */
/*******************************/

double det(Matrix mat);
Matrix inverse(Matrix mat); // Take inverse of matrix
double trace(Matrix mat); // Get trace of matrix
void eigenvalues(Matrix mat, double eigvals[3]); // Get eigenvalues of matrix
Vector eigenvector(Matrix mat, double eigenvalue); // Get eigenvector corresponding to given eigenvalue

Matrix operator*(double scalar, Matrix rhs); // Multiply matrix by scalar, scalar LHS
//Matrix operator*(Matrix lhs, Matrix rhs); // Matrix multiplication
void operator*=(Matrix &lhs, Matrix rhs); // Matrix multiplication
void operator*=(Matrix &lhs, double scalar); // Multiply matrix by scalar
void operator/=(Matrix &lhs, double scalar); // Divide matrix by scalar
void operator+=(Matrix &lhs, Matrix rhs); // Add matrices
void operator-=(Matrix &lhs, Matrix rhs); // Subtract matrices
std::ostream& operator<<(std::ostream& os, const Matrix& vec); // Stream output
std::istream& operator>>(std::istream& is, Matrix& vec); // Stream input

/*****************************/
/* Matrix & Vector Constants */
/*****************************/

// Declare a single identity matrix
extern double identity[3][3];
extern double zero_vec[3];
extern double zero_mat[3][3];
extern const Matrix I;
extern const Vector zero_vector;
extern const Matrix zero_matrix;

#endif /* SRC_SPRING_MATRIX_MATH_H_ */
