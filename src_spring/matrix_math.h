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

// Forward declaration of matrix class (used in Vector class)
class Matrix;

class Vector {
	double array[3];
public:
	/****************/
	/* Constructors */
	/****************/

	// Null constructor
	Vector();
	// Construct from scalar
	Vector(double scalar);
	// Construct from array
	Vector(double array[3]);
	// Construct from initializer list
	Vector(std::initializer_list<double> list);
	// Copy constructor
	Vector(const Vector &vector);

	/************************/
	/* Overloaded operators */
	/************************/

	// Vector addition
	Vector operator+(Vector rhs);
	// Add scalar to each term of vector
	Vector operator+(double scalar);
	// Vector subtraction
	Vector operator-(Vector rhs);
	// Subtract scalar from each term of vector
	Vector operator-(double scalar);
	// Multiply vector by scalar, scalar RHS
	Vector operator*(double scalar);
	// Matrix-vector multiplication
	friend Vector operator*(Matrix lhs, Vector rhs);
	// Divide vector by scalar
	Vector operator/(double scalar);

	/***********/
	/* Getters */
	/***********/

	double getX() const;
	double getY() const;
	double getZ() const;

	/******************/
	/* Misc functions */
	/******************/

	// Length of vector
	double len();
	// Outer product (much easier with access to private var)
	friend Matrix outer(Vector lhs, Vector rhs);
};

class Matrix {
	double array[3][3];
public:
	/****************/
	/* Constructors */
	/****************/

	// Null constructor
	Matrix();
	// Construct from array
	Matrix(double array[3][3]);
	// Construct from initializer list of initializer lists
	Matrix(std::initializer_list<std::initializer_list<double>> list);
	// Construct from initializer list of (row) vectors
	Matrix(std::initializer_list<Vector> list);
	// Copy constructor
	Matrix(const Matrix &matrix);

	/************************/
	/* Overloaded operators */
	/************************/

	// Matrix addition
	Matrix operator+(Matrix rhs);
	// Matrix subtraction
	Matrix operator-(Matrix rhs);
	// Matrix-vector multiplication
	friend Vector operator*(Matrix lhs, Vector rhs);
	// Matrix-matrix multiplication
	friend Matrix operator*(Matrix lhs, Matrix rhs);
	// Multiply matrix by scalar, scalar RHS
	Matrix operator*(double scalar);
	// Divide matrix by scalar
	Matrix operator/(double scalar);

	/***********/
	/* Getters */
	/***********/

	double getXX() const;
	double getXY() const;
	double getXZ() const;
	double getYX() const;
	double getYY() const;
	double getYZ() const;
	double getZX() const;
	double getZY() const;
	double getZZ() const;

	/******************/
	/* Misc functions */
	/******************/

	// Check if matrix is symmetric
	bool isSym();
};

/*******************************/
/* Vector non-member functions */
/*******************************/

// Dot product
double dot(Vector lhs, Vector rhs);
// Cross product
Vector cross(Vector lhs, Vector rhs);
// Outer product
Matrix outer(Vector lhs, Vector rhs);

// Multiply vector by scalar, scalar LHS
Vector operator*(double scalar, Vector rhs);
// Add scalar and assign
// Casts scalar to vector
void operator+=(Vector &lhs, double scalar);
// Add vector and assign
void operator+=(Vector &lhs, Vector rhs);
// Subtract scalar from vector and assign
// Casts scalar to vector
void operator-=(Vector &lhs, double scalar);
// Subtract vector and assign
void operator-=(Vector &lhs, Vector rhs);
// Multiply by scalar and assign
void operator*=(Vector &lhs, double scalar);
// Divide by scalar and assign
void operator/=(Vector &lhs, double scalar);
// Stream output
std::ostream& operator<<(std::ostream &os, const Vector &vec);
// Stream input
std::istream& operator>>(std::istream &is, Vector &vec);

/*******************************/
/* Matrix non-member functions */
/*******************************/

// Get determinant of 3x3 matrix
double det(Matrix mat);
// Get inverse of 3x3 matrix
Matrix inverse(Matrix mat);
// Get trace of 3x3 matrix
double trace(Matrix mat);
// Get eigenvalues of symmetric 3x3 matrix, in descending order
void eigenvalues(Matrix mat, double eigvals[3]);
// Get eigenvector of symmetrix 3x3 matrix that corresponds to given eigenvalue
Vector eigenvector(Matrix mat, double eigenvalue);
// Create rotation matrix about X axis
Matrix getRotMatX(double angle);
// Create rotation matrix about Y axis
Matrix getRotMatY(double angle);
// Create rotation matrix about Z axis
Matrix getRotMatZ(double angle);

// Multiply matrix by scalar, scalar LHS
Matrix operator*(double scalar, Matrix rhs);
// Multiply matrices and assign
void operator*=(Matrix &lhs, Matrix rhs);
// Multiply matrix by scalar and assign
void operator*=(Matrix &lhs, double scalar);
// Divide matrix by scalar and assign
void operator/=(Matrix &lhs, double scalar);
// Add matrices and assign
void operator+=(Matrix &lhs, Matrix rhs);
// Subtract matrices and assign
void operator-=(Matrix &lhs, Matrix rhs);
// Stream output
std::ostream& operator<<(std::ostream &os, const Matrix &vec);
// Stream input
std::istream& operator>>(std::istream &is, Matrix &vec);

/*****************************/
/* Matrix & Vector Constants */
/*****************************/

// Declare a single identity matrix, zero vector, etc
extern double identity[3][3];
extern double zero_vec[3];
extern double zero_mat[3][3];
extern const Matrix I;
extern const Vector zero_vector;
extern const Matrix zero_matrix;

#endif /* SRC_SPRING_MATRIX_MATH_H_ */
