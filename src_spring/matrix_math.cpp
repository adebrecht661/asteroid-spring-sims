/*
 * Matrix.cpp
 *
 * Various matrix and vector methods
 *
 *  Created on: Mar 26, 2020
 *      Author: alex
 */

#include <cmath>
#include <limits>
#include <algorithm>
#include "matrix_math.h"

/*****************************/
/* Matrix & Vector Constants */
/*****************************/

// Declare a single identity matrix
double identity[3][3] = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
double zero_vec[3] = { 0, 0, 0 };
double zero_mat[3][3] = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };
const Matrix I = Matrix(identity);
const Vector zero_vector(zero_vec);
const Matrix zero_matrix(zero_mat);

/****************/
/* Constructors */
/****************/

// Null constructor
Vector::Vector() {
}

// Construct from scalar
Vector::Vector(double scalar) {
	for (int i = 0; i < 3; i++) {
		this->array[i] = scalar;
	}
}

// Construct from array
Vector::Vector(double array[3]) {
	for (int i = 0; i < 3; i++) {
		this->array[i] = array[i];
	}
}

// Construct from initializer list
Vector::Vector(std::initializer_list<double> list) {
	std::copy(list.begin(), list.end(), this->array);
}

// Copy constructor
Vector::Vector(const Vector &vector) {
	*this = vector;
}

// Null constructor
Matrix::Matrix() {
}

// Construct from array
Matrix::Matrix(double array[3][3]) {
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			this->array[i][j] = array[i][j];
		}
	}
}

// Construct from initializer list of initializer lists
Matrix::Matrix(std::initializer_list<std::initializer_list<double>> list) {
	int i = 0, j = 0;
	for (const auto &l : list) {
		for (const auto &v : l) {
			this->array[i][j] = v;
			j++;
		}
		j = 0;
		i++;
	}
}

// Construct from initializer list of (row) vectors
Matrix::Matrix(std::initializer_list<Vector> list) {
	int i = 0;
	for (const auto &l : list) {
		this->array[i][0] = l.getX();
		this->array[i][1] = l.getY();
		this->array[i][2] = l.getZ();
	}
	i++;
}

// Copy constructor
Matrix::Matrix(const Matrix &matrix) {
	*this = matrix;
}

/***********/
/* Getters */
/***********/

double Vector::getX() const {
	return this->array[0];
}

double Vector::getY() const {
	return this->array[1];
}

double Vector::getZ() const {
	return this->array[2];
}

double Matrix::getXX() const {
	return this->array[0][0];
}

double Matrix::getXY() const {
	return this->array[0][1];
}

double Matrix::getXZ() const {
	return this->array[0][2];
}

double Matrix::getYX() const {
	return this->array[1][0];
}

double Matrix::getYY() const {
	return this->array[1][1];
}

double Matrix::getYZ() const {
	return this->array[1][2];
}

double Matrix::getZX() const {
	return this->array[2][0];
}

double Matrix::getZY() const {
	return this->array[2][1];
}

double Matrix::getZZ() const {
	return this->array[2][2];
}

/******************/
/* Misc Functions */
/******************/

// Return length of vector
double Vector::len() {
	return sqrt(dot(*this, *this));
}

// Return dot product of vectors
double dot(Vector lhs, Vector rhs) {
	return lhs.getX() * rhs.getX() + lhs.getY() * rhs.getY()
			+ lhs.getZ() * rhs.getZ();
}

// Return cross product of vectors
Vector cross(Vector lhs, Vector rhs) {
	return Vector(
			{ -lhs.getY() * rhs.getZ() + lhs.getZ() * rhs.getY(), -lhs.getZ()
					* rhs.getX() + lhs.getX() * rhs.getZ(), -lhs.getX()
					* rhs.getY() + lhs.getY() * rhs.getX() });
}

// Return outer product of vectors
Matrix outer(Vector lhs, Vector rhs) {
	double mat[3][3];
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			mat[i][j] = lhs.array[i] * lhs.array[j];
		}
	}
	return Matrix(mat);
}

// Find normalized eigenvector for specified matrix and eigenvalue
Vector eigenvector(Matrix mat, double eigval) {

	// Check if matrix is symmetric
	if (!mat.isSym())
		return zero_vector;

	// Number of iterations
	int N_ITS = 20;

	// Compute inverse of A= (C - lambda I)
	mat -= eigval * I;
	Matrix inv_mat = inverse(mat); // inverse

	// Initialize eigenvector
	Vector eigvec( { 1.0, 1.0, 1.0 });

	// Iterate LHS of eigenvector eqn of inverse
	for (int i = 0; i < N_ITS; i++) {
		double length = (inv_mat * eigvec).len();
		// Normalize length of eigenvector at each loop
		eigvec /= length;
	}

	// Calculate LHS eigenvector eqn
	double diff1 = (mat * eigvec).len();

	// Save eigvec
	Vector eigvec1 = eigvec;

	// Reinitialize eigenvector
	eigvec = Vector( { 1.0, -0.5, 0.5 });

	// Calculate LHS of eigenvector eqn for inverse matrix with new initial eigenvector
	for (int i = 0; i < N_ITS; i++) {
		double length = (inv_mat * eigvec).len();
		// Normalize length of eigenvector at each loop
		eigvec /= length;
	}

	// Calculate new LHS of eigenvector eqn
	double diff2 = (mat * eigvec).len();

	// Set eigvec to correct value
	if (diff1 < diff2) {
		return eigvec1;
	} else {
		return eigvec;
	}
}

// Compute inverse of 3x3 symmetric matrix
Matrix inverse(Matrix mat) {
	return (1.0 / det(mat))
			* (0.5 * (pow(trace(mat), 2.0) - trace(mat * mat)) * I
					- mat * trace(mat) + mat * mat);
}

// Compute determinant of 3x3 matrix
double det(Matrix mat) {

	// Return determinant
	return mat.getXX() * (mat.getYY() * mat.getZZ() - mat.getZY() * mat.getYZ())
			+ mat.getXY()
					* (mat.getYZ() * mat.getZX() - mat.getYX() * mat.getZZ())
			+ mat.getXZ()
					* (mat.getZY() * mat.getYX() - mat.getZX() * mat.getYY());
}

// Return eigenvalues of symmetric matrix
// eigs[0] >= eigs[1] >= eigs[2]
void eigenvalues(Matrix mat, double eigs[3]) {

	// Return 0 eigenvalues if matrix isn't symmetric or other error occurs
	eigs[0] = eigs[1] = eigs[2] = 0;

	// Check if matrix is symmetric
	if (!mat.isSym()) {
		std::cerr
				<< "Matrix isn't symmetric. Can't be diagonalized with this routine."
				<< std::endl;
		return;
	}

	// Recipe from Wikipedia for eigenvalues of a symmetric matrix
	double p1 = mat.getXY() * mat.getYX() + mat.getYZ() * mat.getZY()
			+ mat.getXZ() * mat.getZX();

	// Matrix is diagonal.
	if (p1 == 0) {
		// Make sure eigenvalues are in order
		eigs[0] = std::max( { mat.getXX(), mat.getYY(), mat.getZZ() });
		eigs[2] = std::min( { mat.getXX(), mat.getYY(), mat.getZZ() });
		eigs[1] = trace(mat) - eigs[0] - eigs[2];
		// Matrix isn't diagonal
	} else {
		// Helper variables
		double q = trace(mat) / 3.0;
		double p2 = pow(mat.getXX() - q, 2.0) + pow(mat.getYY() - q, 2.0)
				+ pow(mat.getZZ() - q, 2.0) + 2.0 * p1;
		double p = sqrt(p2 / 6.0);

		Matrix B = (1.0 / p) * (mat - q * I);
		double r = 0.5 * det(B);

		// For a symmetric matrix, -1 <= r <= 1, but computation error can leave it slightly outside this range
		double phi = 0.0;
		if (r <= -1.0) {
			phi = M_PI / 3.0;
		} else if (r >= 1.0) {
			phi = 0;
		} else if ((r < 1.0) && (r > -1.0)) {
			phi = acos(r) / 3.0;
		}

		// The eigenvalues satisfy eig3 <= eig2 <= eig1
		eigs[0] = q + 2.0 * p * cos(phi);
		eigs[2] = q + 2.0 * p * cos(phi + (2.0 * M_PI / 3.0));
		eigs[1] = trace(mat) - eigs[0] - eigs[2]; //  since trace(A) = eig1 + eig2 + eig3
	}
}

// Return trace of matrix
double trace(Matrix mat) {
	return mat.getXX() + mat.getYY() + mat.getZZ();
}

// Check if 3x3 matrix is symmetric
bool Matrix::isSym() {
	return (this->array[0][1] == this->array[1][0]
			&& this->array[0][2] == this->array[2][0]
			&& this->array[1][2] == this->array[2][1]);
}

// Return rotation matrix for rotation about X axis
Matrix getRotMatX(double angle) {
	return Matrix(
			{ { 1, 0, 0 }, { 0, cos(angle), -sin(angle) }, { 0, sin(angle), cos(
					angle) } });
}

// Return rotation matrix for rotation about Y axis
Matrix getRotMatY(double angle) {
	return Matrix( { { cos(angle), 0, sin(angle) }, { 0, 1, 0 }, { -sin(angle),
			0, cos(angle) } });
}

// Return rotation matrix for rotation about Z axis
Matrix getRotMatZ(double angle) {
	return Matrix( { { cos(angle), -sin(angle), 0 },
			{ sin(angle), cos(angle), 0 }, { 0, 0, 1 } });
}

/********************/
/* Vector Operators */
/********************/

// Add two vectors
Vector Vector::operator+(Vector rhs) {
	Vector res;
	for (int i = 0; i < 3; i++) {
		res.array[i] = this->array[i] + rhs.array[i];
	}
	return res;
}

// Add a scalar to a vector
Vector Vector::operator+(double scalar) {
	Vector res;
	for (int i = 0; i < 3; i++) {
		res.array[i] = this->array[i] + scalar;
	}
	return res;
}

// Negate a vector (unary operation)
Vector Vector::operator-() {
	Vector res;
	for (int i = 0; i < 3; i++) {
		res.array[i] = -this->array[i];
	}
	return res;
}

// Subtract two vectors
Vector Vector::operator-(Vector rhs) {
	Vector res;
	for (int i = 0; i < 3; i++) {
		res.array[i] = this->array[i] - rhs.array[i];
	}
	return res;
}

// Subtract a scalar from a vector
Vector Vector::operator-(double scalar) {
	Vector res;
	for (int i = 0; i < 3; i++) {
		res.array[i] = this->array[i] - scalar;
	}
	return res;
}

// Multiply a vector by a scalar (scalar on rhs)
Vector Vector::operator*(double scalar) {
	Vector res;
	for (int i = 0; i < 3; i++) {
		res.array[i] = this->array[i] * scalar;
	}
	return res;
}

// Multiply two vectors (element by element)
Vector Vector::operator*(Vector rhs) {
	Vector res;
	for (int i = 0; i < 3; i++) {
		res.array[i] = this->array[i] * rhs.array[i];
	}
	return res;
}

// Divide a vector by a scalar
Vector Vector::operator/(double scalar) {
	Vector res;
	for (int i = 0; i < 3; i++) {
		res.array[i] = this->array[i] / scalar;
	}
	return res;
}

// Equality
bool Vector::operator==(Vector rhs) const {
	return (this->array[0] == rhs.array[0] && this->array[1] == rhs.array[1] && this->array[2] == rhs.array[2]);
}

// Multiply a matrix by a vector
Vector operator*(Matrix lhs, Vector rhs) {
	Vector res;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			res.array[i] += lhs.array[i][j] * rhs.array[j];
		}
	}
	return res;
}

// Subtract a vector from a scalar
Vector operator-(double scalar, Vector rhs) {
	return -(rhs - scalar);
}

// Multiply a vector by a scalar (scalar on lhs)
Vector operator*(double scalar, Vector rhs) {
	return rhs * scalar;
}

// Add scalar to vector with assignment
void operator+=(Vector &lhs, double scalar) {
	double scalar_vec[3] = { scalar, scalar, scalar };
	Vector rhs(scalar_vec);
	lhs = lhs + rhs;
}

// Add two vectors with assignment
void operator+=(Vector &lhs, Vector rhs) {
	lhs = lhs + rhs;
}

// Subtract scalar from vector with assignment
void operator-=(Vector &lhs, double scalar) {
	double scalar_vec[3] = { scalar, scalar, scalar };
	Vector rhs(scalar_vec);
	lhs = lhs - rhs;
}

// Subtract two vectors with assignment
void operator-=(Vector &lhs, Vector rhs) {
	lhs = lhs - rhs;
}

// Multiply vector by scalar with assignment
void operator*=(Vector &lhs, double scalar) {
	lhs = lhs * scalar;
}

// Divide vector by scalar with assignment
void operator/=(Vector &lhs, double scalar) {
	lhs = lhs / scalar;
}

// Stream output for nicely formatted vector
std::ostream& operator<<(std::ostream &os, const Vector &vec) {
	os << "<" << vec.getX() << ", " << vec.getY() << ", " << vec.getZ() << ">";
	return os;
}

/********************/
/* Matrix Operators */
/********************/

// Add two matrices
Matrix Matrix::operator+(Matrix rhs) {
	Matrix res;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			res.array[i][j] = this->array[i][j] + rhs.array[i][j];
		}
	}
	return res;
}

// Subtract two matrices
Matrix Matrix::operator-(Matrix rhs) {
	Matrix res;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			res.array[i][j] = this->array[i][j] - rhs.array[i][j];
		}
	}
	return res;
}

// Multiply matrix by scalar (scalar on rhs)
Matrix Matrix::operator*(double scalar) {
	Matrix res;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			res.array[i][j] = this->array[i][j] * scalar;
		}
	}
	return res;
}

// Divide matrix by scalar
Matrix Matrix::operator/(double scalar) {
	Matrix res;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			res.array[i][j] = this->array[i][j] / scalar;
		}
	}
	return res;
}

// Multiply matrix by scalar (scalar on lhs)
Matrix operator*(double scalar, Matrix rhs) {
	return rhs * scalar;
}

// Multiply two matrices
Matrix operator*(Matrix lhs, Matrix rhs) {
	Matrix res;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				res.array[i][j] += lhs.array[i][k] * rhs.array[k][j];
			}
		}
	}
	return res;
}

// Multiply matrix by scalar with assignment
void operator*=(Matrix &lhs, double scalar) {
	lhs = lhs * scalar;
}

// Divide matrix by scalar with assignment
void operator/=(Matrix &lhs, double scalar) {
	lhs = lhs / scalar;
}

// Add two matrices with assignment
void operator+=(Matrix &lhs, Matrix rhs) {
	lhs = lhs - rhs;
}

// Subtract two matrices with assignment
void operator-=(Matrix &lhs, Matrix rhs) {
	lhs = lhs - rhs;
}

// Stream output operator for nicely formatted matrix print
std::ostream& operator<<(std::ostream &os, const Matrix &mat) {
	os << "[" << mat.getXX() << ", " << mat.getXY() << ", " << mat.getXZ()
			<< "\n" << mat.getYX() << ", " << mat.getYY() << ", " << mat.getYZ()
			<< "\n" << mat.getZX() << ", " << mat.getZY() << ", " << mat.getZZ()
			<< "]";
	return os;
}
