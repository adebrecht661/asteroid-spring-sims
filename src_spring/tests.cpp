#ifdef __cplusplus
# 	ifdef __GNUC__
#		define restrict __restrict__
#	else
#		define restrict
#	endif
#endif
/*
 * tests.cpp
 *
 *  Created on: Mar 31, 2020
 *      Author: alex
 */

#include <cmath>
#include <vector>
#include "matrix_math.h"
#include "springs.h"
#include <gtest/gtest.h>

using std::vector;

/* These globals are expected in other files */

// Global values
int num_springs = 0;	// Global numbers of springs
vector<spring> springs;	// Global spring array
int num_perts = 0;		// Number of perturbers

// Global scales
double mass_scale, time_scale, length_scale, temp_scale, omega_scale, vel_scale,
		p_scale, L_scale, a_scale, F_scale, E_scale, dEdt_scale, P_scale;

/*********************************************/

/*********************/
/* Vector-only Tests */
/*********************/

TEST(VectorTest, ScalarConstructor) {
	Vector testvec(4);

	EXPECT_EQ(testvec.getX(), 4);
	EXPECT_EQ(testvec.getY(), 4);
	EXPECT_EQ(testvec.getZ(), 4);
}

TEST(VectorTest, ArrayConstructor) {
	Vector testvec(zero_vec);

	EXPECT_EQ(testvec.getX(), 0);
	EXPECT_EQ(testvec.getY(), 0);
	EXPECT_EQ(testvec.getZ(), 0);
}

TEST(VectorTest, ListConstructor) {
	Vector testvec = { 1, 2, 3 };

	EXPECT_EQ(testvec.getX(), 1);
	EXPECT_EQ(testvec.getY(), 2);
	EXPECT_EQ(testvec.getZ(), 3);
}

// Not sure this does what I want it to
TEST(VectorTest, CopyConstructor) {
	Vector testvec;
	{
		Vector smallscopevec = { 1, 2, 3 };
		testvec = smallscopevec;
	}

	EXPECT_EQ(testvec.getX(), 1);
	EXPECT_EQ(testvec.getY(), 2);
	EXPECT_EQ(testvec.getZ(), 3);
}

TEST(VectorTest, Addition) {
	Vector vec1(zero_vec);
	Vector vec2 = { 1, 2, 3 };
	Vector vec3 = vec1 + vec2;
	Vector vec4 = vec2 + 1;
	vec2 += 2;
	vec1 += vec4;

	EXPECT_EQ(vec1, vec4);

	EXPECT_EQ(vec2.getX(), 3);
	EXPECT_EQ(vec2.getY(), 4);
	EXPECT_EQ(vec2.getZ(), 5);

	EXPECT_EQ(vec3.getX(), 1);
	EXPECT_EQ(vec3.getY(), 2);
	EXPECT_EQ(vec3.getZ(), 3);

	EXPECT_EQ(vec4.getX(), 2);
	EXPECT_EQ(vec4.getY(), 3);
	EXPECT_EQ(vec4.getZ(), 4);
}

TEST(VectorTest, Negation) {
	Vector testvec = { 1, 2, 3 };

	EXPECT_EQ(-testvec.getX(), -1);
	EXPECT_EQ(-testvec.getY(), -2);
	EXPECT_EQ(-testvec.getZ(), -3);
}

TEST(VectorTest, Subtraction) {
	Vector vec1(zero_vec);
	Vector vec2 = { 1, 2, 3 };
	Vector vec3 = vec1 - vec2;
	Vector vec4 = vec2 - 1;
	Vector vec5 = 2 - vec2;
	vec2 -= 2;
	vec1 -= vec4;

	EXPECT_EQ(-vec1, vec4);

	EXPECT_EQ(vec2.getX(), -1);
	EXPECT_EQ(vec2.getY(), 0);
	EXPECT_EQ(vec2.getZ(), 1);

	EXPECT_EQ(vec3.getX(), -1);
	EXPECT_EQ(vec3.getY(), -2);
	EXPECT_EQ(vec3.getZ(), -3);

	EXPECT_EQ(vec4.getX(), 0);
	EXPECT_EQ(vec4.getY(), 1);
	EXPECT_EQ(vec4.getZ(), 2);

	EXPECT_EQ(vec5.getX(), 1);
	EXPECT_EQ(vec5.getY(), 0);
	EXPECT_EQ(vec5.getZ(), -1);
}

TEST(VectorTest, BasicMultiplication) {
	Vector vec1 = { 1, 2, 3 };
	Vector vec2 = vec1 * 2;
	Vector vec3 = 2 * vec1;
	Vector vec4 = { 2, 3, 4 };
	Vector vec5 = vec1 * vec4; // element by element
	vec4 *= 3;

	EXPECT_EQ(vec2.getX(), 2);
	EXPECT_EQ(vec2.getY(), 4);
	EXPECT_EQ(vec2.getZ(), 6);

	EXPECT_EQ(vec3.getX(), 2);
	EXPECT_EQ(vec3.getY(), 4);
	EXPECT_EQ(vec3.getZ(), 6);

	EXPECT_EQ(vec5.getX(), 2);
	EXPECT_EQ(vec5.getY(), 6);
	EXPECT_EQ(vec5.getZ(), 12);

	EXPECT_EQ(vec4.getX(), 6);
	EXPECT_EQ(vec4.getY(), 9);
	EXPECT_EQ(vec4.getZ(), 12);
}

TEST(VectorTest, VectorMultiplication) {
	Vector vec1 = { 1, 2, 3 };
	Vector vec2 = { 3, 2, 1 };

	EXPECT_EQ(dot(vec1, vec2), 10);

	EXPECT_EQ(cross(vec1, vec2), Vector( { -4, 8, -4 }));

	EXPECT_EQ(cross(vec1, vec2), -cross(vec2, vec1));
}

TEST(VectorTest, BasicDivision) {
	Vector vec1 = { 1, 2, 3 };
	Vector vec2 = vec1 / 2;
	vec1 /= 3;

	EXPECT_EQ(vec1.getX(), 1. / 3.);
	EXPECT_EQ(vec1.getY(), 2. / 3.);
	EXPECT_EQ(vec1.getZ(), 1);

	EXPECT_EQ(vec2.getX(), 1. / 2.);
	EXPECT_EQ(vec2.getY(), 1);
	EXPECT_EQ(vec2.getZ(), 3. / 2.);
}

TEST(VectorTest, Length) {
	Vector testvec = { 1, 2, 3 };

	EXPECT_EQ(zero_vector.len(), 0.);

	EXPECT_EQ(testvec.len(), sqrt(14));
}

/*********************/
/* Matrix-only Tests */
/*********************/

TEST(MatrixTest, ArrayConstructor) {
	Matrix testmat(zero_mat);

	EXPECT_EQ(testmat.getXX(), 0);
	EXPECT_EQ(testmat.getXY(), 0);
	EXPECT_EQ(testmat.getXZ(), 0);
	EXPECT_EQ(testmat.getYX(), 0);
	EXPECT_EQ(testmat.getYY(), 0);
	EXPECT_EQ(testmat.getYZ(), 0);
	EXPECT_EQ(testmat.getZX(), 0);
	EXPECT_EQ(testmat.getZY(), 0);
	EXPECT_EQ(testmat.getZZ(), 0);
}

TEST(MatrixTest, ListOfListConstructor) {
	Matrix testmat = { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } };

	EXPECT_EQ(testmat.getXX(), 1);
	EXPECT_EQ(testmat.getXY(), 2);
	EXPECT_EQ(testmat.getXZ(), 3);
	EXPECT_EQ(testmat.getYX(), 4);
	EXPECT_EQ(testmat.getYY(), 5);
	EXPECT_EQ(testmat.getYZ(), 6);
	EXPECT_EQ(testmat.getZX(), 7);
	EXPECT_EQ(testmat.getZY(), 8);
	EXPECT_EQ(testmat.getZZ(), 9);
}

TEST(MatrixTest, ListOfVectorConstructor) {
	Matrix testmat = { Vector( { 1, 2, 3 }), Vector( { 4, 5, 6 }), Vector( { 7,
			8, 9 }) };

	EXPECT_EQ(testmat.getXX(), 1);
	EXPECT_EQ(testmat.getXY(), 2);
	EXPECT_EQ(testmat.getXZ(), 3);
	EXPECT_EQ(testmat.getYX(), 4);
	EXPECT_EQ(testmat.getYY(), 5);
	EXPECT_EQ(testmat.getYZ(), 6);
	EXPECT_EQ(testmat.getZX(), 7);
	EXPECT_EQ(testmat.getZY(), 8);
	EXPECT_EQ(testmat.getZZ(), 9);
}

// Not sure this does what I want it to
TEST(MatrixTest, CopyConstructor) {
	Matrix testmat;
	{
		Matrix smallscopemat = { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } };
		testmat = smallscopemat;
	}

	EXPECT_EQ(testmat.getXX(), 1);
	EXPECT_EQ(testmat.getXY(), 2);
	EXPECT_EQ(testmat.getXZ(), 3);
	EXPECT_EQ(testmat.getYX(), 4);
	EXPECT_EQ(testmat.getYY(), 5);
	EXPECT_EQ(testmat.getYZ(), 6);
	EXPECT_EQ(testmat.getZX(), 7);
	EXPECT_EQ(testmat.getZY(), 8);
	EXPECT_EQ(testmat.getZZ(), 9);
}

TEST(MatrixTest, Addition) {
	Matrix mat1(zero_mat);
	Matrix mat2 = { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } };
	Matrix mat3 = mat1 + mat2;
	mat1 += mat2;

	EXPECT_EQ(mat1, mat3);

	EXPECT_EQ(mat3.getXX(), 1);
	EXPECT_EQ(mat3.getXY(), 2);
	EXPECT_EQ(mat3.getXZ(), 3);
	EXPECT_EQ(mat3.getYX(), 4);
	EXPECT_EQ(mat3.getYY(), 5);
	EXPECT_EQ(mat3.getYZ(), 6);
	EXPECT_EQ(mat3.getZX(), 7);
	EXPECT_EQ(mat3.getZY(), 8);
	EXPECT_EQ(mat3.getZZ(), 9);
}

TEST(MatrixTest, Subtraction) {
	Matrix mat1(zero_mat);
	Matrix mat2 = { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } };
	Matrix mat3 = mat1 - mat2;
	mat1 -= mat2;

	EXPECT_EQ(-mat1, mat2);

	EXPECT_EQ(mat3.getXX(), -1);
	EXPECT_EQ(mat3.getXY(), -2);
	EXPECT_EQ(mat3.getXZ(), -3);
	EXPECT_EQ(mat3.getYX(), -4);
	EXPECT_EQ(mat3.getYY(), -5);
	EXPECT_EQ(mat3.getYZ(), -6);
	EXPECT_EQ(mat3.getZX(), -7);
	EXPECT_EQ(mat3.getZY(), -8);
	EXPECT_EQ(mat3.getZZ(), -9);
}

TEST(MatrixTest, BasicMultiplication) {
	Matrix mat1 = { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } };
	Matrix mat2 = mat1 * 2;
	Matrix mat3 = 2 * mat1;
	mat1 *= 2;

	EXPECT_EQ(mat3.getXX(), 2);
	EXPECT_EQ(mat3.getXY(), 4);
	EXPECT_EQ(mat3.getXZ(), 6);
	EXPECT_EQ(mat3.getYX(), 8);
	EXPECT_EQ(mat3.getYY(), 10);
	EXPECT_EQ(mat3.getYZ(), 12);
	EXPECT_EQ(mat3.getZX(), 14);
	EXPECT_EQ(mat3.getZY(), 16);
	EXPECT_EQ(mat3.getZZ(), 18);

	EXPECT_EQ(mat1, mat2);

	EXPECT_EQ(mat2, mat3);
}

TEST(MatrixTest, Operations) {
	Matrix mat1 = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
	Matrix mat2 = { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } };
	Matrix mat3 = transpose(mat2);
	Matrix mat4 = { { 0, 2, 3 }, { 2, 0, 1 }, { 3, 1, 0 } };

	Matrix ans = { { -1. / 12., 1. / 4., 1. / 6. },
			{ 1. / 4., -3. / 4., 1. / 2. }, { 1. / 6., 1. / 2., -1. / 3. } };

	double eigs2[3], eigs4[3];
	EXPECT_ANY_THROW(eigenvalues(mat2, eigs2));
	EXPECT_ANY_THROW(eigenvector(mat2, eigs2[0]));
	eigenvalues(mat4, eigs4);
	Vector eig41 = eigenvector(mat4, eigs4[0]);
	Vector eig42 = eigenvector(mat4, eigs4[1]);
	Vector eig43 = eigenvector(mat4, eigs4[2]);

	EXPECT_TRUE(mat1.isSym());

	EXPECT_FALSE(mat2.isSym());

	EXPECT_TRUE(mat4.isSym());

	EXPECT_EQ(det(mat1), 1);

	EXPECT_EQ(det(mat2), 0);

	EXPECT_EQ(trace(mat1), 3);

	EXPECT_EQ(trace(mat4), 0);

	EXPECT_EQ(trace(mat2), 15);

	EXPECT_EQ(inverse(mat4), ans);

	EXPECT_ANY_THROW(inverse(mat2));

	EXPECT_EQ(transpose(mat1), mat1);

	EXPECT_EQ(mat3.getXX(), 1);
	EXPECT_EQ(mat3.getXY(), 4);
	EXPECT_EQ(mat3.getXZ(), 7);
	EXPECT_EQ(mat3.getYX(), 2);
	EXPECT_EQ(mat3.getYY(), 5);
	EXPECT_EQ(mat3.getYZ(), 8);
	EXPECT_EQ(mat3.getZX(), 3);
	EXPECT_EQ(mat3.getZY(), 6);
	EXPECT_EQ(mat3.getZZ(), 9);

	EXPECT_EQ(transpose(mat3), mat2);

	EXPECT_NEAR(eigs4[0], 4.11309, 0.00001);
	EXPECT_NEAR(eigs4[1], -0.91117, 0.00001);
	EXPECT_NEAR(eigs4[2], -3.20191, 0.00001);

	EXPECT_NEAR(eig41.getX(), 1.11006, 0.00001);
	EXPECT_NEAR(eig41.getY(), 0.78289, 0.00001);
	EXPECT_NEAR(eig41.getZ(), 1, 0.00001);

	EXPECT_NEAR(eig42.getX(), 0.23141, 0.00001);
	EXPECT_NEAR(eig42.getY(), -1.60543, 0.00001);
	EXPECT_NEAR(eig42.getZ(), 1, 0.00001);

	EXPECT_NEAR(eig43.getX(), -1.21648, 0.00001);
	EXPECT_NEAR(eig43.getY(), 0.44753, 0.00001);
	EXPECT_NEAR(eig43.getZ(), 1, 0.00001);
}

TEST(MatrixTest, RotationMatrices) {
	Matrix rotx = getRotMatX(M_PI / 3);
	Matrix roty = getRotMatY(M_PI / 3);
	Matrix rotz = getRotMatZ(M_PI / 3);

	Matrix ansx = { { 1, 0, 0 }, { 0, cos(M_PI / 3), -sin(M_PI / 3) }, { 0, sin(
	M_PI / 3), cos(M_PI / 3) } };
	Matrix ansy = { { cos(M_PI / 3), 0, sin(M_PI / 3) }, { 0, 1, 0 }, { -sin(
	M_PI / 3), 0, cos(M_PI / 3) } };
	Matrix ansz = { { cos(M_PI / 3), -sin(M_PI / 3), 0 }, { sin(M_PI / 3), cos(
	M_PI / 3), 0 }, { 0, 0, 1 } };

	EXPECT_EQ(rotx, ansx);

	EXPECT_EQ(roty, ansy);

	EXPECT_EQ(rotz, ansz);
}

TEST(MatrixTest, MatrixMultiplication) {
	Matrix mat1 = { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } };
	Matrix mat2 = { { 9, 8, 7 }, { 6, 5, 4 }, { 3, 2, 1 } };
	Matrix mat3 = mat1 * mat2;
	Matrix mat4 = mat2;
	mat2 *= mat1;

	Matrix ans1 = { { 30, 24, 18 }, { 84, 69, 54 }, { 138, 114, 90 } };
	Matrix ans2 = { { 90, 114, 138 }, { 54, 69, 84 }, { 18, 24, 30 } };

	EXPECT_EQ(ans1, mat3);

	EXPECT_EQ(ans2, mat2);

	EXPECT_EQ(transpose(mat3), transpose(mat4)*transpose(mat1));
}

TEST(MatrixVectorTest, Multiplication) {
	Vector vec1 = { 1, 2, 3 };
	Vector vec2 = { 3, 2, 1 };
	Matrix mat1 = outer(vec1, vec2);
	Matrix mat2 = outer(vec2, vec1);
	Matrix mat3 = { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } };
	Vector vec3 = mat3 * vec1;

	Matrix ans1 = { { 3, 2, 1 }, { 6, 4, 2 }, { 9, 6, 3 } };
	Vector ans2 = {14, 32, 50};

	EXPECT_EQ(mat1, ans1);

	EXPECT_EQ(transpose(mat1), mat2);

	EXPECT_EQ(vec3, ans2);
}
