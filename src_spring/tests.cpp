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
#include "libconfig.h++"
#include "matrix_math.h"
#include "kepcart.h"
#include "springs.h"
#include "input_spring.h"
#include <gtest/gtest.h>

using namespace libconfig;
using std::vector;
using std::abs;

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

	// Make sure we're comparing to normalized eigenvectors
	Vector ans1 = { 1.11006, 0.78289, 1 };
	ans1 /= ans1.len();
	Vector ans2 = { 0.23141, -1.60543, 1 };
	ans2 /= ans2.len();
	Vector ans3 = { -1.21648, 0.44753, 1 };
	ans3 /= -ans3.len(); // The routine returns the negative normalized eigenvector for this one

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

	EXPECT_NEAR(eig41.getX(), ans1.getX(), 0.00001);
	EXPECT_NEAR(eig41.getY(), ans1.getY(), 0.00001);
	EXPECT_NEAR(eig41.getZ(), ans1.getZ(), 0.00001);

	EXPECT_NEAR(eig42.getX(), ans2.getX(), 0.00001);
	EXPECT_NEAR(eig42.getY(), ans2.getY(), 0.00001);
	EXPECT_NEAR(eig42.getZ(), ans2.getZ(), 0.00001);

	EXPECT_NEAR(eig43.getX(), ans3.getX(), 0.00001);
	EXPECT_NEAR(eig43.getY(), ans3.getY(), 0.00001);
	EXPECT_NEAR(eig43.getZ(), ans3.getZ(), 0.00001);
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

	EXPECT_EQ(transpose(mat3), transpose(mat4) * transpose(mat1));
}

/*****************************/
/* Mixed Matrix-Vector Tests */
/*****************************/

TEST(MatrixVectorTest, Multiplication) {
	Vector vec1 = { 1, 2, 3 };
	Vector vec2 = { 3, 2, 1 };
	Matrix mat1 = outer(vec1, vec2);
	Matrix mat2 = outer(vec2, vec1);
	Matrix mat3 = { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } };
	Vector vec3 = mat3 * vec1;

	Matrix ans1 = { { 3, 2, 1 }, { 6, 4, 2 }, { 9, 6, 3 } };
	Vector ans2 = { 14, 32, 50 };

	EXPECT_EQ(mat1, ans1);

	EXPECT_EQ(transpose(mat1), mat2);

	EXPECT_EQ(vec3, ans2);
}

/****************************************/
/* Keplerian-Cartesian Conversion Tests */
/****************************************/

TEST(KepCartTest, KepEqn) {
	double eccentricity1 = 0.5;
	double eccentricity2 = 1.5;
	double mean_anom = 27 * (M_PI / 180);

	EXPECT_NEAR(eccentric_anomaly(eccentricity1, mean_anom),
			48.43417991487915 * (M_PI/180), 0.001);

	EXPECT_ANY_THROW(eccentric_anomaly_hyperbolic(eccentricity1, mean_anom));

	EXPECT_NEAR(eccentric_anomaly_hyperbolic(eccentricity2, mean_anom),
			0.736899151737, 0.02);

	EXPECT_ANY_THROW(eccentric_anomaly(eccentricity2, mean_anom));
}

TEST(KepCartTest, Conversions) {
	PhaseState ans_state_units;
	ans_state_units.x = { 1.5e8, 2e7, 3e6 };
	ans_state_units.x /= length_scale;
	ans_state_units.v = { 5e5, 1e6, 1e6 };
	ans_state_units.v /= vel_scale;

	OrbitalElements ans_orbel_units;
	ans_orbel_units.a = 132145319.13522287 / length_scale;
	ans_orbel_units.e = 0.45113207590986515;
	ans_orbel_units.i = 0.81984644566308151337;
	ans_orbel_units.long_asc_node = 0.11398192262875016245;
	ans_orbel_units.arg_peri = 3.9982731061823910679 - 2 * M_PI;
	ans_orbel_units.mean_anom = 1.4718545303922809797;

	PhaseState ans_state_unitless;
	ans_state_unitless.x = { -6.03224695, 0.29119389, 5.36967604 };
	ans_state_unitless.v = { -0.11631944, -0.35758130, -0.07805506 };

	OrbitalElements ans_orbel_unitless;
	ans_orbel_unitless.a = 10;
	ans_orbel_unitless.e = 0.2;
	ans_orbel_unitless.i = 45 * (M_PI / 180);
	ans_orbel_unitless.long_asc_node = 60 * (M_PI / 180);
	ans_orbel_unitless.arg_peri = 90 * (M_PI / 180);
	ans_orbel_unitless.mean_anom = 0.229607;

	double G = 6.67428e-8 / F_scale / pow(length_scale, 2.0)
			* pow(mass_scale, 2.0);
	double M = 1; // 1 mass_scale

	OrbitalElements test_orbel_units = cart_to_kep(G * M, ans_state_units);
	OrbitalElements test_orbel_unitless = cart_to_kep(1, ans_state_unitless);

	PhaseState test_state_units = kep_to_cart(G * M, ans_orbel_units);
	PhaseState test_state_unitless = kep_to_cart(1, ans_orbel_unitless);

	EXPECT_NEAR(test_orbel_units.a, ans_orbel_units.a,
			abs(0.001 * ans_orbel_units.a));
	EXPECT_NEAR(test_orbel_units.e, ans_orbel_units.e,
			0.001 * ans_orbel_units.e);
	EXPECT_NEAR(test_orbel_units.i, ans_orbel_units.i,
			abs(0.0025 * ans_orbel_units.i));
	EXPECT_NEAR(test_orbel_units.long_asc_node, ans_orbel_units.long_asc_node,
			abs(0.001 * ans_orbel_units.long_asc_node));
	EXPECT_NEAR(test_orbel_units.arg_peri, ans_orbel_units.arg_peri,
			abs(0.002 * ans_orbel_units.arg_peri));
	EXPECT_NEAR(test_orbel_units.mean_anom, ans_orbel_units.mean_anom,
			abs(0.001 * ans_orbel_units.mean_anom));

	EXPECT_NEAR(test_orbel_unitless.a, ans_orbel_unitless.a,
			abs(0.001 * ans_orbel_unitless.a));
	EXPECT_NEAR(test_orbel_unitless.e, ans_orbel_unitless.e,
			0.001 * ans_orbel_unitless.e);
	EXPECT_NEAR(test_orbel_unitless.i, ans_orbel_unitless.i,
			abs(0.001 * ans_orbel_unitless.i));
	EXPECT_NEAR(test_orbel_unitless.long_asc_node,
			ans_orbel_unitless.long_asc_node,
			abs(0.001 * ans_orbel_unitless.long_asc_node));
	EXPECT_NEAR(test_orbel_unitless.arg_peri, ans_orbel_unitless.arg_peri,
			abs(0.001 * ans_orbel_unitless.arg_peri));
	EXPECT_NEAR(test_orbel_unitless.mean_anom, ans_orbel_unitless.mean_anom,
			abs(0.001 * ans_orbel_unitless.mean_anom));

	EXPECT_NEAR(test_state_units.x.getX(), ans_state_units.x.getX(),
			abs(0.002 * ans_state_units.x.getX()));
	EXPECT_NEAR(test_state_units.x.getY(), ans_state_units.x.getY(),
			abs(0.002 * ans_state_units.x.getY()));
	EXPECT_NEAR(test_state_units.x.getZ(), ans_state_units.x.getZ(),
			abs(0.002 * ans_state_units.x.getZ()));
	EXPECT_NEAR(test_state_units.v.getX(), ans_state_units.v.getX(),
			abs(0.002 * ans_state_units.v.getX()));
	EXPECT_NEAR(test_state_units.v.getY(), ans_state_units.v.getY(),
			abs(0.002 * ans_state_units.v.getY()));
	EXPECT_NEAR(test_state_units.v.getZ(), ans_state_units.v.getZ(),
			abs(0.002 * ans_state_units.v.getZ()));

	EXPECT_NEAR(test_state_unitless.x.getX(), ans_state_unitless.x.getX(),
			abs(0.001 * ans_state_unitless.x.getX()));
	EXPECT_NEAR(test_state_unitless.x.getY(), ans_state_unitless.x.getY(),
			abs(0.002 * ans_state_unitless.x.getY()));
	EXPECT_NEAR(test_state_unitless.x.getZ(), ans_state_unitless.x.getZ(),
			abs(0.001 * ans_state_unitless.x.getZ()));
	EXPECT_NEAR(test_state_unitless.v.getX(), ans_state_unitless.v.getX(),
			abs(0.001 * ans_state_unitless.v.getX()));
	EXPECT_NEAR(test_state_unitless.v.getY(), ans_state_unitless.v.getY(),
			abs(0.001 * ans_state_unitless.v.getY()));
	EXPECT_NEAR(test_state_unitless.v.getZ(), ans_state_unitless.v.getZ(),
			abs(0.001 * ans_state_unitless.v.getZ()));
}

class MyEnvironment: public testing::Environment {
public:
	void SetUp() override {
		Config cfg;
		ASSERT_NO_THROW(cfg.readFile("problem.cfg"));

		ASSERT_NO_THROW(read_scales(&cfg));
	}
};

int main(int argc, char **argv) {
	std::cout << "Running main from " << __FILE__ << "\n";

	testing::InitGoogleTest(&argc, argv);
	testing::AddGlobalTestEnvironment(new MyEnvironment);
	return RUN_ALL_TESTS();
}
