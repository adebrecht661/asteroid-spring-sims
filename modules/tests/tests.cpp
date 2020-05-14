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
extern "C" {
#include "rebound.h"
}
#include "matrix_math.h"
#include "kepcart.h"
#include "springs.h"
#include "input_spring.h"
#include "stress.h"
#include "physics.h"
#include <gtest/gtest.h>

using namespace libconfig;
using std::vector;
using std::abs;

/* These globals are expected in other files */

// Global values
int num_springs = 0;	// Global numbers of springs
vector<spring> springs;	// Global spring array
int num_perts = 0;		// Number of perturbers
extern vector<stress_tensor> stresses;

// Global scales
double mass_scale, time_scale, length_scale, temp_scale, omega_scale, vel_scale,
		p_scale, L_scale, a_scale, F_scale, E_scale, dEdt_scale, P_scale;

/*********************************************/

/************/
/* Fixtures */
/************/

class SpringTest: public testing::Test {
protected:

	// Objects for fixture
	const int num_parts = 3;
	reb_simulation *n_body_sim = (reb_simulation*) malloc(
			sizeof(reb_simulation));
	spring spr;

	// Allocate objects, set properties
	void SetUp() override {
		n_body_sim->particles = (reb_particle*) malloc(
				sizeof(reb_particle) * num_parts);
		n_body_sim->N = num_parts;

		// Particle 1 position
		n_body_sim->particles[0].x = 0;
		n_body_sim->particles[0].y = 0;
		n_body_sim->particles[0].z = 0;

		// Particle 1 velocity
		n_body_sim->particles[0].vx = 0;
		n_body_sim->particles[0].vy = 0;
		n_body_sim->particles[0].vz = 0;

		// Particle 2 position
		n_body_sim->particles[1].x = 1;
		n_body_sim->particles[1].y = 2;
		n_body_sim->particles[1].z = 2;

		// Particle 2 velocity
		n_body_sim->particles[1].vx = 1;
		n_body_sim->particles[1].vy = 2;
		n_body_sim->particles[1].vz = 2;

		// Particle 3 position
		n_body_sim->particles[2].x = 2;
		n_body_sim->particles[2].y = 3;
		n_body_sim->particles[2].z = 6;

		// Particle 3 velocity
		n_body_sim->particles[2].vx = 2;
		n_body_sim->particles[2].vy = 3;
		n_body_sim->particles[2].vz = 6;

		// Particle masses
		n_body_sim->particles[0].m = 1;
		n_body_sim->particles[1].m = 2;
		n_body_sim->particles[2].m = 3;

		// Particles associated with spring
		spr.particle_1 = 0;
		spr.particle_2 = 1;

		// Spring default properties
		spr.rs0 = 1;
		spr.k = 11;
		spr.gamma = 35;
	}

	void TearDown() override {
		// Ensure vectors are reallocated (rather than just changing nominal size)
		vector<spring>().swap(springs);
		num_springs = 0;
		vector<stress_tensor>().swap(stresses);
		free(n_body_sim->particles);
		free(n_body_sim);
	}
};

class PhysicsTest: public SpringTest {
protected:
	void SetUp() override {

		SpringTest::SetUp();

		// Particle 1 acceleration
		n_body_sim->particles[0].ax = 0;
		n_body_sim->particles[0].ay = 0;
		n_body_sim->particles[0].az = 0;

		// Particle 2 acceleration
		n_body_sim->particles[1].ax = 1;
		n_body_sim->particles[1].ay = 2;
		n_body_sim->particles[1].az = 2;

		// Particle 3 acceleration
		n_body_sim->particles[2].ax = 2;
		n_body_sim->particles[2].ay = 3;
		n_body_sim->particles[2].az = 6;
	}

	void TearDown() override {
		SpringTest::TearDown();
	}
};

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

/****************/
/* Spring tests */
/****************/

TEST_F(SpringTest, Displacement) {

	Vector dx = spring_r(n_body_sim, spr);

	EXPECT_EQ(dx.len(), 3);

	EXPECT_EQ(dx, Vector( { -1, -2, -2 }));
}

TEST_F(SpringTest, AddSprings) {
	spring bad_spring;
	bad_spring.particle_1 = 0;
	bad_spring.particle_2 = 0;

	ASSERT_NO_THROW(add_spring_helper(spr));

	std::cout << "Expect warning: \n";
	EXPECT_EQ(add_spring(n_body_sim, 0, 1, spr), 0);

	springs.clear();
	ASSERT_ANY_THROW(add_spring_helper(spr));
	num_springs = 0;

	EXPECT_EQ(add_spring(n_body_sim, 0, 1, spr), 0);
	std::cout << "Expect warning: " << std::endl;
	EXPECT_EQ(add_spring(n_body_sim, 0, 1, spr), 0);
	EXPECT_EQ(add_spring(n_body_sim, 0, 2, spr), 1);

	EXPECT_ANY_THROW(add_spring(n_body_sim, 0, 0, bad_spring));
}

TEST_F(SpringTest, DelSprings) {

	ASSERT_EQ(add_spring(n_body_sim, 0, 1, spr), 0);

	springs.clear();
	EXPECT_ANY_THROW(del_spring(0));
	num_springs = 0;

	ASSERT_EQ(add_spring(n_body_sim, 0, 1, spr), 0);
	EXPECT_EQ(springs[0].particle_1, 0);
	EXPECT_EQ(springs[0].particle_2, 1);
	ASSERT_EQ(add_spring(n_body_sim, 0, 2, spr), 1);

	ASSERT_NO_THROW(del_spring(0));

	EXPECT_EQ(springs[0].particle_1, 0);
	EXPECT_EQ(springs[0].particle_2, 2);
	EXPECT_EQ(springs[0].rs0, 7);

	EXPECT_EQ(num_springs, 1);
}

TEST_F(SpringTest, Gammas) {
	ASSERT_NO_THROW(add_spring(n_body_sim, 0, 1, spr));

	set_gamma(5);

	EXPECT_EQ(springs[0].gamma, 5);

	divide_gamma(2);

	EXPECT_EQ(springs[0].gamma, 2.5);
}

TEST_F(SpringTest, ConnectSprings) {
	connect_springs_dist(n_body_sim, 5, 0, num_parts, spr);

	ASSERT_EQ(num_springs, 2);

	EXPECT_EQ(springs[0].particle_1, 0);
	EXPECT_EQ(springs[1].particle_1, 1);
}

TEST_F(SpringTest, KillSprings) {
	connect_springs_dist(n_body_sim, 5, 0, num_parts, spr);
	set_gamma(5);

	stresses.resize(num_springs);
	stresses[0].failing = true;
	stresses[1].failing = false;

	kill_springs(n_body_sim);
	EXPECT_EQ(springs[0].k, 0);
	EXPECT_EQ(springs[0].gamma, 0);
	EXPECT_EQ(springs[1].k, 11);
	EXPECT_EQ(springs[1].gamma, 5);
}

TEST_F(SpringTest, MeanLength) {
	EXPECT_TRUE(std::isnan(mean_spring_length()));

	add_spring(n_body_sim, 0, 1, spr);

	EXPECT_EQ(mean_spring_length(), 3);

	add_spring(n_body_sim, 0, 2, spr);

	EXPECT_EQ(mean_spring_length(), 5);
}

TEST_F(SpringTest, Midpoint) {
	add_spring(n_body_sim, 0, 2, spr);
	Vector x1 = spr_mid(n_body_sim, springs[0], Vector( { 0, 0, 0 }));
	Vector x2 = spr_mid(n_body_sim, springs[0], Vector( { 1, 2, 4 }));

	// Midpoint from origin is equidistant from particle at origin and particle at <2,3,6>
	EXPECT_EQ(x1, Vector( { 2, 3, 6 }) / 2);

	EXPECT_EQ(x2, Vector( { 0, -0.5, -1 }));
}

TEST_F(SpringTest, Strain) {
	spring zero_spring = spr;
	zero_spring.rs0 = 0;

	EXPECT_TRUE(std::isinf(strain(n_body_sim, zero_spring)));

	EXPECT_EQ(strain(n_body_sim, spr), 2);
}

TEST_F(SpringTest, Force) {
	reb_particle *particles = n_body_sim->particles;
	connect_springs_dist(n_body_sim, 100, 0, num_parts, spr);

	set_gamma(0);

	// Yay nice numbers
	particles[0].x = 1;
	particles[0].y = 1;
	particles[0].z = 4;

	Vector ans_len = spring_r(n_body_sim, springs[1]);
	Vector ans_len_hat = ans_len / ans_len.len();

	Vector test_force1 = spring_i_force_undamped(n_body_sim, 1);
	Vector test_force2 = spring_i_force(n_body_sim, 1);
	spring_forces(n_body_sim);
	Vector ans_force = 44 * ans_len_hat;
	Vector ans_accel0 = (ans_force + spring_i_force(n_body_sim, 0))
			/ particles[0].m;
	Vector ans_accel2 = (-ans_force - spring_i_force(n_body_sim, 2))
			/ particles[2].m;

	EXPECT_NEAR(test_force1.getX(), ans_force.getX(),
			abs(0.000001 * ans_force.getX()));
	EXPECT_NEAR(test_force1.getY(), ans_force.getY(),
			abs(0.000001 * ans_force.getY()));
	EXPECT_NEAR(test_force1.getZ(), ans_force.getZ(),
			abs(0.000001 * ans_force.getZ()));

	EXPECT_NEAR(test_force2.getX(), ans_force.getX(),
			abs(0.000001 * ans_force.getX()));
	EXPECT_NEAR(test_force2.getY(), ans_force.getY(),
			abs(0.000001 * ans_force.getY()));
	EXPECT_NEAR(test_force2.getZ(), ans_force.getZ(),
			abs(0.000001 * ans_force.getZ()));

	EXPECT_NEAR(particles[0].ax, ans_accel0.getX(),
			abs(0.000001 * ans_accel0.getX()));
	EXPECT_NEAR(particles[0].ay, ans_accel0.getY(),
			abs(0.000001 * ans_accel0.getY()));
	EXPECT_NEAR(particles[0].az, ans_accel0.getZ(),
			abs(0.000001 * ans_accel0.getZ()));

	EXPECT_NEAR(particles[2].ax, ans_accel2.getX(),
			abs(0.000001 * ans_accel2.getX()));
	EXPECT_NEAR(particles[2].ay, ans_accel2.getY(),
			abs(0.000001 * ans_accel2.getY()));
	EXPECT_NEAR(particles[2].az, ans_accel2.getZ(),
			abs(0.000001 * ans_accel2.getZ()));

	set_gamma(1);
	zero_accel(n_body_sim);
	spring_forces(n_body_sim);
	Vector test_force3 = spring_i_force(n_body_sim, 1);
	Vector ans_force2 = 44 * ans_len_hat
			- 0.75 * dot(Vector( { -2, -3, -6 }), ans_len_hat) / ans_len.len()
					* ans_len_hat;

	Vector ans_accel01 = (spring_i_force(n_body_sim, 1)
			+ spring_i_force(n_body_sim, 0)) / particles[0].m;
	Vector ans_accel21 = (-spring_i_force(n_body_sim, 1)
			- spring_i_force(n_body_sim, 2)) / particles[2].m;

	EXPECT_NEAR(test_force3.getX(), ans_force2.getX(),
			abs(0.000001 * ans_force2.getX()));
	EXPECT_NEAR(test_force3.getY(), ans_force2.getY(),
			abs(0.000001 * ans_force2.getY()));
	EXPECT_NEAR(test_force3.getZ(), ans_force2.getZ(),
			abs(0.000001 * ans_force2.getZ()));

	EXPECT_NEAR(particles[0].ax, ans_accel01.getX(),
			abs(0.000001 * ans_accel01.getX()));
	EXPECT_NEAR(particles[0].ay, ans_accel01.getY(),
			abs(0.000001 * ans_accel01.getY()));
	EXPECT_NEAR(particles[0].az, ans_accel01.getZ(),
			abs(0.000001 * ans_accel01.getZ()));

	EXPECT_NEAR(particles[2].ax, ans_accel21.getX(),
			abs(0.000001 * ans_accel21.getX()));
	EXPECT_NEAR(particles[2].ay, ans_accel21.getY(),
			abs(0.000001 * ans_accel21.getY()));
	EXPECT_NEAR(particles[2].az, ans_accel21.getZ(),
			abs(0.000001 * ans_accel21.getZ()));
}

TEST_F(SpringTest, AdjustProps) {
	connect_springs_dist(n_body_sim, 100, 0, num_parts, spr);

	adjust_spring_props(n_body_sim, 2, 5, 0, 3);

	EXPECT_EQ(springs[1].k, 2);
	EXPECT_EQ(springs[1].gamma, 5);

	EXPECT_EQ(springs[0].k, spr.k);
	EXPECT_EQ(springs[0].gamma, spr.gamma);
}

/*****************/
/* Physics tests */
/*****************/

TEST_F(PhysicsTest, ZeroAccel) {
	reb_particle *particles = n_body_sim->particles;

	ASSERT_EQ(particles[1].ax, 1);
	ASSERT_EQ(particles[1].ay, 2);
	ASSERT_EQ(particles[1].az, 2);

	zero_accel(n_body_sim);

	EXPECT_EQ(particles[0].ax, 0);
	EXPECT_EQ(particles[0].ay, 0);
	EXPECT_EQ(particles[0].az, 0);

	EXPECT_EQ(particles[1].ax, 0);
	EXPECT_EQ(particles[1].ay, 0);
	EXPECT_EQ(particles[1].az, 0);

	EXPECT_EQ(particles[2].ax, 0);
	EXPECT_EQ(particles[2].ay, 0);
	EXPECT_EQ(particles[2].az, 0);
}

TEST_F(PhysicsTest, MoveResolved) {
	reb_particle *particles = n_body_sim->particles;
	move_resolved(n_body_sim, -1, 1, 0, num_parts);

	EXPECT_EQ(particles[0].x, -1);
	EXPECT_EQ(particles[0].y, -1);
	EXPECT_EQ(particles[0].z, -1);

	EXPECT_EQ(particles[1].x, 0);
	EXPECT_EQ(particles[1].y, 1);
	EXPECT_EQ(particles[1].z, 1);

	EXPECT_EQ(particles[0].vx, 1);
	EXPECT_EQ(particles[0].vy, 1);
	EXPECT_EQ(particles[0].vz, 1);

	EXPECT_EQ(particles[1].vx, 2);
	EXPECT_EQ(particles[1].vy, 3);
	EXPECT_EQ(particles[1].vz, 3);
}

TEST_F(PhysicsTest, CoM) {
	reb_particle *particles = n_body_sim->particles;
	Vector orig_CoM = { 4. / 3., 13. / 6., 11. / 3. };

	EXPECT_EQ(compute_com(n_body_sim, 0, num_parts), orig_CoM);

	ASSERT_NO_THROW(subtract_com(n_body_sim, 0, num_parts));

	EXPECT_EQ(particles[0].x, -orig_CoM.getX());
	EXPECT_EQ(particles[0].y, -orig_CoM.getY());
	EXPECT_EQ(particles[0].z, -orig_CoM.getZ());

	EXPECT_EQ(particles[1].x, 1 - orig_CoM.getX());
	EXPECT_EQ(particles[1].y, 2 - orig_CoM.getY());
	EXPECT_EQ(particles[1].z, 2 - orig_CoM.getZ());

	EXPECT_EQ(particles[0].vx, 0);
	EXPECT_EQ(particles[0].vy, 0);
	EXPECT_EQ(particles[0].vz, 0);

	Vector new_CoM = compute_com(n_body_sim, 0, num_parts);

	EXPECT_NEAR(new_CoM.getX(), 0, 1e-15);
	EXPECT_NEAR(new_CoM.getY(), 0, 1e-15);
	EXPECT_NEAR(new_CoM.getZ(), 0, 1e-15);
}

TEST_F(PhysicsTest, CoV) {
	reb_particle *particles = n_body_sim->particles;
	Vector orig_CoV = { 5. / 3., 13. / 6., 11. / 3. };
	particles[1].vx = 2;

	EXPECT_EQ(compute_cov(n_body_sim, 0, num_parts), orig_CoV);

	ASSERT_NO_THROW(subtract_cov(n_body_sim, 0, num_parts));

	EXPECT_EQ(particles[0].vx, -orig_CoV.getX());
	EXPECT_EQ(particles[0].vy, -orig_CoV.getY());
	EXPECT_EQ(particles[0].vz, -orig_CoV.getZ());

	EXPECT_EQ(particles[1].vx, 2 - orig_CoV.getX());
	EXPECT_EQ(particles[1].vy, 2 - orig_CoV.getY());
	EXPECT_EQ(particles[1].vz, 2 - orig_CoV.getZ());

	EXPECT_EQ(particles[0].x, 0);
	EXPECT_EQ(particles[0].y, 0);
	EXPECT_EQ(particles[0].z, 0);

	Vector new_CoV = compute_cov(n_body_sim, 0, num_parts);

	EXPECT_NEAR(new_CoV.getX(), 0, 1e-15);
	EXPECT_NEAR(new_CoV.getY(), 0, 1e-15);
	EXPECT_NEAR(new_CoV.getZ(), 0, 1e-15);
}

TEST_F(PhysicsTest, CenterAll) {
	reb_particle *particles = n_body_sim->particles;
	Vector orig_CoM = { 4. / 3., 13. / 6., 11. / 3. };
	ASSERT_EQ(orig_CoM, compute_com(n_body_sim, 0, num_parts));

	ASSERT_NO_THROW(center_sim(n_body_sim, 0, 2));

	EXPECT_EQ(particles[0].x, -2. / 3.);
	EXPECT_EQ(particles[0].y, -4. / 3.);
	EXPECT_EQ(particles[0].z, -4. / 3.);

	EXPECT_EQ(particles[2].x, 2 - 2. / 3.);
	EXPECT_EQ(particles[2].y, 3 - 4. / 3.);
	EXPECT_EQ(particles[2].z, 6 - 4. / 3.);
}

TEST_F(PhysicsTest, AdjustMass) {
	reb_particle *particles = n_body_sim->particles;
	particles[3].m = 0;

	ASSERT_NO_THROW(adjust_mass_side(n_body_sim, 1. / 2., 1.5));

	EXPECT_EQ(particles[0].m, 2. / 9.);
	EXPECT_EQ(particles[1].m, 4. / 9.);
	EXPECT_EQ(particles[2].m, 1. / 3.);
}

TEST_F(PhysicsTest, Inertia) {
	Matrix inertia_mat = mom_inertia(n_body_sim, 0, num_parts);

	EXPECT_DOUBLE_EQ(inertia_mat.getXX(), 253. / 6.);
	EXPECT_DOUBLE_EQ(inertia_mat.getYY(), 116. / 3.);
	EXPECT_DOUBLE_EQ(inertia_mat.getZZ(), 61. / 6.);
	EXPECT_DOUBLE_EQ(inertia_mat.getXY(), -14. / 3.);
	EXPECT_DOUBLE_EQ(inertia_mat.getXZ(), -32. / 3.);
	EXPECT_DOUBLE_EQ(inertia_mat.getYZ(), -43. / 3.);
}

TEST_F(PhysicsTest, GPE) {
	double G = 6.67428e-8 / F_scale / pow(length_scale, 2.0)
			* pow(mass_scale, 2.0);
	n_body_sim->G = G;

	EXPECT_EQ(grav_potential_energy(n_body_sim, 0, num_parts),
			-G * (23. / 21. + sqrt(2.)));
}

TEST_F(PhysicsTest, MeasureAngularMomentum) {
	// Add perturbing particle
	num_perts = 1;
	reb_particle temp_particles[num_parts];
	for (int i = 0; i < num_parts; i++) {
		temp_particles[i] = n_body_sim->particles[i];
	}

	n_body_sim->particles = (reb_particle*) malloc(
			sizeof(reb_particle) * (num_parts + num_perts));
	n_body_sim->N = num_parts + num_perts;
	for (int i = 0; i < num_parts; i++) {
		n_body_sim->particles[i] = temp_particles[i];
	}

	// Perturber location
	n_body_sim->particles[3].x = 12;
	n_body_sim->particles[3].y = 16;
	n_body_sim->particles[3].z = 21;

	// Perturber velocity
	n_body_sim->particles[3].vx = 12;
	n_body_sim->particles[3].vy = 16;
	n_body_sim->particles[3].vz = 21;

	// Perturber mass
	n_body_sim->particles[3].m = 200;

	Vector test_L = measure_L(n_body_sim, 0, num_parts);
	Vector test_Lo = compute_Lorb(n_body_sim, 0, num_parts);
	Vector test_L_origin = measure_L_origin(n_body_sim, 0, num_parts);

	EXPECT_NEAR(test_L.getX(), 0., 1e-15);
	EXPECT_NEAR(test_L.getY(), 0., 1e-15);
	EXPECT_NEAR(test_L.getZ(), 0., 1e-15);

	EXPECT_NEAR(test_Lo.getX(), 0., 1e-13);
	EXPECT_NEAR(test_Lo.getY(), 0., 1e-13);
	EXPECT_NEAR(test_Lo.getZ(), 0., 1e-13);

	EXPECT_NEAR(test_L_origin.getX(), 0., 1e-13);
	EXPECT_NEAR(test_L_origin.getY(), 0., 1e-13);
	EXPECT_NEAR(test_L_origin.getZ(), 0., 1e-13);
}

TEST_F(PhysicsTest, Spin) {
	Vector ans_omega = { 1, 12, 12 };
	Vector orig_CoV = compute_cov(n_body_sim, 0, num_parts);

	ASSERT_NO_THROW(spin_body(n_body_sim, 0, num_parts, ans_omega));

	Vector test_omega = body_spin(n_body_sim, 0, num_parts);
	Vector test_CoV = compute_cov(n_body_sim, 0, num_parts);

	EXPECT_NEAR(test_omega.getX(), ans_omega.getX(), 1e-13);
	EXPECT_NEAR(test_omega.getY(), ans_omega.getY(), 1e-13);
	EXPECT_NEAR(test_omega.getZ(), ans_omega.getZ(), 1e-13);

	EXPECT_NEAR(test_CoV.getX(), orig_CoV.getX(), 1e-15);
	EXPECT_NEAR(test_CoV.getY(), orig_CoV.getY(), 1e-15);
	EXPECT_NEAR(test_CoV.getZ(), orig_CoV.getZ(), 1e-15);
}

TEST_F(PhysicsTest, RotateBody) {
	reb_particle *particles = n_body_sim->particles;
	rotate_body(n_body_sim, 0, num_parts, 1, 1, 1);

	EXPECT_NEAR(particles[0].x, 0.375236, 0.000001);
	EXPECT_NEAR(particles[0].y, 4.0924, 0.0001);
	EXPECT_NEAR(particles[0].z, -0.243612, 0.000001);
}

TEST_F(PhysicsTest, RotateOrigin) {
	reb_particle *particles = n_body_sim->particles;
	rotate_origin(n_body_sim, 0, num_parts, 1, 1, 1);

	EXPECT_NEAR(particles[1].x, -0.0750932, 0.0000001);
	EXPECT_NEAR(particles[1].y, -1.30969, 0.00001);
	EXPECT_NEAR(particles[1].z, 2.69798, 0.00001);
}

TEST_F(PhysicsTest, RotateToPrinciple) {
	Matrix orig_inertia = mom_inertia(n_body_sim, 0, num_parts);
	double orig_eigs[3];
	eigenvalues(orig_inertia, orig_eigs);

	rotate_to_principal(n_body_sim, 0, num_parts);

	Matrix test_inertia = mom_inertia(n_body_sim, 0, num_parts);
	double test_eigs[3];
	eigenvalues(test_inertia, test_eigs);

	EXPECT_NEAR(test_inertia.getXY(), 0, 1e-13);
	EXPECT_NEAR(test_inertia.getXZ(), 0, 1e-13);
	EXPECT_NEAR(test_inertia.getYX(), 0, 1e-13);
	EXPECT_NEAR(test_inertia.getYZ(), 0, 1e-13);
	EXPECT_NEAR(test_inertia.getZX(), 0, 1e-13);
	EXPECT_NEAR(test_inertia.getZY(), 0, 1e-13);

	EXPECT_NEAR(test_eigs[0], orig_eigs[0], 1e-13);
	EXPECT_NEAR(test_eigs[1], orig_eigs[1], 1e-13);
	EXPECT_NEAR(test_eigs[2], orig_eigs[2], 1e-13);
}

TEST_F(PhysicsTest, SpringPower) {
	reb_particle *particles = n_body_sim->particles;
	connect_springs_dist(n_body_sim, 100, 0, num_parts, spr);
	set_gamma(1);

	// Yay nice numbers
	particles[0].x = 1;
	particles[0].y = 1;
	particles[0].z = 4;

	ASSERT_EQ(springs[1].particle_1, 0);
	ASSERT_EQ(springs[1].particle_2, 2);

	Vector ans_len = spring_r(n_body_sim, springs[1]);
	Vector ans_len_hat = ans_len / ans_len.len();
	double ans_pow = 0.75
			* pow(dot(Vector( { -2, -3, -6 }), ans_len_hat) / ans_len.len(),
					2.);
	EXPECT_NEAR(dEdt(n_body_sim, springs[1]), ans_pow, 1e-5);

	double ans_total_pow = dEdt(n_body_sim, springs[0])
			+ dEdt(n_body_sim, springs[1]) + dEdt(n_body_sim, springs[2]);

	EXPECT_EQ(dEdt_total(n_body_sim), ans_total_pow);
}

TEST_F(PhysicsTest, SPE) {
	reb_particle *particles = n_body_sim->particles;
	connect_springs_dist(n_body_sim, 100, 0, num_parts, spr);

	// Yay nice numbers
	particles[0].x = 1;
	particles[0].y = 1;
	particles[0].z = 4;

	EXPECT_NEAR(spring_potential_energy(n_body_sim), 88+3.20975, 0.0001);
}

TEST_F(PhysicsTest, NonTransKE) {
	EXPECT_EQ(compute_rot_kin(n_body_sim, 0, num_parts), 22.75);
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
