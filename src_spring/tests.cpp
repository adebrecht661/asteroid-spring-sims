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
