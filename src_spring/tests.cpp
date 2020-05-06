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
	Vector testvec = {1,2,3};

	EXPECT_EQ(testvec.getX(), 1);
	EXPECT_EQ(testvec.getY(), 2);
	EXPECT_EQ(testvec.getZ(), 3);
}
