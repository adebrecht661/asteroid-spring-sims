/*
 * tests.cpp
 *
 *  Created on: Mar 31, 2020
 *      Author: alex
 */

#include "matrix_math.h"
#include <gtest/gtest.h>

TEST(VectorTest, ScalarConstructor) {
	Vector testvec(4);

	EXPECT_EQ(testvec.getX(), 4);
	EXPECT_EQ(testvec.getY(), 4);
	EXPECT_EQ(testvec.getZ(), 4);
}
