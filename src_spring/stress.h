/*
 * stress.h
 *
 * Operations relating to tensile strength and stress on nodes
 *
 *  Created on: Mar 26, 2020
 *      Author: alex
 */

#ifndef SRC_SPRING_STRESS_H_
#define SRC_SPRING_STRESS_H_

#include "matrix_math.h"

// Stress tensor for each node
struct stress_tensor {
	// Stress tensor
	Matrix stress;

	// Eigenvalues of stress tensor, from largest to smallest
	double eigs[3];

	double max_force;  // Max force on a spring
	int max_force_spring;  // Index of spring giving max force
	bool failing;  // Is there material failure?
};

// Update stress tensor for each node
void update_stress(reb_simulation *const n_body_sim);
// Mark nodes that have surpassed tensile strength (internal and surface nodes differ in strength)
int mark_failed_nodes(reb_simulation *const n_body_sim, double mass_div,
		double tens_str_int, double tens_str_surf);
// Calculate Young's modulus of springs with midpoints inside radius range [r_min, r_max] from center of mass of particles in set [i_low, i_high)
double Young_mesh(reb_simulation *const n_body_sim, int i_low,
		int i_high, double r_min, double r_max);
// Calculate Young's modulus of complete simulation
// Caution: Assumes radius = 1
double Young_full_mesh();

#endif /* SRC_SPRING_STRESS_H_ */
