/**
 * @file        spring.h
 * @brief       springs
 * @author      Alice Quillen
 */

#ifndef _SPRING_H
#define _SPRING_H

extern "C" {
#include "rebound.h"
}
#include "matrix_math.h"

// Spring structure
struct spring {
	double k;		// Spring constant
	double rs0;		// Rest length
	double gamma; 	// Damping coefficient
	int particle_1;	// Attached to particle 1
	int particle_2;	// Attached to particle 2
};

/*********************/
/* Spring operations */
/*********************/

// Delete spring i
void del_spring(int i);
// Add a spring with properties of spr between particle 1 and particle 2. Returns index of spring connecting these particles.
int add_spring(reb_simulation *const n_body_sim, int particle_1, int particle_2,
		spring spr);
// Helper function to add a spring at end of current array, expand if required (no sanity checking)
void add_spring_helper(spring spr);
// Connects springs between all particles closer than max_dist in index range i_low -> i_high-1
void connect_springs_dist(reb_simulation *const n_body_sim, double max_dist,
		int i_low, int i_high, spring spr);
// Set damping coefficient of all springs
void set_gamma(double new_gamma);
// Set damping coefficient of all springs
void divide_gamma(double gamma_fac);
// Modify spring constant, damping coefficient, heat diffusion coefficient of springs with midpoints inside [r_min, r_max]
void adjust_spring_props(reb_simulation *const n_body_sim, double new_k,
		double new_gamma, double r_min, double r_max);
// Kill springs that have failed
void kill_springs(reb_simulation *const n_body_sim);

/*********************/
/* Spring properties */
/*********************/

// Get length vector of passed spring (r_ij)
Vector spring_r(reb_simulation *const n_body_sim, spring spr);
// Get mean rest length of springs
double mean_spring_length();
// Compute spring midpoint location from arbitrary center in Cartesian coordinates
Vector spr_mid(reb_simulation *const n_body_sim, spring spr, Vector center);
// Returns strain on given spring
double strain(reb_simulation *const n_body_sim, spring spr);
// Apply spring forces to simulation
void spring_forces(reb_simulation *n_body_sim);
// Get force from spring i (with damping)
Vector spring_i_force(reb_simulation *const n_body_sim, int i);
// Get force from spring i (without damping)
Vector spring_i_force_undamped(reb_simulation *const n_body_sim, int i);

/*************/
/* Operators */
/*************/

// Stream output
std::ostream& operator<<(std::ostream &os, const spring &spr);
// Equality
bool operator==(const spring lhs, const spring rhs);

#endif // _SPRING_H
