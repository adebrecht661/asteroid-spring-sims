/**
 * @file        spring.h
 * @brief       springs 
 * @author      Alice Quillen
 */

#ifndef _SPRING_H
#define _SPRING_H

#include "matrix_math.h"

// Spring structure
typedef struct spring {
	double k;		// Spring constant
	double rs0;		// Rest length
	double gamma; 	// Damping coefficient
	double k_heat;  // Heat diffusion coefficient (actually thermal conductance????? W/K vs W/m*K)
	int particle_1;	// Attached to particle 1
	int particle_2;	// Attached to particle 2
} spring;

extern struct spring* springs;	// Array to store springs
extern int num_springs;			// Number of springs
extern int num_perts;			// Number of perturbing point masses
extern double min_sep;			// Minimum separation between particles in simulation
extern double max_spring_dist;	// Maximum separation for particles to be connected by spring
extern double gamma_all;		// Gamma of all springs

/*********************/
/* Spring operations */
/*********************/

// Delete spring i
void del_spring(struct reb_simulation *const n_body_sim, int i);
// Add a spring with properties of spr between particle 1 and particle 2. Returns index of spring connecting these particles.
int add_spring(struct reb_simulation *const n_body_sim, int particle_1,
		int particle_2, spring spr);
// Helper function to add a spring at end of current array, expand if required (no sanity checking)
void add_spring_helper(spring spr);
// Connects springs between all particles closer than max_dist in index range i_low -> i_high-1
void connect_springs_dist(struct reb_simulation *const n_body_sim,
		double max_dist, int i_low, int i_high, spring spr);
// Connects springs between all particles closer than max_dist and farther than min_dist, no more than nodemax springs per particle
void connect_springs_dist_nodemax(struct reb_simulation *const n_body_sim,
		double min_dist, double max_dist, int i_low, int i_high, spring spr,
		int nodemax);
// Set damping coefficient of all springs
void set_gamma(double new_gamma);
// Modify spring constant, damping coefficient, heat diffusion coefficient of springs with midpoints inside [r_min, r_max]
void adjust_spring_props(struct reb_simulation *const n_body_sim, double new_k,
		double new_gamma, double new_k_heat, double r_min, double r_max);
// Kill springs that have failed
void kill_springs(struct reb_simulation *const n_body_sim);
// Make a binary system with a spring connecting them, orbiting with vector omega
void make_binary_spring(struct reb_simulation *const n_body_sim, double m1,
		double m2, double sep, Vector omega, spring spring_vals);

/*********************/
/* Spring properties */
/*********************/

// Get length of passed spring
double spring_length(struct reb_simulation *const n_body_sim, spring spr);
// Get mean rest length of springs
double mean_spring_length(struct reb_simulation *const r);
// Mean spring constant
double mean_ks(struct reb_simulation *const n_body_sim);
// Compute spring midpoint location from arbitrary center in Cartesian coordinates
Vector spr_mid(struct reb_simulation *const n_body_sim, spring spr,
		Vector center);
// Returns strain on given spring
double strain(struct reb_simulation *const n_body_sim, spring spr);
// Apply spring forces to simulation
void spring_forces(struct reb_simulation *n_body_sim);
// Get force from spring part_1 (with damping)
Vector spring_i_force(struct reb_simulation *const n_body_sim, int i);
// Get force from spring part_1 (without damping)
Vector spring_i_force_undamped(struct reb_simulation *const n_body_sim, int i);

#endif // _SPRING_H

