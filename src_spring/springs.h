/**
 * @file        spring.h
 * @brief       springs 
 * @author      Alice Quillen
 */

#ifndef _SPRING_H
#define _SPRING_H

#include "matrix_math.h"

// Spring structure
struct spring {
	double k;		// Spring constant
	double rs0;		// Rest length
	double gamma; 	// Damping coefficient
	double k_heat;// Heat diffusion coefficient (actually thermal conductance????? W/K vs W/m*K)
	int particle_1;	// Attached to particle 1
	int particle_2;	// Attached to particle 2
};

/*********************/
/* Spring operations */
/*********************/

// Delete spring i
void del_spring(reb_simulation *const n_body_sim, int i);
// Add a spring with properties of spr between particle 1 and particle 2. Returns index of spring connecting these particles.
int add_spring(reb_simulation *const n_body_sim, int particle_1, int particle_2,
		spring spr);
// Helper function to add a spring at end of current array, expand if required (no sanity checking)
void add_spring_helper(spring spr);
// Connects springs between all particles closer than max_dist in index range i_low -> i_high-1
void connect_springs_dist(reb_simulation *const n_body_sim, double max_dist,
		int i_low, int i_high, spring spr);
// Connects springs between all particles closer than max_dist and farther than min_dist, no more than nodemax springs per particle
void connect_springs_dist_nodemax(reb_simulation *const n_body_sim,
		double min_dist, double max_dist, int i_low, int i_high, spring spr,
		int nodemax);
// Set damping coefficient of all springs
void set_gamma(double new_gamma);
// Modify spring constant, damping coefficient, heat diffusion coefficient of springs with midpoints inside [r_min, r_max]
void adjust_spring_props(reb_simulation *const n_body_sim, double new_k,
		double new_gamma, double new_k_heat, double r_min, double r_max);
// Adjust spring properties either inside or outside designated ellipse (center x0, semiaxes a,b,c)
void adjust_spring_props_ellipsoid(reb_simulation *const n_body_sim,
		double new_k, double new_gamma, double new_k_heat, double a, double b,
		double c, Vector x0, bool inside);
// Adjust spring properties either inside or outside designated ellipse (center x0, semiaxes a,b,c) dependent on angle in xy plane
void adjust_spring_props_ellipsoid_phase(reb_simulation *const n_body_sim,
		double k_amp, double gamma_amp, double k_heat_amp, double freq,
		double phi_0, double a, double b, double c, Vector x0, bool inside);
// Kill springs that have failed
void kill_springs(reb_simulation *const n_body_sim);
// Make a binary system with a spring connecting them, orbiting with vector omega
void make_binary_spring(reb_simulation *const n_body_sim, double m1, double m2,
		double sep, Vector omega, spring spring_vals);

/*********************/
/* Spring properties */
/*********************/

// Get length vector of passed spring (r_ij)
Vector spring_r(reb_simulation *const n_body_sim, spring spr);
// Get mean rest length of springs
double mean_spring_length();
// Mean spring constant
double mean_ks(reb_simulation *const n_body_sim);
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

#endif // _SPRING_H

