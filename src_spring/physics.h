/*
 * physics.h
 *
 * Various physics-related routines
 *
 *  Created on: Apr 2, 2020
 *      Author: alex
 */

#ifndef SRC_SPRING_PHYSICS_H_
#define SRC_SPRING_PHYSICS_H_

#include "matrix_math.h"
#include "springs.h"

/*******************/
/* Center routines */
/*******************/

// Center of mass of particles in set [i_low, i_high)
Vector compute_com(reb_simulation *const n_body_sim, int i_low,
		int i_high);
// Center of velocity of particles in set [i_low, i_high)
Vector compute_cov(reb_simulation *const n_body_sim, int i_low,
		int i_high);
// Move center of mass of particles in set [i_low, i_high)
void subtract_com(reb_simulation *const n_body_sim, int i_low,
		int i_high);
// Move center of velocity of particles in set [i_low, i_high)
void subtract_cov(reb_simulation *const n_body_sim, int i_low,
		int i_high);
// Center simulation on center of mass of particles in set [i_low, i_high)
void center_sim(reb_simulation *const n_body_sim, int i_low, int i_high);

/***********************/
/* Rotational routines */
/***********************/

// Calculate moment of inertia tensor for particles in set [i_low, i_high)
Matrix mom_inertia(reb_simulation *const n_body_sim, int i_low,
		int i_high);
// Calculate (spin) angular momentum about center of mass of particles in set [i_low, i_high)
Vector measure_L(reb_simulation *const n_body_sim, int i_low,
		int i_high);
// Calculate spin angular momentum of particles in set [i_low, i_high) about center of mass using inverse moment of inertia
Vector body_spin(reb_simulation *const n_body_sim, int i_low, int i_high,
		double eigs[3]);
// Start particles in set [i_low, i_high) spinning around center of mass with spin vector omega
void spin_body(reb_simulation *const n_body_sim, int i_low, int i_high,
		Vector omega);
// Calculate orbital angular momentum
Vector compute_Lorb(reb_simulation *const n_body_sim, int i_low, int i_high);
// Calculate non-translational kinetic energy of particles in set [i_low, i_high)
double compute_rot_kin(reb_simulation *const n_body_sim, int i_low,
		int i_high);
// Compute (spin) angular momentum vector of particles in range [i_low, i_high) with respect to origin
Vector measure_L_origin(reb_simulation *const n_body_sim, int i_low,
		int i_high);
// Rotate particles in set [i_low, i_high) about their center of mass and velocity by Euler angles alpha, beta, gamma
void rotate_body(reb_simulation *const n_body_sim, int i_low, int i_high,
		double alpha, double beta, double gamma);
// Rotate particles in set [i_low, i_high) about the origin by Euler angles alpha, beta, gamma
void rotate_origin(reb_simulation *const n_body_sim, int i_low,
		int i_high, double alpha, double beta, double gamma);
// Rotate particles in set [i_low, i_high) about origin so that z is along largest eigenvalue, x is along smallest
void rotate_to_principal(reb_simulation *const n_body_sim, int i_low,
		int i_high);

/*******************/
/* Linear routines */
/*******************/

// Shift particles in range [i_low, i_high) by position dx, velocity dv
void move_resolved(reb_simulation *const n_body_sim, Vector dx, Vector dv,
		int i_low, int i_high);

/********************/
/* Energy and power */
/********************/

// Compute power lost to damping for a particular spring
double dEdt(reb_simulation *const n_body_sim, spring spr);
// Compute power lost to damping over all springs
double dEdt_total(reb_simulation *const n_body_sim);
// Calculate total potential energy in spring network
double spring_potential_energy(reb_simulation *const n_body_sim);
// Calculate total potential energy in particles in range [i_ow, i_high)
double grav_potential_energy(reb_simulation *const n_body_sim, int i_low,
		int i_high);

/********/
/* Misc */
/********/

// Zero all accelerations in Rebound n-body simulation
void zero_accel(reb_simulation *n_body_sim);
// Multiply masses to the right of x_min by m_fac and renormalize
void adjust_mass_side(reb_simulation *const n_body_sim, double m_fac,
		double x_min);
// Multiply masses either inside or outside of ellipse defined by center x0, semiaxes a,b,c by m_fac and renormalize
void adjust_mass_ellipsoid(reb_simulation *const n_body_sim, double m_fac, double a, double b, double c, Vector x0, bool inside);

#endif /* SRC_SPRING_PHYSICS_H_ */
