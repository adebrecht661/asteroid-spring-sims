/*
 * orb.h
 *
 *  Created on: Apr 15, 2020
 *      Author: alex
 */

#ifndef SRC_SPRING_ORB_H_
#define SRC_SPRING_ORB_H_

#include "kepcart.h"

/************************/
/* Perturber operations */
/************************/

// Add binary perturbers with primary mass m_prim, secondary mass m_prim * m_ratio, display radius r_prim
// Resolved body orbits center of mass of binary system in orbit defined by orb_el
double add_bin_kep(reb_simulation *const n_body_sim, double m_prim,
		double r_prim, double m_ratio, double sep, OrbitalElements orb_el);
// Add perturbing mass of mass and radius
// if i_p < 0, mass is placed at origin and resolved body orbits it in orbit defined by orb_el
// if i_p >= 0, mass is placed into orbit defined by orb_el around mass at i_p
double add_pt_mass_kep(reb_simulation *const n_body_sim, int i_low, int i_high,
		int i_p, double mass, double radius, OrbitalElements orb_el);
// Drift orbits of binary particles 1 and 2
void drift_bin(reb_simulation *const n_body_sim, double timestep,
		double inv_tau_a, double inv_tau_e, int part_1, int part_2);
// Drift orbits of particle and resolved body
void drift_resolved(reb_simulation *const n_body_sim, double timestep,
		double inv_tau_a, double inv_tau_e, int i_part, int i_low, int i_high);
// Apply quadrupole moment acceleration from particle i_p (with dimensionless quadrupole coefficient J2_p)
void quadrupole_accel(reb_simulation *const n_body_sim, double J2_p, double R_p,
		double phi_p, double theta_p, int i_p);

/***************************************/
/* Orbital properties of resolved body */
/***************************************/

// Compute orbital properties of resolved body with respect to center of mass of perturbing masses
void compute_orb(reb_simulation *const n_body_sim, int i_low, int i_high,
		double *a, double *mean_motion, double *e, double *i, double *L);

/********/
/* Misc */
/********/

// Get total mass of particles in range [i_low, i_high)
double sum_mass(reb_simulation *const n_body_sim, int i_low, int i_high);

#endif /* SRC_SPRING_ORB_H_ */
