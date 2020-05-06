/*
 * heat.h
 *
 * Apply tidal heating to particles
 *
 *  Created on: Apr 15, 2020
 *      Author: alex
 */

#ifndef SRC_SPRING_HEAT_H_
#define SRC_SPRING_HEAT_H_

#include <vector>
using std::vector;

// Node structure for spring heating
struct node {
	bool is_surf;	// Is node a surface node?
	double temp;	// Temperature of node
	double cv;		// Specific heat of node
};

// Set up node vector
void init_nodes(reb_simulation *const n_body_sim, double Cv, double T0);
// Initialize temperature of surface nodes
void init_surf_temp(double T_surf);
// Calculate average thermal conductivity of springs with a radius from center of mass of particles in range [i_low, i_high) between r_min and r_max
double therm_cond_mesh(reb_simulation* const n_body_sim, int i_low, int i_high, double r_min, double r_max);
// Transport heat along springs
void transport_heat(reb_simulation *const n_body_sim, double dt);
// Apply spring heating to nodes
void rec_power(reb_simulation *const n_body_sim);
// Apply tidal heating to internal nodes
void heat_int_nodes_tidal(reb_simulation *const n_body_sim, double dt);

#endif /* SRC_SPRING_HEAT_H_ */
