/*
 * heat.h
 *
 * Apply tidal or radiogenic heating to particles
 *
 *  Created on: Apr 15, 2020
 *      Author: alex
 */

#ifndef SRC_SPRING_HEAT_H_
#define SRC_SPRING_HEAT_H_

// Node structure for spring heating
typedef struct node {
	bool is_surf;	// Is node a surface node?
	double temp;	// Temperature of node
	double cv;		// Specific heat of node (integrated over mass of node?????)
} node;

// Transport heat along springs
void transport_heat(struct reb_simulation *const n_body_sim, node nodes[],
		double dt);
// Apply tidal heating to internal nodes
void heat_nodes_tidal(struct reb_simulation *const n_body_sim, node nodes[],
		double dt);
// Apply radiogenic heating to internal nodes
void heat_nodes_radiogenic(struct reb_simulation *n_body_sim, node nodes[],
		double dot_E_rad);

#endif /* SRC_SPRING_HEAT_H_ */
