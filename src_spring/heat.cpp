#ifdef __cplusplus
# 	ifdef __GNUC__
#		define restrict __restrict__
#	else
#		define restrict
#	endif
#endif

#include <cmath>
#include <vector>
extern "C" {
#include "rebound.h"
}
#include "physics.h"
#include "springs.h"
#include "heat.h"

using std::vector;

extern vector<spring> springs; // Spring array
extern int num_springs; // Number of springs
extern int num_perts; // Number of perturbers
vector<node> nodes; // Node array

// Set up node vector
void init_nodes(struct reb_simulation *const n_body_sim, double Cv, double T0) {
	// Get default values
	node def_node;
	def_node.cv = Cv;
	def_node.is_surf = false;
	def_node.temp = T0;

	// Init n copies of default node
	nodes = vector<node>(n_body_sim->N, def_node);
}

// Conductive heat transport over spring network
// Update temperature on each node for the timestep with size dt
void transport_heat(reb_simulation* const n_body_sim, double dt) {
	// Allocate and initialize array for new temperatures
	double* delta_T = (double*) malloc(sizeof(double) * n_body_sim->N);
	for (int k = 0; k < n_body_sim->N; k++)
		delta_T[k] = 0.0;

	// Get energy transfer along each spring
	for (int i = 0; i < num_springs; i++) {
		int ii = springs[i].particle_1;
		int jj = springs[i].particle_2;

		// Heat flux
		double Q_ij = (nodes[jj].temp - nodes[ii].temp) * springs[i].k_heat;

		// Update temperature changes
		delta_T[ii] += Q_ij * dt / (n_body_sim->particles[ii].m * nodes[ii].cv);
		delta_T[jj] -= Q_ij * dt / (n_body_sim->particles[jj].m * nodes[jj].cv);
	}

	// Update node temperatures after all fluxes have been added up
	for (int k = 0; k < n_body_sim->N - num_perts; k++) {
		if (!nodes[k].is_surf) { // only do this to interior nodes
			nodes[k].temp += delta_T[k];
		}
	}

	// Free malloc'd array
	free(delta_T);
}

// Apply tidal heating to internal nodes
// Should only raise temperature
void heat_nodes_tidal(reb_simulation* const n_body_sim, double dt) {
	for (int k = 0; k < num_springs; k++) {

		// Get heat from each spring
		double power = dEdt(n_body_sim, springs[k]);

		// Get change in temperature of each node
		int ii = springs[k].particle_1;
		int jj = springs[k].particle_2;
		double dT_ii = 0.5 * power * dt / (n_body_sim->particles[ii].m * nodes[ii].cv);
		double dT_jj = 0.5 * power * dt / (n_body_sim->particles[jj].m * nodes[jj].cv);

		// Update temp of internal nodes
		if (!nodes[ii].is_surf)	nodes[ii].temp += dT_ii;
		if (!nodes[jj].is_surf) nodes[jj].temp += dT_jj;
	}
}

// Add a constant heating rate to each node (radiogenic)
// dot_E_rad is energy per unit mass per unit time
void heat_nodes_radiogenic(reb_simulation* n_body_sim,
		double dot_E_rad) {

	// Node for each non-perturbing particle
	for (int i = 0; i < n_body_sim->N - num_perts; i++) {

		// Get temperature change
		double dT_i = dot_E_rad * n_body_sim->dt / nodes[i].cv;

		// Update interior nodes
		// Radiogenic heating should affect all nodes??????
		if (!nodes[i].is_surf) {
			nodes[i].temp += dT_i;
		}
	}
}
