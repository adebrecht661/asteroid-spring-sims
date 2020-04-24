#ifdef __cplusplus
# 	ifdef __GNUC__
#		define restrict __restrict__
#	else
#		define restrict
#	endif
#endif

/**
 * resolved mass spring model
 * using the leap frog integrator.
 */
#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
extern "C" {
#include "rebound.h"
}
#include "matrix_math.h"
#include "input_spring.h"
#include "output_spring.h"
#include "shapes.h"
#include "springs.h"
#include "physics.h"
#include "stress.h"

using std::string;
using std::vector;

// Global values
int num_springs = 0;				// Number of springs
vector<spring> springs;				// Array of springs

// Default global values
int num_perts = 0;					// Number of perturbers
double def_gamma = 1.0;				// Default spring damping coefficient
double t_damp = 1.0;				// Time to end increased damping
double print_interval = 100000.0;	// Interval to print output
string fileroot = "t1";				// Root name for output files
double surf_dist = 0.1;	// Distance from surface to count as a surface particle
double pressure = 0.0;				// Pressure from above

// Forward declarations
void reb_springs(reb_simulation *const n_body_sim);
void table_top(reb_simulation *const n_body_sim);
void top_push(reb_simulation *const n_body_sim);
void heartbeat(reb_simulation *const n_body_sim);
void additional_forces(reb_simulation *n_body_sim);

int main(int argc, char *argv[]) {
	// Set up rebound n-body simulation
	struct reb_simulation *const n_body_sim = reb_create_simulation();

	// Set up rebound constants
	n_body_sim->integrator = reb_simulation::REB_INTEGRATOR_LEAPFROG;
	n_body_sim->gravity = reb_simulation::REB_GRAVITY_NONE;
	n_body_sim->boundary = reb_simulation::REB_BOUNDARY_NONE;
	n_body_sim->G = 0;
	n_body_sim->additional_forces = additional_forces;

	// Set defaults of other values
	double cube_mass = 1.0;				// Total mass of all particles
	double cube_side_len = 1.0;			// Side length of cube
	double t_max = 0.0;		// Max integration time - if 0, integrate forever
	double dt = 1e-3;					// Integration timestep
	double min_part_dist = 0.15;		// Minimum interparticle spacing
	double max_short_spring_fac = 2.3;// fac*min_part_distance is maximum spring length
	double max_med_spring_fac = 2.3;// fac*min_part_distance is maximum spring length
	double max_long_spring_fac = 2.3;// fac*min_part_distance is maximum spring length
	double gamma_fac = 0.0;	// Factor by which to increase damping until t_damp
	double k_short = 8e-2;				// Spring constant for short springs
	double k_med = 8e-2;				// Spring constant for medium springs
	double k_long = 8e-2;				// Spring constant for long springs

	// Read in a parameter file, if one is provided
	// Caution: no real error checking
	if (argc > 1) {
		std::ifstream param_file(argv[1], std::ios::in);
		param_file >> fileroot;
		param_file >> dt;
		param_file >> t_max;
		param_file >> print_interval;
		param_file >> min_part_dist;
		param_file >> max_short_spring_fac;
		param_file >> max_med_spring_fac;
		param_file >> max_long_spring_fac;
		param_file >> k_short;
		param_file >> k_med;
		param_file >> k_long;
		param_file >> def_gamma;
		param_file >> gamma_fac;
		param_file >> t_damp;
		param_file >> surf_dist;
		param_file >> pressure;

		std::cout << "Read in " << argv[1] << std::endl;
	}

	// Set simulation parameters based on input
	n_body_sim->dt = dt;								// Integration timestep
	const double box_size = 1.1 * cube_side_len;		// Display window size
	reb_configure_box(n_body_sim, box_size, 1, 1, 1);// Last three terms are number of root boxes in x,y,z directions
	n_body_sim->softening = min_part_dist / 100.0;// Gravitational softening length

	// Set default spring parameters
	spring def_spring;
	def_spring.gamma = def_gamma;
	def_spring.k = k_short;
	def_spring.k_heat = 1.0;
	double short_spring_dist = min_part_dist * max_short_spring_fac;
	double med_spring_dist = min_part_dist * max_med_spring_fac;
	double long_spring_dist = min_part_dist * max_long_spring_fac;

	// Filename for info output
	string filename;
	filename = fileroot + "_run.txt";

	// Create rectangular particle distribution
	rand_rectangle(n_body_sim, min_part_dist, cube_side_len, cube_side_len,
			cube_side_len, cube_mass);

	// Connect particles within interparticle distance by springs
	// Short springs
	connect_springs_dist(n_body_sim, short_spring_dist, 0, n_body_sim->N,
			def_spring);
	// Medium springs
	def_spring.k = k_med;
	connect_springs_dist(n_body_sim, med_spring_dist, 0, n_body_sim->N,
			def_spring);
	// Long springs
	def_spring.k = k_long;
	connect_springs_dist(n_body_sim, long_spring_dist, 0, n_body_sim->N,
			def_spring);

	// Pass springs to display
	reb_springs(n_body_sim);

	// Start with enhanced damping
	set_gamma(def_gamma * gamma_fac);

	// Set function called at every integration timestep
	n_body_sim->heartbeat = heartbeat;

	// Integrate simulation
	if (t_max == 0.0) {
		reb_integrate(n_body_sim, INFINITY);
	} else {
		reb_integrate(n_body_sim, t_max);
	}
}

// Things to happen every integration timestep
void heartbeat(struct reb_simulation *const n_body_sim) {

	// Output rebound timing info every 10 timesteps
	if (reb_output_check(n_body_sim, 10.0 * n_body_sim->dt)) {
		reb_output_timing(n_body_sim, 0);
	}

	// Turn off damping at t_damp
	if (abs(n_body_sim->t - t_damp) < 0.9 * n_body_sim->dt) {
		set_gamma(def_gamma);
	}

	// Output particle info every print interval
	if (reb_output_check(n_body_sim, print_interval)) {
		write_particles(n_body_sim, fileroot,
				(int) n_body_sim->dt / print_interval);
	}
}

// Initialize accelerations (because not using rebound gravity) and add forces other than gravity
void additional_forces(reb_simulation *n_body_sim) {
	// Init accelerations
	zero_accel(n_body_sim);

	// Springs
	spring_forces(n_body_sim);

	// Pressure from top
	top_push(n_body_sim);

	// Sitting on a table
	table_top(n_body_sim);
}

// Make a spring index list to pass to viewer
void reb_springs(struct reb_simulation *const n_body_sim) {
	// Set rebound spring info
	n_body_sim->NS = num_springs;
	n_body_sim->springs_i = (int*) malloc(num_springs * sizeof(int));
	n_body_sim->springs_j = (int*) malloc(num_springs * sizeof(int));

	// Add each spring to arrays
	for (int i = 0; i < num_springs; i++) {
		n_body_sim->springs_i[i] = springs[i].particle_1;
		n_body_sim->springs_j[i] = springs[i].particle_2;
	}
}

// Apply pressure to top of cube
// Mark top surface in the first call, since it can move
void top_push(struct reb_simulation *const n_body_sim) {
	// Get particle info
	reb_particle *particles = n_body_sim->particles;

	// List and number of top surface particles
	static vector<int> surf_list(n_body_sim->N, -1);
	static int N_surf = 0;

	// Mark surface particles at first run
	static bool first = true;
	if (first) {
		first = false;
		for (int i = 0; i < n_body_sim->N; i++) {
			if (particles[i].y > 0.5 - surf_dist) {
				surf_list[N_surf] = i;
				N_surf++;
			}
		}
		std::cout << "There are " << N_surf
				<< " particles on the top of the cube." << std::endl;
	}

	// Apply pressure to surface particles
	for (int j = 0; j < N_surf; j++) {
		int i = surf_list[j];
		particles[i].ay -= pressure / N_surf / particles[i].m;
	}
}

// Cube can't move through table
void table_top(struct reb_simulation *const n_body_sim) {
	// Get particle info
	reb_particle *particles = n_body_sim->particles;

	// List and number of top surface particles
	static vector<int> surf_list(n_body_sim->N, -1);
	static int N_surf = 0;

	// Mark surface particles at first run
	static bool first = true;
	if (first) {
		first = false;
		for (int i = 0; i < n_body_sim->N; i++) {
			if (particles[i].y < -0.5 + surf_dist) {
				surf_list[N_surf] = i;
				N_surf++;
			}
		}
		std::cout << "There are " << N_surf
				<< " particles on the bottom of the cube." << std::endl;
	}

	// Stop particles from accelerating
	for (int j = 0; j < N_surf; j++) {
		int i = surf_list[j];
		particles[i].ay = 0.0;
	}
}
