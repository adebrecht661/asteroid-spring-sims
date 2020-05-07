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
using std::abs;

// Global values
int num_springs = 0;	// Number of springs
vector<spring> springs;	// Array of springs

// Default global values
int num_perts = 0;					// Number of perturbers
double def_gamma, t_damp, print_interval, pulse_period, pulse_wait, pulse_amp,
		chosen_print_interval, grav_accel, base_pos;	// Read in - see cfg
int chosen_particle; 				// printing out info on this particle
string fileroot;					// Root name for output files
std::ofstream pfile;				// Output file for chosen particle

// Global scales
double mass_scale, time_scale, length_scale, temp_scale, omega_scale, vel_scale,
		p_scale, L_scale, a_scale, F_scale, E_scale, dEdt_scale, P_scale;

// Forward declarations
void heartbeat(reb_simulation *const n_body_sim);
void sphere_forces(reb_simulation *const n_body_sim);
void pulse_base(reb_simulation *const n_body_sim);
void additional_forces(reb_simulation *n_body_sim);

int main(int argc, char *argv[]) {

	// Read in configuration file
	Config cfg;
	cfg.readFile("problem.cfg");

	// Get scales
	read_scales(&cfg);

	// Vars read in
	double r_cube, t_max, dt, max_spring_dist, gamma_fac, k;
	int num_parts;

	if (!(cfg.lookupValue("fileroot", fileroot) && cfg.lookupValue("dt", dt)
			&& cfg.lookupValue("t_max", t_max)
			&& cfg.lookupValue("r_cube", r_cube)
			&& cfg.lookupValue("print_interval", print_interval)
			&& cfg.lookupValue("max_spring_dist", max_spring_dist)
			&& cfg.lookupValue("k_short", k)
			&& cfg.lookupValue("gamma", def_gamma)
			&& cfg.lookupValue("damp_fac", gamma_fac)
			&& cfg.lookupValue("t_damp", t_damp)
			&& cfg.lookupValue("num_particles", num_parts)
			&& cfg.lookupValue("chosen_particle", chosen_particle)
			&& cfg.lookupValue("chosen_particle_interval",
					chosen_print_interval)
			&& cfg.lookupValue("pulse_amp", pulse_amp)
			&& cfg.lookupValue("grav_accel", grav_accel)
			&& cfg.lookupValue("pulse_period", pulse_period)
			&& cfg.lookupValue("pulse_wait", pulse_wait)
			&& cfg.lookupValue("base_pos", base_pos))) {
		throw "Failed to read in problem.cfg. Exiting.";
	} else {
		std::cout << "Read in problem.cfg." << std::endl;
	}

	// Set up rebound n-body simulation
	reb_simulation *const n_body_sim = reb_create_simulation();

	// Set up rebound constants
	n_body_sim->integrator = reb_simulation::REB_INTEGRATOR_LEAPFROG;
	n_body_sim->gravity = reb_simulation::REB_GRAVITY_NONE;
	n_body_sim->boundary = reb_simulation::REB_BOUNDARY_NONE;
	n_body_sim->G = 0;
	n_body_sim->additional_forces = additional_forces;

	// Generate particles
	uniform_line(n_body_sim, num_parts, base_pos);

	// Configure more rebound stuff
	n_body_sim->dt = dt;
	const double boxsize = 1.1 * r_cube;
	reb_configure_box(n_body_sim, boxsize, 1, 2, 1);

	// Set up spring
	spring def_spring;
	def_spring.gamma = def_gamma * gamma_fac;
	def_spring.k_heat = 0.0;
	def_spring.k = k;

	// Connect springs within L_overlap
	double L_overlap = 1.0 / num_parts;
	connect_springs_dist(n_body_sim, L_overlap, 0, n_body_sim->N, def_spring);

	// Extra particle info file
	string pfilename = fileroot + "_p" + std::to_string((int) num_parts / 2)
			+ ".txt";
	pfile.open(fileroot, std::ios::out | std::ios::app);

	// Function called every timestep
	n_body_sim->heartbeat = heartbeat;

	// Start integration
	if (t_max == 0.0) {
		reb_integrate(n_body_sim, INFINITY);
	} else {
		reb_integrate(n_body_sim, t_max);
	}
}

// Things done every timestep
void heartbeat(reb_simulation *const n_body_sim) {
	// Get particle infp
	reb_particle *particles = n_body_sim->particles;

	// Every 10 timesteps, output simulation timing info
	if (reb_output_check(n_body_sim, 10.0 * n_body_sim->dt)) {
		reb_output_timing(n_body_sim, 0);
	}

	// Turn off damping at t_damp
	if (abs(n_body_sim->t - t_damp) < 0.9 * n_body_sim->dt) {
		set_gamma(def_gamma);
	}

	// Output particle info every print_interval
	if (reb_output_check(n_body_sim, print_interval)) {
		write_particles(n_body_sim, fileroot,
				(int) n_body_sim->t / print_interval);
	}

	// Output extra info on one particle in the middle once damping ends
	if (n_body_sim->t - t_damp > 0) {
		if (reb_output_check(n_body_sim, chosen_print_interval)) {
			pfile << std::setprecision(6) << n_body_sim->t - t_damp << "\t"
					<< particles[chosen_particle].y << "\t"
					<< particles[chosen_particle].vy << "\t"
					<< particles[chosen_particle].ay << "\n";
		}
	}

	// Apply boundary condition - sinusoidal pulse on bottom particle
	pulse_base(n_body_sim);

}

// Apply gravitational force
void sphere_forces(reb_simulation *const n_body_sim) {

	// Get particle info
	reb_particle *particles = n_body_sim->particles;

	// Acceleration of gravity
	for (int i = 0; i < n_body_sim->N; i++) {
		particles[i].ay -= grav_accel;
	}

	// First particle is fixed
	particles[0].ay = 0.0;
}

// Pulse base particle
void pulse_base(reb_simulation *const n_body_sim) {

	// Get particle info
	reb_particle *particles = n_body_sim->particles;

	// Start pulses at t_damp
	double t = n_body_sim->t - t_damp;
	double freq = 2.0 * M_PI / pulse_period;

	// Continue pulse if required
	t = fmod(t, (pulse_period + pulse_wait));
	if (t < pulse_period) {

		// Goes from 0 to pulse_amp back to 0
		double pos = pulse_amp * (1.0 - cos(t * freq)) / 2.0;
		// Derivative = velocity
		double pos_dot = pulse_amp * freq * sin(t * freq) / 2.0;

		// Set particle
		particles[0].y = base_pos + pos;
		particles[0].vy = pos_dot;
		particles[0].ay = 0.0;
	}
}

// Initialize accelerations (because not using rebound gravity) and apply forces
void additional_forces(reb_simulation *n_body_sim) {
	// Init accelerations
	zero_accel(n_body_sim);

	// Spring forces
	spring_forces(n_body_sim);

	// Gravity forces
	sphere_forces(n_body_sim);
}
