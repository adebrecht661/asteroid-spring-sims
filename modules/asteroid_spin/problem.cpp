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
#include <libconfig.h++>
extern "C" {
#include "rebound.h"
}
#include "matrix_math.h"
#include "input_spring.h"
#include "output_spring.h"
#include "springs.h"
#include "physics.h"
#include "stress.h"

using namespace libconfig;

using std::string;
using std::vector;

// Global values
int num_springs = 0;	// Global numbers of springs
vector<spring> springs;	// Global spring array

// Default global values
int num_perts = 0;		// Number of perturbers
double def_gamma, t_damp, print_interval;	// Read in - see cfg
string fileroot;   		// Output file base name

// Global scales
double mass_scale, time_scale, length_scale, temp_scale, omega_scale, vel_scale,
		p_scale, L_scale, a_scale, F_scale, E_scale, dEdt_scale, P_scale;

// Forward declarations
void additional_forces(reb_simulation *n_body_sim);
void heartbeat(reb_simulation *const n_body_sim);
void reb_springs(reb_simulation *const n_body_sim);

int main(int argc, char *argv[]) {

	// Read in configuration file
	Config cfg;
	cfg.readFile("problem.cfg");

	// Get scales
	read_scales(cfg);

	// Vars read in
	double r_ball, t_max, dt, max_spring_dist, gamma_fac, k, omega_in[3];

	cfg.lookupValue("fileroot", fileroot);
	cfg.lookupValue("dt", dt);
	cfg.lookupValue("t_max", t_max);
	cfg.lookupValue("r_ball", r_ball);
	cfg.lookupValue("print_interval", print_interval);
	cfg.lookupValue("max_spring_dist", max_spring_dist);
	cfg.lookupValue("k_short", k);
	cfg.lookupValue("gamma", def_gamma);
	cfg.lookupValue("damp_fac", gamma_fac);
	cfg.lookupValue("t_damp", t_damp);
	cfg.lookupValue("omega_x", omega_in[0]);
	cfg.lookupValue("omega_y", omega_in[1]);
	cfg.lookupValue("omega_z", omega_in[2]);
	Vector omega(omega_in);

	std::cout << "Read in problem.cfg." << std::endl;

	// Create rebound simulation
	reb_simulation *const n_body_sim = reb_create_simulation();

	// Set rebound constants
	n_body_sim->integrator = reb_simulation::REB_INTEGRATOR_LEAPFROG;
	n_body_sim->gravity = reb_simulation::REB_GRAVITY_NONE;
	n_body_sim->boundary = reb_simulation::REB_BOUNDARY_NONE;
	n_body_sim->G = 6.67428e-8 / F_scale * pow(length_scale, 2.0)
			/ pow(mass_scale, 2.0);
	n_body_sim->additional_forces = additional_forces;

	// Set more rebound parameters
	n_body_sim->dt = dt;								// Integration timestep
	const double boxsize = 3.2 * r_ball;				// Size of display box
	reb_configure_box(n_body_sim, boxsize, 1, 1, 1);// Configure rebound box (last three arguments are number of root boxes in x,y,z directions)
	n_body_sim->softening = 1e-6;			// Gravitational softening length

	// Set up default spring parameters
	spring default_spring;
	default_spring.gamma = gamma_fac * def_gamma;// initial damping coefficient
	default_spring.k = k;							// spring constant

	// Open output file - overwrites any existing files
	std::ofstream outfile(fileroot + "_run.txt",
			std::ios::out | std::ios::trunc);

	// Read in particles here
	read_particles(n_body_sim, fileroot, 0);

	// Get index range for resolved body
	int i_low = 0;
	int i_high = n_body_sim->N;

	// Move reference frame to resolved body
	subtract_com(n_body_sim, i_low, i_high);

	// Make relative velocity zero
	subtract_cov(n_body_sim, i_low, i_high);

	// Add spin to resolved body
	spin_body(n_body_sim, i_low, i_high, omega);

	// Should not have added any velocity??????
	//subtract_cov(n_body_sim, i_low, i_high);

	// Print spin period to file
	double spin_period = abs(2.0 * M_PI / omega.len());
	std::cout << "Spin period: " << std::setprecision(6) << spin_period << "\n";
	outfile << "spin period " << std::setprecision(6) << spin_period << "\n";

	// Make springs
	// Connect all particles within max_spring_dist
	connect_springs_dist(n_body_sim, max_spring_dist, i_low, i_high,
			default_spring);

	// Print Young's modulus
	double r = 0.4; 					// Radius for computing Young modulus
	double E_mesh = Young_mesh(n_body_sim, i_low, i_high, 0.0, r);
	double E_mesh_all = Young_full_mesh();
	std::cout << std::setprecision(3) << "r = " << r << " mesh_distance = "
			<< max_spring_dist << "\n";
	std::cout << std::setprecision(6) << "Youngs_modulus " << E_mesh
			<< "\nYoungs_modulus_all " << E_mesh_all << "\n";
	outfile << std::setprecision(6) << "Youngs_modulus " << E_mesh
			<< "\nYoungs_modulus_big " << E_mesh_all << "\n";

	// Print max spring distance and mean spring length
	outfile << std::setprecision(4) << "mesh_distance " << max_spring_dist
			<< "\n";
	double L = mean_spring_length();
	std::cout << std::setprecision(4) << "mean_L " << L << "\n";
	outfile << std::setprecision(4) << "mean_L " << L << "\n";

	// Ratio of numbers of particles to numbers of springs for resolved body
	double Nratio = (double) num_springs / (i_high - i_low);
	std::cout << "N = " << n_body_sim->N << ", NS = " << num_springs
			<< ", NS/N = " << Nratio << std::endl;
	outfile << "N " << n_body_sim->N << "\n";
	outfile << "NS " << num_springs << "\n";
	outfile << std::setprecision(1) << "NS/N " << Nratio << "\n";
	outfile.close();

	// Pass springs to rebound for display
	reb_springs(n_body_sim);

	// Set up calls for every integration timestep
	n_body_sim->heartbeat = heartbeat;

	// Move reference frame to resolved body
	center_sim(n_body_sim, i_low, i_high);

	// Integrate rebound simulation
	if (t_max == 0.0) {
		reb_integrate(n_body_sim, INFINITY);
	} else {
		reb_integrate(n_body_sim, t_max);
	}
}

// Do this every rebound timestep
void heartbeat(reb_simulation *const n_body_sim) {
	// Get filename once
	static string extendedfile = fileroot + "_ext.txt";

	// Every 10 timesteps, output rebound timing info
	if (reb_output_check(n_body_sim, 10.0 * n_body_sim->dt)) {
		reb_output_timing(n_body_sim, 0);
	}

	// Damp initial bounce only
	// Reset gamma at t near t_damp
	if (abs(n_body_sim->t - t_damp) < 0.9 * n_body_sim->dt)
		set_gamma(def_gamma);

	// Move reference frame to resolved body for display
	center_sim(n_body_sim, 0, n_body_sim->N - num_perts);

	// Check if we've hit a print time
	if (reb_output_check(n_body_sim, print_interval)) {
		// If so, set file index and write out spring and particle info
		int index = (int) (n_body_sim->t / print_interval);
		if (index > 0) {
			write_springs(n_body_sim, fileroot, index);
			write_particles(n_body_sim, fileroot, index);
		}
	}
}

// Make a spring index list for display
void reb_springs(reb_simulation *const n_body_sim) {
	// Set number of springs and initialize arrays for spring-particle interface
	n_body_sim->NS = num_springs;
	n_body_sim->springs_i = (int*) malloc(num_springs * sizeof(int));
	n_body_sim->springs_j = (int*) malloc(num_springs * sizeof(int));

	// For each spring, note particle indices
	for (int i = 0; i < num_springs; i++) {
		n_body_sim->springs_i[i] = springs[i].particle_1;
		n_body_sim->springs_j[i] = springs[i].particle_2;
	}
}

// Passed to rebound so you can have forces other than gravity
void additional_forces(reb_simulation *n_body_sim) {
	// If you're not using gravity routines in rebound (i.e. REB_GRAVITY_NONE), need to initialize particle accelerations
	zero_accel(n_body_sim);

	// Apply spring forces
	spring_forces(n_body_sim);
}
