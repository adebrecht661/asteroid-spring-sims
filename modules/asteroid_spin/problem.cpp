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
#include "springs.h"
#include "physics.h"
#include "stress.h"

using std::string;
using std::vector;

int num_springs;	// Global numbers of springs
vector<spring> springs;	// Global spring array

double def_gamma;		// Default damping coefficient for all springs
double t_damp;			// Time to end stronger damping coefficient
double print_interval;	// Interval to print out spring and particle info
string fileroot;   		// Output file base name
int num_perts = 0;		// Number of perturbers

// Forward declarations
void additional_forces(reb_simulation *n_body_sim);
void heartbeat(reb_simulation *const n_body_sim);
void reb_springs(reb_simulation *const n_body_sim);

int main(int argc, char *argv[]) {
	// Create rebound simulation
	reb_simulation *const n_body_sim = reb_create_simulation();

	// Set rebound constants
	n_body_sim->integrator = reb_simulation::REB_INTEGRATOR_LEAPFROG;
	n_body_sim->gravity = reb_simulation::REB_GRAVITY_NONE;
	n_body_sim->boundary = reb_simulation::REB_BOUNDARY_NONE;
	n_body_sim->G = 1;
	n_body_sim->additional_forces = additional_forces;

	// Max integration time
	// (0 = integrate forever)
	double tmax = 0.0;

	// Values to set
	double dt = 0.000;
	Vector omega = 0.0;
	double k = 0.0;
	double gamma_fac = 0;
	double max_spring_dist = 0.0;

	// 1 argument, no filename passed
	if (argc == 1) {
		fileroot = "a2";			// Output file basename
		dt = 1e-2;					// Integration timestep
		tmax = 0.0;					// Max integration time
		print_interval = 100.0;     // Interval at which to print results
		k = 0.005;					// Spring constant
		def_gamma = 1.0;			// Base spring damping coefficient
		gamma_fac = 5.0;// Factor by which initial damping is higher that gamma
		t_damp = 1.0;				// Turn of damping at this time
		omega = { 0, 0, 0.8 };		// Initial spin
		max_spring_dist = 0.1;// Max distance to connect particles with springs
		// More than one argument, filename to read from
	} else {
		throw "Error: Reading in from file not currently implemented.";
	}

	// Set more rebound parameters
	n_body_sim->dt = dt;								// Integration timestep
	double rball = 1.0;								// Radius of resolved body
	const double boxsize = 3.2 * rball;					// Size of display box
	reb_configure_box(n_body_sim, boxsize, 1, 1, 1);// Configure rebound box (last three arguments are number of root boxes in x,y,z directions)
	n_body_sim->softening = 1e-6;			// Gravitational softening length

	// Set up default spring parameters
	spring default_spring;
	default_spring.gamma = gamma_fac * def_gamma; // initial damping coefficient
	default_spring.k = k; // spring constant

	// Open output file - overwrites any existing files
	std::ofstream outfile(fileroot + "_run.txt",
			std::ios::out | std::ios::trunc);

	// Start with no springs
	num_springs = 0;

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
	double speriod = abs(2.0 * M_PI / omega.len());
	std::cout << "Spin period: " << std::setprecision(6) << speriod << "\n";
	outfile << "spin period " << std::setprecision(6) << speriod << "\n";

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
	if (tmax == 0.0) {
		reb_integrate(n_body_sim, INFINITY);
	} else {
		reb_integrate(n_body_sim, tmax);
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
