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
#include "libconfig.h++"
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
#include "orb.h"

using namespace libconfig;

using std::string;
using std::vector;
using std::abs;

// Global values
int num_springs = 0;	// Global numbers of springs
vector<spring> springs;	// Global spring array

// Default global values
int num_perts = 0;		// Number of perturbers
double def_gamma, t_damp, print_interval, J2_plus, R_plus, theta_plus, phi_plus;// Read in - see cfg
string fileroot;		// Output file base name

// Central mass location
// Should be set to -1 initially
int i_central = -1;

// Lattice selectors
const int RAND_ELLIPSE = 0, HCP_ELLIPSE = 1, CUBIC_ELLIPSE = 2, BENNU = 4;

// Global scales
double mass_scale, time_scale, length_scale, temp_scale, omega_scale, vel_scale,
		p_scale, L_scale, a_scale, F_scale, E_scale, dEdt_scale, P_scale;

// Forward declarations
void heartbeat(reb_simulation *const n_body_sim);
void reb_springs(reb_simulation *const r);
void additional_forces(reb_simulation *n_body_sim);

int main(int argc, char *argv[]) {

	// Read in configuration file
	Config cfg;
	cfg.readFile("problem.cfg");

	// Get scales
	read_scales(&cfg);

	// Vars read in
	double t_max, dt, max_spring_dist, gamma_fac, k, omega_in[3], min_part_dist,
			ratio1, ratio2, m_ball, r_ball, R_shell, k_int, gamma_int;
	int surface_type;

	if (!(cfg.lookupValue("fileroot", fileroot) && cfg.lookupValue("dt", dt)
			&& cfg.lookupValue("t_max", t_max)
			&& cfg.lookupValue("print_interval", print_interval)
			&& cfg.lookupValue("max_spring_dist", max_spring_dist)
			&& cfg.lookupValue("k", k) && cfg.lookupValue("gamma", def_gamma)
			&& cfg.lookupValue("damp_fac", gamma_fac)
			&& cfg.lookupValue("t_damp", t_damp)
			&& cfg.lookupValue("mass_ball", m_ball)
			&& cfg.lookupValue("radius_ball", r_ball)
			&& cfg.lookupValue("omega_x", omega_in[0])
			&& cfg.lookupValue("omega_y", omega_in[1])
			&& cfg.lookupValue("omega_z", omega_in[2])
			&& cfg.lookupValue("surface_type", surface_type)
			&& cfg.lookupValue("min_particle_distance", min_part_dist)
			&& cfg.lookupValue("ratio1", ratio1)
			&& cfg.lookupValue("ratio2", ratio2)
			&& cfg.lookupValue("R_shell", R_shell)
			&& cfg.lookupValue("k_int", k_int)
			&& cfg.lookupValue("gamma_int", gamma_int))) {
		throw "Failed to read in problem.cfg. Exiting.";
	} else {
		std::cout << "Read in problem.cfg." << std::endl;
	}
	Vector omega(omega_in);

	// Create rebound simulation
	reb_simulation *const n_body_sim = reb_create_simulation();

	// Set rebound constants
	n_body_sim->integrator = reb_simulation::REB_INTEGRATOR_LEAPFROG;
	n_body_sim->gravity = reb_simulation::REB_GRAVITY_NONE;
	n_body_sim->boundary = reb_simulation::REB_BOUNDARY_NONE;
	n_body_sim->G = 6.67428e-8 / F_scale / pow(length_scale, 2.0)
			* pow(mass_scale, 2.0);
	n_body_sim->additional_forces = additional_forces;

	// Set more rebound parameters
	n_body_sim->dt = dt;							// Integration timestep
	n_body_sim->softening = min_part_dist / 100.0;// Gravitational softening length

	// Set display window size
	const double boxsize = 3.2 * r_ball;
	reb_configure_box(n_body_sim, boxsize, 1, 1, 1);

	// Spring properties
	spring def_spring;
	def_spring.gamma = gamma_fac * def_gamma;
	def_spring.k = k;

	// Get output filename
	string filename = fileroot + "_run.txt";
	std::ofstream outfile(filename, std::ios::out | std::ios::app);

	// Adjust volume to be that of sphere with given radius
	double volume_ratio = pow(r_ball, 3.0) * ratio1 * ratio2;
	double vol_radius = pow(volume_ratio, 1.0 / 3.0);
	r_ball /= vol_radius;

	// Output semi-axes
	outfile << std::setprecision(3) << "a " << r_ball << "\n";
	outfile << std::setprecision(3) << "b " << r_ball * ratio1 << "\n";
	outfile << std::setprecision(3) << "c " << r_ball * ratio2 << "\n";

	// Should be 1 now
	volume_ratio = pow(r_ball, 3.0) * ratio1 * ratio2;
	outfile << std::setprecision(6) << "vol_ratio " << volume_ratio << "\n";

	// Create resolved body particle distribution
	switch (surface_type) {
	case RAND_ELLIPSE: {
		rand_ellipsoid(n_body_sim, min_part_dist, r_ball, r_ball * ratio1,
				r_ball * ratio2, m_ball);
		break;
	}
	case HCP_ELLIPSE: {
		hcp_ellipsoid(n_body_sim, min_part_dist, r_ball, r_ball * ratio1,
				r_ball * ratio2, m_ball);
		break;
	}
	case CUBIC_ELLIPSE: {
		cubic_ellipsoid(n_body_sim, min_part_dist, r_ball, r_ball * ratio1,
				r_ball * ratio2, m_ball);
		break;
	}
	case BENNU: {

		// Read in vertex file
		read_vertices(n_body_sim, "101955bennu.tab");

		// Note final shape vertex
		int N_bennu = n_body_sim->N;
		std::cout << "Bennu shape model read in.\n";

		// Correct units from km to length_scale
		stretch(n_body_sim, 0, n_body_sim->N, 1.0 / length_scale);

		// Fill shape with particles
		rand_shape(n_body_sim, min_part_dist, 1.0);
		std::cout << "Shape filled.\n";

		// Clean up shape vertices
		rm_particles(n_body_sim, 0, N_bennu);
		std::cout << "Shape vertices deleted.\n";
		break;
	}
	default: {
		throw "No shape of requested type available. Exiting.";
	}
	}

	// Select indices of resolved body
	int i_low = 0;
	int i_high = n_body_sim->N;

	// Move reference frame to resolved body and rotate
	subtract_com(n_body_sim, i_low, i_high);
	subtract_cov(n_body_sim, i_low, i_high);

	// Spin the body
	spin_body(n_body_sim, i_low, i_high, omega);

	// Connect all particles within max_spring_dist by hot springs
	connect_springs_dist(n_body_sim, max_spring_dist, 0, n_body_sim->N,
			def_spring);

	// Output Young's modulus of mesh
	double ddr = 0.4;
	print_run_double(max_spring_dist, "max spring length", &outfile);
	print_run_double(ddr, "ddr", &outfile);
	print_run_double(Young_mesh(n_body_sim, i_low, i_high, 0.0, ddr),
			"Young's modulus hot", &outfile);
	double E_mesh = Young_full_mesh();
	print_run_double(Young_full_mesh(), "Young's modulus big hot", &outfile);
	print_run_double(mean_spring_length(), "Mean spring length", &outfile);

	// Calculate Kelvin-Voigt relaxation time and viscosity
	// Factor of 0.5 is due to reduced mass being used in calculation
	double tau_relax = 1.0 * def_gamma * 0.5
			* (m_ball / (n_body_sim->N - num_perts)) / def_spring.k;
	print_run_double(tau_relax, "relaxation time", &outfile);
	print_run_double(tau_relax * E_mesh / 2.5, "viscosity", &outfile);

	// Change properties of interior of shell
	std::cout << "Rshell = " << R_shell << "\n";
	adjust_spring_props(n_body_sim, k_int, gamma_int, 0.0, R_shell);

	// Get eigenvalues of moment of inertia
	Matrix inertia_mat = mom_inertia(n_body_sim, 0, n_body_sim->N);
	double eigs[3];
	eigenvalues(inertia_mat, eigs);
	print_run_double(sqrt(eigs[2] / eigs[0]), " should be an axis ratio \n",
			&outfile);
	print_run_double(sqrt(eigs[2] / eigs[1]), " should be an axis ratio \n",
			&outfile);

	// Get current angular momentum
	Vector L = measure_L(n_body_sim, i_low, i_high);

	// Compute angle between angular momentum and principal axis and precession rate
	double theta = 0.0;
	double omega_3 = 0.0;
	double omega_prec = 0.0;
	// For oblate shape
	if (ratio1 == 1) {
		double hsquare = pow(ratio2, 2.0); // is (c/a)^2
		if (surface_type == BENNU) {
			hsquare = eigs[2] / eigs[0];
		}
		theta = acos(L.getZ() / L.len()); // axis of symmetry is z axis
		omega_3 = L.len() / eigs[0]; // eig1 = I3 is largest moment of inertia, corresponding to smallest body axis
		omega_prec = (1.0 - hsquare) / (1.0 + hsquare);
		omega_prec *= omega_3 * cos(theta);
		std::cout << "Oblate\n";
		// For prolate shape
	} else if (ratio1 == ratio2) {
		double hsquare = 1.0 / (ratio2 * ratio2); // is a/c
		theta = acos(L.getX() / L.len()); // axis of symmetry is x axis
		omega_3 = L.len() / eigs[2]; // eig3 is smallest moment,  is I parallel, largest body axis
		omega_prec = (1.0 - hsquare) / (1.0 + hsquare);
		omega_prec *= omega_3 * cos(theta);
		std::cout << "Prolate\n";
	}

	std::cout << "theta (deg)= " << theta * 180.0 / M_PI << " (radians)= "
			<< theta << "\n";
	outfile << "theta (deg)= " << theta * 180.0 / M_PI << " (radians)= "
			<< theta << "\n";
	print_run_double(omega_prec, "omega_prec= ", &outfile);
	print_run_double(omega_3, "omega_3= ", &outfile);

	// Rotate resolved body so that angular momentum is up
	// Replace with rotate_to_principal???????
	L = measure_L(n_body_sim, i_low, i_high);
	rotate_body(n_body_sim, i_low, i_high, 0.0, -atan2(L.getZ(), L.getY()), 0);
	L = measure_L(n_body_sim, i_low, i_high);
	rotate_body(n_body_sim, i_low, i_high, -atan2(L.getY(), L.getX()), 0, 0);
	L = measure_L(n_body_sim, i_low, i_high);
	rotate_body(n_body_sim, i_low, i_high, 0, 0, M_PI / 2.0);
	L = measure_L(n_body_sim, i_low, i_high);
	std::cout << "llx lly llz = " << L.getX() << " " << L.getY() << " "
			<< L.getZ() << "\n";

	// Bar chi???????
	// this is probably not correct if obliquity is greater than pi/2
	double barchi = 1.0 * abs(omega_prec) * tau_relax; // initial value of barchi
	print_run_double(barchi, "barchi", &outfile);

	// ratio of numbers of particles to numbers of springs for resolved body
	double Nratio = (double) num_springs / (i_high - i_low);
	std::cout << "N= " << n_body_sim->N << " NS=" << num_springs << " NS/N="
			<< Nratio << "\n";
	outfile << "N= " << n_body_sim->N << " NS=" << num_springs << " NS/N="
			<< Nratio << std::endl;
	outfile.close();

	// Set final rebound info
	reb_springs(n_body_sim); // pass spring index list to display
	n_body_sim->heartbeat = heartbeat; // set up integration

	// Recenter simulation
	subtract_cov(n_body_sim, i_low, i_high);
	center_sim(n_body_sim, i_low, i_high);

	// Integrate simulation
	if (t_max == 0.0) {
		reb_integrate(n_body_sim, INFINITY);
	} else {
		reb_integrate(n_body_sim, t_max);
	}
}

// Things to do every timestep
void heartbeat(reb_simulation *const n_body_sim) {
	static string extendedfile = fileroot + "_ext.txt";

	// Output simulation timing info every 10 timesteps
	if (reb_output_check(n_body_sim, 10.0 * n_body_sim->dt)) {
		reb_output_timing(n_body_sim, 0);
	}

	// After t_damp, end extra damping
	if (abs(n_body_sim->t - t_damp) < 0.9 * n_body_sim->dt) {
		set_gamma(def_gamma);
	}

	// Recenter simulation
	center_sim(n_body_sim, 0, n_body_sim->N - num_perts);

	// Write orbital info of resolved body
	if (reb_output_check(n_body_sim, print_interval)) {
		write_resolved_with_E(n_body_sim, 0, n_body_sim->N, extendedfile, 0.0);
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

void additional_forces(reb_simulation *n_body_sim) {
	// Spring forces
	spring_forces(n_body_sim);
}
