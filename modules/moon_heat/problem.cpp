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
#include "springs.h"
#include "physics.h"
#include "stress.h"
#include "shapes.h"
#include "heat.h"
#include "orb.h"

using namespace libconfig;

using std::string;
using std::vector;

// Global values
int num_springs = 0;	// Global numbers of springs
vector<spring> springs;	// Global spring array

extern vector<double> tot_power;

// Default global values
int num_perts = 0;		// Number of perturbers
double gamma_fac, t_damp, print_interval, heat_print_interval, power_fac;// Read in - see cfg
string fileroot;		// Output file base name

// Central mass location
// Should be set to -1 initially
int i_central = -1;

// Lattice selectors
const int RAND = 0, HCP = 1;

// Global scales
double mass_scale, time_scale, length_scale, temp_scale, omega_scale, vel_scale,
		p_scale, L_scale, a_scale, F_scale, E_scale, dEdt_scale, P_scale;

// Forward declarations
void heartbeat(reb_simulation *const n_body_sim);
void print_run_double(double quantity, string label, std::ofstream *file);
void reb_springs(reb_simulation *const r);
void additional_forces(reb_simulation *n_body_sim);

int main(int argc, char *argv[]) {

	// Read in configuration file
	Config cfg;
	cfg.readFile("problem.cfg");

	// Get scales
	read_scales(&cfg);

	// Vars read in
	double t_max, dt, max_spring_dist, gamma_hot, gamma_cold, gamma_fac, k_hot,
			k_cold, k_heat_hot, k_heat_cold, omega_in[3], min_part_dist,
			surf_dist, obliquity, ratio1, ratio2, m_ball, r_ball, freq, phi0,
			T_surf, T_int, R_shell, ratio1_shell, ratio2_shell, x_shell,
			y_shell, z_shell, cv_hot, cv_cold, k_fac, rel_dens;
	int lattice_type, num_point_masses;

	if (!(cfg.lookupValue("fileroot", fileroot) && cfg.lookupValue("dt", dt)
			&& cfg.lookupValue("t_max", t_max)
			&& cfg.lookupValue("print_interval", print_interval)
			&& cfg.lookupValue("heat_print_interval", heat_print_interval)
			&& cfg.lookupValue("max_spring_dist", max_spring_dist)
			&& cfg.lookupValue("k_hot", k_hot)
			&& cfg.lookupValue("k_cold", k_cold)
			&& cfg.lookupValue("k_heat_hot", k_heat_hot)
			&& cfg.lookupValue("k_heat_cold", k_heat_cold)
			&& cfg.lookupValue("gamma_hot", gamma_hot)
			&& cfg.lookupValue("gamma_cold", gamma_cold)
			&& cfg.lookupValue("damp_fac", gamma_fac)
			&& cfg.lookupValue("power_fac", power_fac)
			&& cfg.lookupValue("t_damp", t_damp)
			&& cfg.lookupValue("mass_ball", m_ball)
			&& cfg.lookupValue("radius_ball", r_ball)
			&& cfg.lookupValue("omega_x", omega_in[0])
			&& cfg.lookupValue("omega_y", omega_in[1])
			&& cfg.lookupValue("omega_z", omega_in[2])
			&& cfg.lookupValue("lattice_type", lattice_type)
			&& cfg.lookupValue("min_particle_distance", min_part_dist)
			&& cfg.lookupValue("surface_dist", surf_dist)
			&& cfg.lookupValue("obliquity", obliquity)
			&& cfg.lookupValue("max_strength_coeff", freq)
			&& cfg.lookupValue("max_strength_location", phi0)
			&& cfg.lookupValue("T_int", T_int)
			&& cfg.lookupValue("T_surf", T_surf)
			&& cfg.lookupValue("ratio1", ratio1)
			&& cfg.lookupValue("ratio2", ratio2)
			&& cfg.lookupValue("R_shell", R_shell)
			&& cfg.lookupValue("ratio1_shell", ratio1_shell)
			&& cfg.lookupValue("ratio2_shell", ratio2_shell)
			&& cfg.lookupValue("x_shell", x_shell)
			&& cfg.lookupValue("y_shell", y_shell)
			&& cfg.lookupValue("z_shell", z_shell)
			&& cfg.lookupValue("Cv_hot", cv_hot)
			&& cfg.lookupValue("Cv_cold", cv_cold)
			&& cfg.lookupValue("k_factor", k_fac)
			&& cfg.lookupValue("relative_density", rel_dens))) {
		throw "Failed to read in problem.cfg. Exiting.";
	} else {
		std::cout << "Read in problem.cfg." << std::endl;
	}
	Vector omega(omega_in);
	Vector shell_cent = { x_shell, y_shell, z_shell };

	// Read in point masses
	cfg.readFile("point_mass.cfg");

	if (!cfg.lookupValue("num_point_masses", num_point_masses)) {
		throw "Failed to read in number of point masses.";
	} else {
		std::cout << "Reading " << num_point_masses
				<< " point masses from point_mass.cfg." << std::endl;
	}

	// Properties of point masses
	double radii[num_point_masses], masses[num_point_masses];
	OrbitalElements orb_els[num_point_masses];

	Setting &point_masses = cfg.getRoot()["point_masses"];

	// Load into vectors
	for (int i = 0; i < num_point_masses; i++) {
		point_masses[i].lookupValue("radius", radii[i]);
		point_masses[i].lookupValue("mass", masses[i]);
		point_masses[i].lookupValue("semi-major_axis", orb_els[i].a);
		point_masses[i].lookupValue("eccentricity", orb_els[i].e);
		point_masses[i].lookupValue("inclination", orb_els[i].i);
		point_masses[i].lookupValue("long_asc_node", orb_els[i].long_asc_node);
		point_masses[i].lookupValue("arg_periapsis", orb_els[i].arg_peri);
		point_masses[i].lookupValue("mean_anom", orb_els[i].mean_anom);
	}

	// Correct to radians, limit power_fac
	obliquity *= (M_PI / 180.0);
	if (power_fac > 1.0) {
		power_fac = 1.0;
	} else if (power_fac < 0.0) {
		power_fac = 1.0;
	}

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
	n_body_sim->dt = dt;							// Integration timestep
	n_body_sim->softening = min_part_dist / 100.0;// Gravitational softening length

	// Set display window size
	const double boxsize = 3.2 * r_ball;
	reb_configure_box(n_body_sim, boxsize, 1, 1, 1);

	// Spring properties
	spring def_spring_hot;
	def_spring_hot.gamma = gamma_hot * gamma_fac;
	def_spring_hot.k = k_hot;
	def_spring_hot.k_heat = k_heat_hot;

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
	outfile << std::setprecision(6) << "fmfac " << rel_dens << "\n"; //

	// Create particle distribution
	switch (lattice_type) {
	case RAND: {
		rand_ellipsoid(n_body_sim, min_part_dist, r_ball, r_ball * ratio1,
				r_ball * ratio2, m_ball);
		init_nodes(n_body_sim, cv_hot, T_int); // temperatures on nodes
		mark_surf_shrink_int_ellipsoid(n_body_sim, surf_dist, r_ball,
				r_ball * ratio1, r_ball * ratio2);
		init_surf_temp(T_surf);
		break;
	}
	case HCP: {
		hcp_ellipsoid(n_body_sim, min_part_dist, r_ball, r_ball * ratio1,
				r_ball * ratio2, m_ball);
		init_nodes(n_body_sim, cv_hot, T_int); // temperatures on nodes
		mark_surf_shrink_int_ellipsoid(n_body_sim, surf_dist, r_ball,
				r_ball * ratio1, r_ball * ratio2);
		init_surf_temp(T_surf);
		break;
	}
	default: {
		throw "No lattice of type selected. Exiting.";
	}
	}

	// Select indices of resolved body
	int i_low = 0;
	int i_high = n_body_sim->N;

	// Move reference frame to resolved body
	subtract_com(n_body_sim, i_low, i_high);
	subtract_cov(n_body_sim, i_low, i_high);

	// Connect all particles within max_spring_dist by hot springs
	connect_springs_dist(n_body_sim, max_spring_dist, 0, n_body_sim->N,
			def_spring_hot);

	// Set size of total power vector
	tot_power.resize(num_springs, 0.0);

	// Output Young's modulus of mesh
	double ddr = 0.5;
	double Emush = Young_mesh(n_body_sim, i_low, i_high, 0.0, ddr);
	double Emush_big = Young_full_mesh();
	print_run_double(max_spring_dist, "max spring length", &outfile);
	print_run_double(ddr, "ddr", &outfile);
	print_run_double(Emush, "Young's modulus hot", &outfile);
	print_run_double(Emush_big, "Young's modulus big hot", &outfile);
	print_run_double(Emush * k_cold / k_hot, "Young's modulus cold", &outfile);
	print_run_double(Emush_big * k_cold / k_hot, "Young's modulus big cold",
			&outfile);

	// Output thermal conductivity of mesh
	double K_T = therm_cond_mesh(n_body_sim, i_low, i_high, 0.0, ddr);
	print_run_double(K_T, "Thermal conductivity hot", &outfile);
	print_run_double(K_T * k_heat_cold / k_heat_hot,
			"Thermal conductivity cold", &outfile);
	print_run_double(mean_spring_length(), "Mean spring length", &outfile);

	// Set up point masses
	double omega_orb = 0.0;
	if (num_point_masses > 0) {
		for (int i = 0; i < num_point_masses; i++) {

			// Add point mass (returns mean motion)
			omega_orb = add_pt_mass_kep(n_body_sim, i_low, i_high, i_central,
					masses[i], radii[i], orb_els[i]);

			// Output info after adding central mass
			if (i_central == -1) {
				double a = omega_orb * orb_els[i].a;
				double a_dot = 3.0 * masses[i] * a / pow(orb_els[i].a, 5.0); // should approximately be adot
				outfile << std::setprecision(3) << "a_dot " << a_dot << "\n";
				i_central = i_high;
			}

			// Output mean motion, period
			outfile << std::setprecision(3) << "resbody mm= " << omega_orb
					<< std::setprecision(2) << " period= "
					<< 2.0 * M_PI / omega_orb << "\n";
			std::cout << "resbody mm=" << omega_orb << " period="
					<< 2.0 * M_PI / omega_orb << "\n";
		}

		// Update number of perturbers
		num_perts = num_point_masses;
	}

	// Calculate Kelvin-Voigt relaxation time
	// Factor of 0.5 is due to reduced mass being used in calculation
	double tau_relax_hot = 1.0 * gamma_hot * 0.5
			* (m_ball / (n_body_sim->N - num_perts)) / def_spring_hot.k;
	double tau_relax_cold = 1.0 * gamma_cold * 0.5
			* (m_ball / (n_body_sim->N - num_perts)) / k_cold;
	print_run_double(tau_relax_hot, "relaxation time hot", &outfile);
	print_run_double(tau_relax_cold, "relaxation time cold", &outfile);

	// Bar chi??????
	double barchi_hot = 2.0 * abs(omega_orb) * tau_relax_hot;
	double barchi_cold = 2.0 * abs(omega_orb) * tau_relax_cold;
	print_run_double(barchi_hot, "barchi hot", &outfile);
	print_run_double(barchi_cold, "barchi cold", &outfile);

	// Spin period
	double speriod = abs(2.0 * M_PI / omega.getZ());
	print_run_double(speriod, "spin period", &outfile);
	print_run_double(omega.getZ(), "omegaz", &outfile);

	// Ratio of springs to particles
	double Nratio = (double) num_springs / (double) n_body_sim->N;
	std::cout << "N= " << n_body_sim->N << " NS=" << num_springs << " NS/N="
			<< Nratio << "\n";
	outfile << "N= " << n_body_sim->N << " NS=" << num_springs << " NS/N="
			<< Nratio << "\n";
	outfile.close();

	// Spin the body
	spin_body(n_body_sim, i_low, i_high, omega);

	// Rotate orbit of body
	Vector L = measure_L_origin(n_body_sim, 0, n_body_sim->N);
	std::cout << std::setprecision(6) << "llx lly llz = " << L.getX() << " "
			<< L.getY() << " " << L.getZ() << "\n";

	rotate_origin(n_body_sim, 0, n_body_sim->N, 0, -atan2(L.getY(), L.getZ()),
			0);
	L = measure_L_origin(n_body_sim, 0, n_body_sim->N);
	std::cout << std::setprecision(6) << "llx lly llz = " << L.getX() << " "
			<< L.getY() << " " << L.getZ() << "\n";

	rotate_origin(n_body_sim, 0, n_body_sim->N, 0, 0, M_PI / 2);
	L = measure_L_origin(n_body_sim, 0, n_body_sim->N);
	std::cout << std::setprecision(6) << "llx lly llz = " << L.getX() << " "
			<< L.getY() << " " << L.getZ() << "\n";

	rotate_origin(n_body_sim, 0, n_body_sim->N, 0, -atan2(L.getY(), L.getZ()),
			0);
	L = measure_L_origin(n_body_sim, 0, n_body_sim->N);
	std::cout << std::setprecision(6) << "llx lly llz = " << L.getX() << " "
			<< L.getY() << " " << L.getZ() << "\n";

	// Tilt body by obliquity
	if (obliquity != 0.0) {
		rotate_body(n_body_sim, i_low, i_high, 0.0, obliquity, 0.0);
	}

	// Set properties of exterior of shell
	if (R_shell > 0) {
		// Set spring properties of exterior of shell to cold spring props
		adjust_spring_props_ellipsoid(n_body_sim, k_cold,
				gamma_cold * gamma_fac, k_heat_cold, R_shell,
				ratio1_shell * R_shell, ratio2_shell * R_shell, shell_cent,
				false);

		// Set shell to new relative density
		if (abs(rel_dens) > 1e-3)
			adjust_mass_ellipsoid(n_body_sim, rel_dens, R_shell,
					ratio1_shell * R_shell, ratio2_shell * R_shell, shell_cent,
					false);

		// Soften springs lopsidedly
		if (abs(k_fac) > 1e-3) {
			adjust_spring_props_ellipsoid_phase(n_body_sim, k_fac, 0.0, 0.0,
					freq, phi0, R_shell, ratio1_shell * R_shell,
					ratio2_shell * R_shell, shell_cent, false);
		}
	}

	// Set final rebound info
	reb_springs(n_body_sim); // pass spring index list to display
	n_body_sim->heartbeat = heartbeat; // set up integration

	// Recenter reference frame again
	subtract_com(n_body_sim, 0, n_body_sim->N - num_perts);

	// Integrate simulation
	if (t_max == 0.0) {
		reb_integrate(n_body_sim, INFINITY);
	} else {
		reb_integrate(n_body_sim, t_max);
	}
}

// Things to do every time step
void heartbeat(reb_simulation *const n_body_sim) {
	// Save file names
	static string extendedfile = fileroot + "_ext.txt";
	static string twopfile = fileroot + "_2p.txt";
	static vector<string> pointmassfiles;

	// First time through, make file names
	static bool first = true;
	if (first) {
		first = false;
		for (int i = 0; i < num_perts; i++) {
			pointmassfiles.push_back(
					fileroot + "_pm" + std::to_string(i) + ".txt");
		}

		// Recenter frame on center of mass
		subtract_com(n_body_sim, 0, n_body_sim->N - num_perts);
	}

	// Output simulation timing info every 10 timesteps
	if (reb_output_check(n_body_sim, 10.0 * n_body_sim->dt)) {
		reb_output_timing(n_body_sim, 0);
	}

	// After t_damp, end extra damping
	if (abs(n_body_sim->t - t_damp) < 0.9 * n_body_sim->dt) {
		divide_gamma(gamma_fac);
	}

	// Recenter frame on center of mass
	subtract_com(n_body_sim, 0, n_body_sim->N - num_perts);

	// After damping ends, do heat stuff
	if (n_body_sim->t > t_damp) {

		// Store heat accumulated during each timestep
		rec_power(n_body_sim);

		// Heat nodes with tidal heating (doesn't make use of recorded heat above)
		heat_int_nodes_tidal(n_body_sim, n_body_sim->dt);

		// Transport heat across springs
		transport_heat(n_body_sim, n_body_sim->dt);
	}

	// At every heat_print_interval, print out heating info
	if (reb_output_check(n_body_sim, heat_print_interval)) {
		// Get number of timesteps since last print
		int ndt = (int) (heat_print_interval / n_body_sim->dt);

		// Make filename
		string hfile = heat_filename(n_body_sim, fileroot, heat_print_interval);

		// Write info to file
		write_heat(n_body_sim, hfile, ndt, power_fac);

		// Make filename
		string nfile = node_filename(n_body_sim, fileroot, heat_print_interval);

		// Write node info to file
		write_nodes(n_body_sim, nfile);
	}

	// Every print_interval, write out info about resolved body and point masses
	if (reb_output_check(n_body_sim, print_interval)) {

		// Orbital info
		write_resolved_no_E(n_body_sim, 0, n_body_sim->N - num_perts,
				extendedfile); // orbital info and stuff

		// Nodes with largest x and z values
		write_resolved_2nodes(n_body_sim, 0, n_body_sim->N - num_perts,
				twopfile);

		// Point mass info
		if (num_perts > 0) {
			for (int i = 0; i < num_perts; i++) {
				int i_part = i_central + i;
				write_pt_mass(n_body_sim, i_part, i, pointmassfiles[i]);
			}
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

void additional_forces(reb_simulation *n_body_sim) {
	// Spring forces
	spring_forces(n_body_sim);
}
