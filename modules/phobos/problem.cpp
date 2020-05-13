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
const int RAND = 0, HCP = 1, ELLIPSE = 2;

// Migration timescales
vector<double> inv_tau_a, inv_tau_e, migr_time;

// Global scales
double mass_scale, time_scale, length_scale, temp_scale, omega_scale, vel_scale,
		p_scale, L_scale, a_scale, F_scale, E_scale, dEdt_scale, P_scale;

// Forward declarations
void heartbeat(reb_simulation *const n_body_sim);
void reb_springs(reb_simulation *const n_body_sim);
void additional_forces(reb_simulation *n_body_sim);

int main(int argc, char *argv[]) {

	// Read in configuration file
	Config cfg;
	cfg.readFile("problem.cfg");

	// Get scales
	read_scales(&cfg);

	// Vars read in
	double t_max, dt, max_spring_dist, gamma_fac, k, omega_in[3], min_part_dist,
			obliquity, ratio1, ratio2, m_ball, r_ball;
	int lattice_type, num_point_masses;

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
			&& cfg.lookupValue("lattice_type", lattice_type)
			&& cfg.lookupValue("min_particle_distance", min_part_dist)
			&& cfg.lookupValue("obliquity", obliquity)
			&& cfg.lookupValue("ratio1", ratio1)
			&& cfg.lookupValue("ratio2", ratio2)
			&& cfg.lookupValue("J2_p", J2_plus)
			&& cfg.lookupValue("R_p", R_plus)
			&& cfg.lookupValue("theta_p", theta_plus)
			&& cfg.lookupValue("phi_p", phi_plus))) {
		throw "Failed to read in problem.cfg. Exiting.";
	} else {
		std::cout << "Read in problem.cfg." << std::endl;
	}
	Vector omega(omega_in);

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
	inv_tau_a.resize(num_point_masses, 0.0);
	inv_tau_e.resize(num_point_masses, 0.0);
	migr_time.resize(num_point_masses, 0.0);

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
		point_masses[i].lookupValue("inv_axis_drift", inv_tau_a[i]);
		point_masses[i].lookupValue("in_e_drift", inv_tau_e[i]);
		point_masses[i].lookupValue("drift_fold_time", migr_time[i]);
	}

	obliquity *= (M_PI / 180.0);

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
	switch (lattice_type) {
	case RAND: {
		// rand_football_from_sphere(r,b_distance,r_ball,r_ball*ratio1, r_ball*ratio2,mball );
		rand_ellipsoid(n_body_sim, min_part_dist, r_ball, r_ball * ratio1,
				r_ball * ratio2, m_ball);
		break;
	}
	case HCP: {
		hcp_ellipsoid(n_body_sim, min_part_dist, r_ball, r_ball * ratio1,
				r_ball * ratio2, m_ball);
		break;
	}
	case ELLIPSE: {
		cubic_ellipsoid(n_body_sim, min_part_dist, r_ball, r_ball * ratio1,
				r_ball * ratio2, m_ball);
		break;
	}
	default: {
		throw "No lattice of type selected. Exiting.";
	}
	}

	// Select indices of resolved body
	int i_low = 0;
	int i_high = n_body_sim->N;

	// Move reference frame to resolved body and rotate
	subtract_com(n_body_sim, i_low, i_high);
	rotate_to_principal(n_body_sim, i_low, i_high);
	subtract_cov(n_body_sim, i_low, i_high);

	// Spin the body
	spin_body(n_body_sim, i_low, i_high, omega);

	// Spin period
	double speriod = abs(2.0 * M_PI / omega.getZ());
	print_run_double(speriod, "spin period", &outfile);

	// Tilt body by obliquity
	rotate_body(n_body_sim, i_low, i_high, 0.0, obliquity, 0.0);

	// Connect all particles within max_spring_dist by hot springs
	connect_springs_dist(n_body_sim, max_spring_dist, 0, n_body_sim->N,
			def_spring);

	// Output Young's modulus of mesh
	double ddr = 0.4;
	print_run_double(max_spring_dist, "max spring length", &outfile);
	print_run_double(ddr, "ddr", &outfile);
	print_run_double(Young_mesh(n_body_sim, i_low, i_high, 0.0, ddr),
			"Young's modulus hot", &outfile);
	print_run_double(Young_full_mesh(), "Young's modulus big hot", &outfile);
	print_run_double(mean_spring_length(), "Mean spring length", &outfile);

	// Calculate Kelvin-Voigt relaxation time
	// Factor of 0.5 is due to reduced mass being used in calculation
	double tau_relax = 1.0 * def_gamma * 0.5
			* (m_ball / (n_body_sim->N - num_perts)) / def_spring.k;
	print_run_double(tau_relax, "relaxation time", &outfile);

	// Set up all point masses
	double omega_orb = 0.0;
	if (num_point_masses > 0) {
		for (int i = 0; i < num_point_masses; i++) {
			omega_orb = add_pt_mass_kep(n_body_sim, i_low, i_high, i_central,
					masses[i], radii[i], orb_els[i]);
			outfile << std::setprecision(3) << "resbody mm= " << omega_orb
					<< std::setprecision(2) << " period= "
					<< 2.0 * M_PI / omega_orb << "\n";
			std::cout << "resbody mm=" << omega_orb << " period="
					<< 2.0 * M_PI / omega_orb << "\n";
			i_central = i_high;
		}
		num_perts = num_point_masses;
	}
	std::cout << std::setprecision(3) << "varpi precession fac = "
			<< 1.5 * J2_plus * pow(R_plus / orb_els[0].a, 2.0) << "\n";

	// Bar chi???????
	// this is probably not correct if obliquity is greater than pi/2
	double barchi = 2.0 * abs(omega_orb - omega.getZ()) * tau_relax; // initial value of barchi
	double posc = 0.5 * 2.0 * M_PI / abs(omega_orb - omega.getZ()); // for oscillations!
	print_run_double(barchi, "barchi", &outfile);
	print_run_double(posc, "posc", &outfile);

	// ratio of numbers of particles to numbers of springs for resolved body
	double Nratio = (double) num_springs / (i_high - i_low);
	std::cout << "N= " << n_body_sim->N << " NS=" << num_springs << " NS/N="
			<< Nratio << "\n";
	outfile << "N= " << n_body_sim->N << " NS=" << num_springs << " NS/N="
			<< Nratio << "\n";
	outfile.close();

	// Set final rebound info
	reb_springs(n_body_sim); // pass spring index list to display
	n_body_sim->heartbeat = heartbeat; // set up integration

	// Recenter reference frame
	center_sim(n_body_sim, i_low, i_high);

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
	static vector<string> pointmassfiles;

	// First time through, make file names
	static bool first = true;
	if (first) {
		first = false;
		for (int i = 0; i < num_perts; i++) {
			pointmassfiles.push_back(
					fileroot + "_pm" + std::to_string(i) + ".txt");
		}
	}

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

	// Drift the orbits
	// Get resolved body index range
	int i_low = 0;
	int i_high = n_body_sim->N - num_perts;
	if (num_perts > 0) {
		for (int i = 0; i < num_perts; i++) {
			// Note rebound number of drifting particle
			int i_part = i_central + i;

			// Get current state of migration
			double migfac = exp(-n_body_sim->t / migr_time[i]);

			// Drift resolved body around central mass
			if (i == 0) {
				drift_resolved(n_body_sim, n_body_sim->dt,
						inv_tau_a[i] * migfac, inv_tau_e[i] * migfac, i_central,
						i_low, i_high);
				// Otherwise drift point mass with respect to central mass
			} else {
				drift_bin(n_body_sim, n_body_sim->dt, inv_tau_a[i] * migfac,
						inv_tau_e[i] * migfac, i_central, i_part);
			}
		}
	}

	// Save total power dissipation across calls
	static double dEdt_sum = 0.0;
	dEdt_sum += dEdt_total(n_body_sim);

	// Output resolved body info every print_interval with energy info
	if (reb_output_check(n_body_sim, print_interval)) {
		// Get average power dissipation
		double dEdt_ave = dEdt_sum * n_body_sim->dt / print_interval;

		// Write out
		write_resolved_with_E(n_body_sim, 0, n_body_sim->N - num_perts,
				extendedfile, dEdt_ave);

		// Reset total power
		dEdt_sum = 0.0;

		// Write out point masses
		if (num_perts > 0) {
			for (int i = 0; i < num_perts; i++) {
				int ip = i_central + i;
				write_pt_mass(n_body_sim, ip, i, pointmassfiles[i]);
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

// Additional accelerations
void additional_forces(reb_simulation *n_body_sim) {
	// Spring forces
	spring_forces(n_body_sim);

	// Quadrupole force from primary
	quadrupole_accel(n_body_sim, J2_plus, R_plus, phi_plus, theta_plus,
			n_body_sim->N - num_perts);
}
