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
 * This lets you read in the Bennu shape model and give it an impact
 * no external point mass particles in this routine
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
#include "gravity.h"
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

extern vector<node> nodes;

// Default global values
int num_perts = 0;		// Number of perturbers
double def_gamma, t_damp, print_interval, lat_pulse, long_pulse, dt_pulse,
		force_pulse, dangle_pulse;		// Read in - see cfg
string fileroot;		// Output file base name

// Shape model selectors
const int BENNU = 0, ELLIPSOID = 1, CONE = 2;

// Global scales
double mass_scale, time_scale, length_scale, temp_scale, omega_scale, vel_scale,
		p_scale, L_scale, a_scale, F_scale, E_scale, dEdt_scale, P_scale;

// Forward declarations
void reset_surface_radii(reb_simulation *const n_body_sim, double b_distance);
void heartbeat(reb_simulation *const n_body_sim);
void reb_springs(reb_simulation *const n_body_sim);
void apply_impact(reb_simulation *const n_body_sim, double lat, double longi,
		double dangle, double famp, double tau);
void additional_forces(reb_simulation *n_body_sim);

int main(int argc, char *argv[]) {

	// Read in configuration file
	Config cfg;
	cfg.readFile("problem.cfg");

	// Get scales
	read_scales(&cfg);

	// Vars read in
	double t_max, dt, max_spring_dist, gamma_fac, k, omega_in[3], lat_pulse,
			long_pulse, min_part_dist, dangle_pulse, force_pulse, dt_pulse,
			surf_dist, obliquity, ratio1, ratio2, m_ball;
	int surface_type;

	if (!(cfg.lookupValue("fileroot", fileroot) && cfg.lookupValue("dt", dt)
			&& cfg.lookupValue("t_max", t_max)
			&& cfg.lookupValue("print_interval", print_interval)
			&& cfg.lookupValue("max_spring_dist", max_spring_dist)
			&& cfg.lookupValue("k", k) && cfg.lookupValue("gamma", def_gamma)
			&& cfg.lookupValue("damp_fac", gamma_fac)
			&& cfg.lookupValue("t_damp", t_damp)
			&& cfg.lookupValue("mass_ball", m_ball)
			&& cfg.lookupValue("omega_x", omega_in[0])
			&& cfg.lookupValue("omega_y", omega_in[1])
			&& cfg.lookupValue("omega_z", omega_in[2])
			&& cfg.lookupValue("surface_type", surface_type)
			&& cfg.lookupValue("min_particle_distance", min_part_dist)
			&& cfg.lookupValue("surface_dist", surf_dist)
			&& cfg.lookupValue("obliquity", obliquity)
			&& cfg.lookupValue("lat_impact", lat_pulse)
			&& cfg.lookupValue("long_impact", long_pulse)
			&& cfg.lookupValue("dangle_impact", dangle_pulse)
			&& cfg.lookupValue("force_impact", force_pulse)
			&& cfg.lookupValue("dt_impact", dt_pulse))) {
		throw "Failed to read in problem.cfg. Exiting.";
	} else {
		std::cout << "Read in problem.cfg." << std::endl;
	}
	Vector omega(omega_in);

	if (surface_type == 1) {
		if (!(cfg.lookupValue("ratio1", ratio1)
				&& cfg.lookupValue("ratio2", ratio2))) {
			throw "Failed to read ellipse parameters from problem.cfg. Exiting.";
		} else {
			std::cout << "Read ellipse parameters from problem.cfg.";
		}
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

	// Convert input in degrees to radians
	obliquity *= (M_PI / 180.0);
	lat_pulse *= (M_PI / 180.0);
	long_pulse *= (M_PI / 180.0);
	dangle_pulse *= (M_PI / 180.0);

	// Set up default spring parameters
	spring def_spring;
	def_spring.gamma = gamma_fac * def_gamma; // initial damping coefficient
	def_spring.k = k; // spring constant

	// Output spring and particle distances
	std::cout << std::setprecision(2) << "Max spring length = "
			<< max_spring_dist << " Min interparticle distance = "
			<< min_part_dist << "\n";

	// Use Bennu shape model
	switch (surface_type) {
	case BENNU: {

		// Read in shape file (in units of km)
		string filename = "101955bennu.tab";
		read_vertices(n_body_sim, filename);

		// Number of particles that make up Bennu model
		int N_bennu = n_body_sim->N;
		std::cout << filename << " read in.\n";

		// Correct units from km to length_scale
		stretch(n_body_sim, 0, n_body_sim->N, 1.0 / length_scale);

		// Fill shape with particles
		rand_shape(n_body_sim, min_part_dist, 1.0);
		std::cout << "Shape filled.\n";

		// Find surface particles
		mark_surf_shrink_int_shape(n_body_sim, 0, N_bennu, surf_dist);

		// Clean up shape vertices
		rm_particles(n_body_sim, 0, N_bennu);
		std::cout << "Shape vertices deleted.\n";
		break;
	}
		// Random ellipsoid
	case ELLIPSOID: {
		// Scale radius
		double r_ball = 1.0;

		// Neglecting 4pi/3 factor
		// Vol of a triax ellipsoid is abc*4pi/3
		double volume_ratio = pow(r_ball, 3.0) * ratio1 * ratio2;

		// Determine equivalent spherical radius
		double vol_radius = pow(volume_ratio, 1.0 / 3.0);

		// Compute actual semi-major axis
		r_ball /= vol_radius;

		// Create random ellipsoide and mark its surface
		rand_ellipsoid(n_body_sim, min_part_dist, r_ball, r_ball * ratio1,
				r_ball * ratio2, m_ball);
		init_nodes(n_body_sim, 0, 0);
		mark_surf_shrink_int_ellipsoid(n_body_sim, surf_dist, r_ball,
				r_ball * ratio1, r_ball * ratio2);
		break;
	}
		// Random cone
	case CONE: {
		// Set radius of cone
		double r_cone = 1.0;

		// Vol of a cone is 2 pi/3 r^2 h
		double volume_ratio = pow(r_cone, 3.0) * ratio1 * 0.5;

		// Get actual radius of cone
		double vol_radius = pow(volume_ratio, 1.0 / 3.0);
		r_cone /= vol_radius;

		// Create random cone particles and mark surface particles
		rand_cone(n_body_sim, min_part_dist, r_cone, r_cone * ratio1, m_ball);
		init_nodes(n_body_sim, 0, 0);
		mark_surf_shrink_int_cone(n_body_sim, surf_dist, r_cone,
				r_cone * ratio1);
		break;
	}
	default: {
		throw "No existing shape model selected. Exiting.";
	}
	}

	// Get resolved body indices
	int i_low = 0;
	int i_high = n_body_sim->N;

	// Move reference frame to center of mass frame
	subtract_com(n_body_sim, i_low, i_high);
	subtract_cov(n_body_sim, i_low, i_high);

	// Spin the body
	spin_body(n_body_sim, i_low, i_high, omega);

	// Required??????
	//subtract_cov(n_body_sim, i_low, i_high); // subtract center of velocity

	// Set output filename
	string filename = fileroot + "_run.txt";
	std::ofstream outfile(filename, std::ios::out | std::ios::app);

	// Output spin info
	double speriod = abs(2.0 * M_PI / omega.getZ());
	std::cout << std::setprecision(6) << "spin period " << speriod << "\n";
	outfile << std::setprecision(6) << "spin period " << speriod << "\n";
	outfile << std::setprecision(6) << "omegaz " << omega.getZ() << "\n";

	// Tilt by obliquity
	if (obliquity != 0.0)
		rotate_body(n_body_sim, i_low, i_high, 0.0, obliquity, 0.0);

	// Connect particles within max_spring_dist by springs
	connect_springs_dist(n_body_sim, max_spring_dist, 0, n_body_sim->N,
			def_spring);
	std::cout << "Springs connected.\n";

	// Set display properties
	double max_r = max_radius(n_body_sim, i_low, i_high);
	const double boxsize = 2.5 * max_r;
	reb_configure_box(n_body_sim, boxsize, 1, 1, 1);

	// Output Young's modulus
	double ddr = 0.6 * max_r;
	double E_mesh = Young_mesh(n_body_sim, i_low, i_high, 0.0, ddr);
	std::cout << std::setprecision(6) << "Young's modulus " << E_mesh << "\n";
	outfile << std::setprecision(6) << "Young's modulus " << E_mesh << "\n";
	std::cout << std::setprecision(3) << "ddr = " << ddr << " mesh_distance = "
			<< max_spring_dist << "\n";
	outfile << "max spring length " << max_spring_dist << "\n";

	// Output mean spring length
	double LL = mean_spring_length();
	std::cout << std::setprecision(4) << "mean L = " << LL << "\n";
	outfile << std::setprecision(4) << "mean_L " << LL << "\n";

	// Calculate Kelvin-Voigt relaxation time
	// Factor of 0.5 is due to reduced mass being used in calculation
	double tau_relax = 1.0 * def_gamma * 0.5 * (m_ball / (n_body_sim->N - 1))
			/ def_spring.k;
	std::cout << std::setprecision(3) << "relaxation time " << tau_relax
			<< "\n";
	outfile << std::setprecision(3) << "relaxation_time " << tau_relax << "\n";

	// Get ratio of springs to particles
	double N_ratio = (double) num_springs / (double) n_body_sim->N;
	std::cout << "N=" << n_body_sim->N << " NS=" << num_springs << " NS/N="
			<< N_ratio << "\n";
	outfile << "N=" << n_body_sim->N << " NS=" << num_springs << " NS/N="
			<< N_ratio << "\n";

	// Get number of surface particles
	int Nsurf = 0;
	for (int i = 0; i < n_body_sim->N; i++) {
		if (nodes[i].is_surf) {
			Nsurf++;
		}
	}
	std::cout << "Nsurf=" << Nsurf << std::endl;
	outfile << "Nsurf=" << Nsurf << "\n";

	outfile.close();

	// Set final rebound info
	reb_springs(n_body_sim); // pass spring index list to display
	n_body_sim->heartbeat = heartbeat; // set up integration

	// Calculate surface gravity alone at beginning of simulation
	reb_calculate_acceleration(n_body_sim);

	// Output surface particle info at first timestep
	filename = fileroot + "_surf_nosprings.txt";
	write_surf_part(n_body_sim, 0, n_body_sim->N - num_perts, filename);

	// Integrate simulation
	if (t_max == 0.0) {
		reb_integrate(n_body_sim, INFINITY);
	} else {
		reb_integrate(n_body_sim, t_max);
	}
}

// Stuff done every timestep
void heartbeat(reb_simulation *const n_body_sim) {

	// Output timing info every 10 timesteps
	if (reb_output_check(n_body_sim, 10.0 * n_body_sim->dt)) {
		reb_output_timing(n_body_sim, 0);
	}

	// Turn off damping at t_damp
	if (abs(n_body_sim->t - t_damp) < 0.9 * n_body_sim->dt) {
		set_gamma(def_gamma);
	}

	// Output after t_damp
	if (n_body_sim->t > t_damp) {

		// Output file every print_interval
		if (reb_output_check(n_body_sim, print_interval * n_body_sim->dt)) {
			string filename = fileroot
					+ zero_pad_int(5, (int) n_body_sim->dt / print_interval);
			write_surf_part(n_body_sim, 0, n_body_sim->N - num_perts, filename);
		}
	}

	// Recenter body
	subtract_com(n_body_sim, 0, n_body_sim->N - num_perts);

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

// Simulate impact
// Push on surface particles radially within dangle of center at lat,long (in radians)
// Total force distributed among all particles in the angular distance
// Meant to be called for a specific length of time (see additional_forces)
// Currently applies cosine for force profile over time
void apply_impact(reb_simulation *const n_body_sim, double lat, double longi,
		double dangle, double force, double tau) {

	// Get particle info
	reb_particle *particles = n_body_sim->particles;

	// Store start time of pulse and initial display radius
	static bool first = true;
	static double t0 = 0.0;
	static double r0 = 0.0;
	if (first) {
		first = false;
		t0 = n_body_sim->t;
		for (int i = 0; i < n_body_sim->N - num_perts; i++) {
			if (nodes[i].is_surf) {
				r0 = n_body_sim->particles[i].r;
				break;
			}
		}
	}

	// Get direction pointing at center of pulse
	Vector x0_hat = { cos(longi) * cos(lat), sin(longi) * cos(lat), sin(lat) };

	// Time since start of pulse
	double t = n_body_sim->t - t0;

	// Figure out how many particles we will push, so we can normalize force
	int num_pushed = 0;
	for (int i = 0; i < n_body_sim->N - num_perts; i++) {

		// If particle is on surface
		if (nodes[i].is_surf) {
			// Get particle location
			Vector x = { particles[i].x, particles[i].y, particles[i].z };
			Vector x_hat = x / x.len();

			// Get angle between particle and center of pulse
			double projection = dot(x0_hat, x_hat);
			double ang_dist = acos(projection); // Is positive in [0,pi]

			// If within dangle, count
			if (ang_dist < dangle) {
				num_pushed++;
			}
		}
	}

	std::cout << "apply_impact: Pushing on " << num_pushed << " particles."
			<< std::endl;

	// Apply force
	for (int i = 0; i < n_body_sim->N - num_perts; i++) {

		// If particle is on surface
		if (nodes[i].is_surf) {
			// Get particle info
			double m = particles[i].m;
			Vector x = { particles[i].x, particles[i].y, particles[i].z };
			Vector x_hat = x / x.len();

			// Get angle between particle and center of pulse
			double projection = dot(x0_hat, x_hat);
			double ang_dist = acos(projection); // Is positive in [0,pi]

			// Get acceleration
			Vector a = force * x_hat / (m * num_pushed * x.len()); // dv = F*dt/m  with total force applied F =

			// Smooth pulse by cos
			a *= 1.0 - cos(2.0 * M_PI * t / tau);

			// Apply force radially
			if (ang_dist < dangle) {
				n_body_sim->particles[i].ax -= a.getX();
				n_body_sim->particles[i].ay -= a.getY();
				n_body_sim->particles[i].az -= a.getZ();

				// Set display radius so we can see where we push
				if (t < tau - 2.1 * n_body_sim->dt) {
					n_body_sim->particles[i].r = 0.09;
				} else {
					n_body_sim->particles[i].r = r0;
				}
			}
		}
	}

	// Ensure display radius is restored at end of pulse
	if (t >= tau - 2.1 * n_body_sim->dt) {
		for (int i = 0; i < n_body_sim->N - num_perts; i++) {
			if (nodes[i].is_surf)
				particles[i].r = r0;
		}
	}
}

void additional_forces(reb_simulation *n_body_sim) {
	// Spring forces
	spring_forces(n_body_sim);

	// Impact
	if ((n_body_sim->t > t_damp) && (n_body_sim->t <= t_damp + dt_pulse)) {
		apply_impact(n_body_sim, lat_pulse, long_pulse, dangle_pulse,
				force_pulse, dt_pulse); // pulse on surface
	}
}

// Reset surface particle radii post-impact
void reset_surface_radii(reb_simulation *const n_body_sim, double b_distance) {
	for (int i = 0; i < n_body_sim->N; i++) {
		if (nodes[i].is_surf)
			n_body_sim->particles[i].r = b_distance / 4;
	}
}
