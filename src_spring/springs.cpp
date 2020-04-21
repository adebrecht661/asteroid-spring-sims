#include <cmath>
#include <iostream>
extern "C" {
#include "rebound.h"
}
#include "stress.h"
#include "physics.h"
#include "matrix_math.h"
#include "springs.h"

extern int num_springs; // Current number of springs
extern stress_tensor *stresses; // Array of stresses for each particle
int NSmax = 0; // Max number of springs (size of array)

#define L_EPS 1e-6; // Softening for spring length

/********************/
/* Spring functions */
/********************/

// Delete spring at spring_index
void del_spring(struct reb_simulation *const n_body_sim, int spring_index) {
	// Overwrite spring at index with last spring, and reduce count (last spring is still present at original location)
	if (num_springs > 0) {
		springs[spring_index] = springs[num_springs - 1];
		num_springs--;
	}
}

// Add a spring connecting particle_1 and particle_2
// Check that a spring doesn't already exist between these two particles
// Set the natural distance of the spring to the current interparticle distance
// Spring constant is not scaled
// Return -1 if no spring created because you're trying to connect a particle to itself
int add_spring(struct reb_simulation *const n_body_sim, int particle_1,
		int particle_2, spring spr) {
	int particle_low, particle_high;

	// Don't add spring if vertices are the same
	if (particle_1 == particle_2)
		return -1;

	// Place indices in order, for convenience
	if (particle_2 < particle_1) {
		particle_low = particle_2;
		particle_high = particle_1; // order of indices
	} else {
		particle_low = particle_1;
		particle_high = particle_2;
	}

	// Check if these two particles are already connected.
	for (int spring_index = 0; spring_index < num_springs; spring_index++) {
		if ((springs[spring_index].particle_1 == particle_low)
				&& (springs[spring_index].particle_2 == particle_high))
			return spring_index;
	}

	// No spring connects these two indices. Create one.
	spr.particle_1 = particle_low;
	spr.particle_2 = particle_high;
	spr.rs0 = spring_r(n_body_sim, spr).len(); // rest spring length
	add_spring_helper(spr);
	return num_springs - 1; // index of new spring!
}

// Helper function that adds a spring to array
// Caution: no sanity checking
void add_spring_helper(spring spr) {
	// If max size is smaller than current size, increase
	while (NSmax <= num_springs) {
		NSmax += 128;
		springs = (spring*) realloc(springs, sizeof(spring) * NSmax);
	}

	// Add spring to end of array, update count
	springs[num_springs] = spr;
	num_springs++;
}

// Connect springs to all particles with interparticle distances less than max_dist apart for particle index range [i_low,i_high-1]
// Spring is added with rest length at current length
void connect_springs_dist(struct reb_simulation *const n_body_sim,
		double max_dist, int i_low, int i_high, spring spr) {
	// Get particle info from sim
	reb_particle *particles = n_body_sim->particles;

	// Return if low index is above high index
	if (i_high <= i_low)
		return;

	// Create all springs for near neighbors
	// Scan through particles from i_low to i_high-1 (i_high can't pair with itself)
	for (int i = i_low; i < i_high - 1; i++) {
		Vector x_i = { particles[i].x, particles[i].y, particles[i].z };

		// Scan through all possible pairs, without repetition
		for (int j = i + 1; j < i_high; j++) {
			Vector x_j = { particles[j].x, particles[j].y, particles[j].z };

			// Calculate distance between ith and jth particle
			double dist = (x_i - x_j).len();

			// Add spring if particles are close enough
			if (dist < max_dist) {
				add_spring(n_body_sim, i, j, spr); // Will not be added if there is already a spring present
			}
		}
	}
}

// Set the spring damping coefficient for all springs
void set_gamma(double new_gamma) {
	for (int i = 0; i < num_springs; i++) {
		springs[i].gamma = new_gamma;
	}
	std::cout.precision(2);
	std::cout << "\n gamma set to " << new_gamma << std::endl;
}

// Adjust spring constant, damping coefficient, and heat diffusion coefficient for springs with midpoints between r_min and r_max from center of mass
void adjust_spring_props(struct reb_simulation *const n_body_sim, double new_k,
		double new_gamma, double new_k_heat, double r_min, double r_max) {
	// Get particle and index info
	struct reb_particle *particles = n_body_sim->particles;

	// Search all particles
	int i_low = 0;
	int i_high = n_body_sim->N - num_perts;

	// Computer center of mass of all particles
	Vector CoM = compute_com(n_body_sim, i_low, i_high);

	// Count number of springs modified
	int NC = 0;
	for (int i = 0; i < num_springs; i++) {
		// Compute spring mid point from central position
		Vector x_mid = spr_mid(n_body_sim, springs[i], CoM);
		double r_mid = x_mid.len();

		// Modify spring if within radius bounds
		if ((r_mid >= r_min) && (r_mid <= r_max)) {
			springs[i].k = new_k;
			springs[i].gamma = new_gamma;
			springs[i].k_heat = new_k_heat;
			NC++;
		}
	}
	std::cout.precision(3);
	std::cout << "adjust_spring_props: \n\t Number of springs changed: " << NC
			<< "\n\t Fraction of springs changed: "
			<< (double) NC / (double) num_springs << std::endl;
}

// Kill springs that have failed
void kill_springs(struct reb_simulation *const n_body_sim) {
	// Set particle index boundaries
	int i_low = 0;
	int i_high = n_body_sim->N - num_perts;

	// Search particles for failing springs
	for (int i = i_low; i < i_high; i++) {
		if (stresses[i].failing) {
			int spring_index = stresses[i].max_force_spring; // spring with max force amplitude on it
			springs[spring_index].k = 0.0;   // spring dead
			springs[spring_index].gamma = 0.0;
		}
	}
}

// Make a binary with two masses mass_1, mass_2 spinning with vector omega, separated by distance sep and connected with spring spring_vals
// Center of mass is set to origin
void make_binary_spring(struct reb_simulation *const n_body_sim, double mass_1,
		double mass_2, double sep, Vector omega, spring spring_vals) {
	// Get particle info
	int i_low = n_body_sim->N;
	int i_high = i_low + 2;

	// Declare new particle
	struct reb_particle particle;

	// Set particle defaults
	particle.ax = 0.0;
	particle.ay = 0.0;
	particle.az = 0.0;
	particle.vx = 0.0;
	particle.vy = 0.0;
	particle.vz = 0.0;
	particle.y = 0.0;
	particle.z = 0.0;

	// Add particle 1
	particle.m = mass_1;
	particle.x = sep * mass_2 / (mass_1 + mass_2);
	particle.r = sep * 0.3; // How was this determined????
	reb_add(n_body_sim, particle);

	// Add particle 2
	particle.m = mass_2;
	particle.x = -sep * mass_1 / (mass_1 + mass_2);
	particle.r *= pow(mass_1 / mass_2, 0.33333); // How was this determined?????
	reb_add(n_body_sim, particle);

	// Spin the two particles
	spin_body(n_body_sim, i_low, i_high, omega);

	// Connect particles with a spring
	connect_springs_dist(n_body_sim, sep * 1.1, i_low, i_high, spring_vals);
}

/*********************/
/* Spring properties */
/*********************/

// Return spring length vector
Vector spring_r(struct reb_simulation *const n_body_sim, spring spr) {
	// Get simulation and spring information
	struct reb_particle *particles = n_body_sim->particles;
	int i = spr.particle_1;
	int j = spr.particle_2;

	// Get locations of endpoints of spring
	Vector x_i = { particles[i].x, particles[i].y, particles[i].z };
	Vector x_j = { particles[j].x, particles[j].y, particles[j].z };

	// Return length of spring
	return (x_i - x_j);
}

// Compute mean rest length of springs
double mean_spring_length() {
	// Initialize total length
	double total_length = 0.0;

	// Sum rest lengths of springs
	for (int i = 0; i < num_springs; i++) {
		total_length += springs[i].rs0;
	}

	// Return average spring length
	return total_length / num_springs;
}

// Return mean spring constant of all springs
double mean_k(struct reb_simulation *const n_body_sim) {
	// Init vars
	double k_tot = 0.0;

	// Sum spring constants
	for (int i = 0; i < num_springs; i++) {
		k_tot += springs[i].k;
	}

	// Return average
	return k_tot / num_springs;
}

// Compute spring midpoint location from arbitrary center (Cartesian coordinates)
Vector spr_mid(struct reb_simulation *const n_body_sim, spring spr,
		Vector center) {
	// Get simulation and spring information
	struct reb_particle *particles = n_body_sim->particles;
	int i = spr.particle_1;
	int j = spr.particle_2;

	// Get locations of endpoints of spring
	Vector x_i = { particles[i].x, particles[i].y, particles[i].z };
	Vector x_j = { particles[j].x, particles[j].y, particles[j].z };

	// Calculate location of midpoint of spring from arbitrary center position
	return 0.5 * (x_i + x_j) - center;
}

// Return spring strain
double strain(struct reb_simulation *const n_body_sim, spring spr) {
	double len = spring_r(n_body_sim, spr).len(); // Spring length
	return (len - spr.rs0) / spr.rs0; // Positive under extension
}

// Calculate spring forces using viscoelastic model
// Clavet et al. '05, "Particle-based Viscoelastic Fluid Simulation"
// Eurographics/ACM SIGGRAPH Symposium on Computer Animation
void spring_forces(struct reb_simulation *const n_body_sim) {
	// Get particle info
	reb_particle *particles = n_body_sim->particles;

	// Get spring forces
	for (int spring_index = 0; spring_index < num_springs; spring_index++) {
		// Get spring lengths
		double len = spring_r(n_body_sim, springs[spring_index]).len() + L_EPS
		;
		double len0 = springs[spring_index].rs0;
		double k = springs[spring_index].k;

		// Get info for spring endpoint particles
		int i = springs[spring_index].particle_1;
		int j = springs[spring_index].particle_2;
		Vector x_i = { particles[i].x, particles[i].y, particles[i].z };
		Vector x_j = { particles[j].x, particles[j].y, particles[j].z };
		Vector dx = (x_i - x_j);
		Vector len_hat = dx / len;
		double m_i = particles[i].m;
		double m_j = particles[j].m;

		// Calculate force from spring
		double spring_force = -k * (len - len0);

		// Apply elastic accelerations
		particles[i].ax += (spring_force * len_hat / m_i).getX();
		particles[j].ax -= (spring_force * len_hat / m_j).getX();
		particles[i].ay += (spring_force * len_hat / m_i).getY();
		particles[j].ay -= (spring_force * len_hat / m_i).getY();
		particles[i].az += (spring_force * len_hat / m_i).getZ();
		particles[j].az -= (spring_force * len_hat / m_i).getZ();

		// Apply damping, dependent on strain rate
		double gamma = springs[spring_index].gamma;
		if (gamma > 0.0) {
			// Get speed of spread of endpoints
			Vector v_i = { particles[i].vx, particles[i].vy, particles[i].vz };
			Vector v_j = { particles[j].vx, particles[j].vy, particles[j].vz };
			Vector dv = v_i - v_j;

			// Strain rate
			double strain_rate = dot(dv, len_hat) / len;

			// Reduced mass
			double m_red = m_i * m_j / (m_i + m_j);

			// Apply damping force
			Vector damp_force = gamma * m_red * strain_rate * len_hat;
			particles[i].ax -= (damp_force / m_i).getX();
			particles[j].ax += (damp_force / m_j).getX();
			particles[i].ay -= (damp_force / m_i).getY();
			particles[j].ay += (damp_force / m_j).getY();
			particles[i].az -= (damp_force / m_i).getZ();
			particles[j].az += (damp_force / m_j).getZ();
		}
	}
}

// Return the spring force vector for spring spring_index
// Divide by mass to get acceleration of particle i or j
Vector spring_i_force(struct reb_simulation *const n_body_sim,
		int spring_index) {
	// Get particle info
	reb_particle *particles = n_body_sim->particles;

	// Get spring lengths
	double len = spring_r(n_body_sim, springs[spring_index]).len() + L_EPS
	;
	double len0 = springs[spring_index].rs0;
	double k = springs[spring_index].k;

	// Get info for spring endpoint particles
	int i = springs[spring_index].particle_1;
	int j = springs[spring_index].particle_2;
	Vector x_i = { particles[i].x, particles[i].y, particles[i].z };
	Vector x_j = { particles[j].x, particles[j].y, particles[j].z };
	Vector dx = (x_i - x_j);
	Vector len_hat = dx / len;
	double m_i = particles[i].m;
	double m_j = particles[j].m;

	// Calculate force from spring
	Vector F = -k * (len - len0) * len_hat;

	// Apply damping, dependent on strain rate
	double gamma = springs[spring_index].gamma;
	if (gamma > 0.0) {
		// Get speed of spread of endpoints
		Vector v_i = { particles[i].vx, particles[i].vy, particles[i].vz };
		Vector v_j = { particles[j].vx, particles[j].vy, particles[j].vz };
		Vector dv = v_i - v_j;

		// Strain rate
		double strain_rate = dot(dv, len_hat) / len;

		// Reduced mass
		double m_red = m_i * m_j / (m_i + m_j);

		// Apply damping force
		Vector damp_force = gamma * m_red * strain_rate * len_hat;
		F -= damp_force;
	}

	return F;
}

// Return the spring force vector for spring i, ignoring damping
Vector spring_i_force_undamped(struct reb_simulation *const n_body_sim,
		int spring_index) {
	// Get particle info
	reb_particle *particles = n_body_sim->particles;

	// Get spring lengths
	double len = spring_r(n_body_sim, springs[spring_index]).len() + L_EPS
	;
	double len0 = springs[spring_index].rs0;
	double k = springs[spring_index].k;

	// Get info for spring endpoint particles
	int i = springs[spring_index].particle_1;
	int j = springs[spring_index].particle_2;
	Vector x_i = { particles[i].x, particles[i].y, particles[i].z };
	Vector x_j = { particles[j].x, particles[j].y, particles[j].z };
	Vector dx = (x_i - x_j);
	Vector len_hat = dx / len;
	double m_i = particles[i].m;
	double m_j = particles[j].m;

	// Calculate force from spring
	return -k * (len - len0) * len_hat;
}
