#ifdef __cplusplus
# 	ifdef __GNUC__
#		define restrict __restrict__
#	else
#		define restrict
#	endif
#endif

#include <cmath>
#include <vector>
#include <iostream>
extern "C" {
#include "rebound.h"
}
#include "matrix_math.h"
#include "physics.h"
#include "springs.h"
#include "stress.h"

using std::vector;

extern vector<spring> springs;
extern int num_springs; // Number of springs
extern int num_perts;

vector<stress_tensor> stresses; // Stress at each node (particle)

// Convention is that tensile stress is negative
// Update the stress tensor
// Create it if it does not exist
// Caution: Assumes total volume = 4/3 pi is a constant (used to normalize)
void update_stress(reb_simulation *const n_body_sim) {
	// Allocate stress if unallocated
	if (!stresses.size()) {
		stresses.resize(n_body_sim->N);
	}

	// Initialize stresses for each node
#pragma omp parallel for
	for (int i = 0; i < n_body_sim->N; i++) {
		stresses[i].stress = zero_matrix;
		stresses[i].eigs[0] = 0.0;
		stresses[i].eigs[1] = 0.0;
		stresses[i].eigs[2] = 0.0;
		stresses[i].max_force = 0.0;
		stresses[i].max_force_spring = -1;
		stresses[i].failing = false;
	}

	// Calculate stress from each spring
#pragma omp parallel for
	for (int k = 0; k < num_springs; k++) {
		// Get particles
		int i = springs[k].particle_1;
		int j = springs[k].particle_2;

		// Calculate force and length vector of spring k
		Vector F = spring_i_force(n_body_sim, k);
		Vector L = spring_r(n_body_sim, springs[k]);

		// Add stresses to each node
		// Tensile stress should be negative
		stresses[i].stress += outer(F, L);
		stresses[j].stress -= outer(F, L);

		// Keep track of which spring is applying maximum force
		double F_mag = F.len();
		if (F_mag > stresses[i].max_force) {
			stresses[i].max_force = F_mag;
			stresses[i].max_force_spring = k;
		}
		if (F_mag > stresses[j].max_force) {
			stresses[j].max_force = F_mag;
			stresses[j].max_force_spring = k;
		}
	}

	// Normalize the stress
	// We set R = 1 for our masses, so the total volume is 4/3 pi
	// Volume of a single node depends on number of particles
	double vol = 4.0 * M_PI / (3.0 * (double) n_body_sim->N);
#pragma omp parallel for
	for (int i = 0; i < n_body_sim->N; i++) {
		stresses[i].stress /= vol;
	}

	// Compute and store eigenvalues for each node
#pragma omp parallel for
	for (int i = 0; i < n_body_sim->N; i++) {
		try {
			eigenvalues(stresses[i].stress, stresses[i].eigs);
		} catch (char *str) {
			std::cerr << "Error: " << str << " Exiting." << std::endl;
			exit(1);
		}
	}
}

// Check for failure at all nodes
// tens_str_int is tensile strength of interior (these are negative)
// tens_str_surf is tensile strength of shell
// Note: we need to know if the particle was originally in shell or interior
// Decide if particle is in interior by whether particle mass is below pmass_div
// Caution: we can't use current radius as particles are going to move!
// Return number of failed nodes
int mark_failed_nodes(reb_simulation *const n_body_sim, double mass_div,
		double tens_str_int, double tens_str_surf) {
	// Get particle indices
	int i_low = 0;
	int i_high = n_body_sim->N - num_perts;

	// Init number of failures
	int n_fail = 0;

	// Determine if node fails
#pragma omp parallel for
	for (int i = i_low; i < i_high; i++) {
		// Get particle mass and tensile strength
		double m = n_body_sim->particles[i].m;
		double tens_str;

		// Particle is interior if mass is above pmass_div
		if (m > mass_div) {
			tens_str = tens_str_int;
		} else {
			tens_str = tens_str_surf;
		}

		// Compare stress eigenvalues to failure conditions
		double tau1 = stresses[i].eigs[0];
		double tau3 = stresses[i].eigs[2];
		double tau13 = tau1 - tau3;
		if (((tau1 < -3.0 * tau3) && (tau3 < tens_str))
				|| ((tau1 >= -3.0 * tau3)
						&& (tau13 * tau13 + 8.0 * tens_str * tau13 > 0))) {
			stresses[i].failing = true;
			n_fail++;
		}
	}

	// Return number of failed nodes
	if (n_fail > 0)
		std::cout << "markfailure: " << n_fail << " nodes failed." << std::endl;

	return n_fail;
}

// Compute Young's modulus of springs using midpoints of springs from center of mass in radial range [r_min,r_max]
// Equation 20 by Kot et al. 2014 -> sum_i [k_i L_i^2]/(6V)
// Uses rest lengths
// Only computes center of mass using particles in range [i_low,i_high)
double Young_mesh(reb_simulation *const n_body_sim, int i_low, int i_high,
		double r_min, double r_max) {
	// Initialize sum
	double sum = 0.0;

	// Find center of mass of requested particles
	Vector CoM = compute_com(n_body_sim, i_low, i_high);

	// Calculate sum
#pragma omp parallel for
	for (int i = 0; i < num_springs; i++) {
		Vector x_mid = spr_mid(n_body_sim, springs[i], CoM);

		// If center of spring is within requested shell, update sum
		double dist_from_CoM = x_mid.len();
		if ((dist_from_CoM < r_max) && (dist_from_CoM > r_min)) {
			double k = springs[i].k;
			double len = springs[i].rs0;
			sum += k * len * len;
		}
	}

	// Calculate volume of shell
	double volume = (4.0 * M_PI / 3.0) * (pow(r_max, 3.0) - pow(r_min, 3.0));

	// Return Young's modulus
	return sum / (6.0 * volume);
}

// Calculate Young's modulus of entire simulation
// Equation 20 by Kot et al. 2014 -> sum_i [k_i L_i^2]/(6V)
// Uses rest lengths
double Young_full_mesh() {
	// Initialize sum
	double sum = 0.0;

	// Update sum for each spring
#pragma omp parallel for
	for (int i = 0; i < num_springs; i++) {
		double k = springs[i].k;
		double len = springs[i].rs0;
		sum += k * len * len;
	}

	// Assume sphere of radius 1
	double volume = 4.0 * M_PI / 3.0;

	// Return Young's modulus
	return sum / (6.0 * volume);
}

/*************/
/* Operators */
/*************/

// Stream output
std::ostream& operator<<(std::ostream &os, const stress_tensor &stress) {
	os << "Stress tensor: " << stress.stress << "\nEigenvalues: "
			<< stress.eigs[0] << ", " << stress.eigs[1] << ", "
			<< stress.eigs[2] << "\nMax force: " << stress.max_force
			<< "\nSpring with max force: " << stress.max_force_spring
			<< "\nFailing?: " << stress.failing;
	return os;
}
