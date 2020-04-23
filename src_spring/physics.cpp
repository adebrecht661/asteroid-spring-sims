#ifdef __cplusplus
# 	ifdef __GNUC__
#		define restrict __restrict__
#	else
#		define restrict
#	endif
#endif

/*
 * physics.cpp
 *
 * Various physics routines
 *
 *  Created on: Apr 2, 2020
 *      Author: alex
 */

#include <cmath>
#include <iostream>
extern "C" {
#include "rebound.h"
}
#include "matrix_math.h"
#include "springs.h"
#include "physics.h"

extern int num_springs;
extern spring springs[];
extern int num_perts;
extern const double L_EPS;

/*******************/
/* Center routines */
/*******************/

// Compute center of mass coordinates in particle range [i_low,i_high)
Vector compute_com(reb_simulation *const n_body_sim, int i_low,
		int i_high) {
	// Get simulation information and initialize sums
	reb_particle *particles = n_body_sim->particles;
	double m_tot = 0.0;
	Vector x_times_m(zero_vec);

	// Sum masses and mass-weighted locations
	for (int i = i_low; i < i_high; i++) {
		x_times_m += particles[i].m * Vector( { particles[i].x, particles[i].y,
				particles[i].z });
		m_tot += particles[i].m;
	}

	// Calculate center of mass
	return x_times_m / m_tot;
}

// Compute center of velocity of particles in particle range [i_low,i_high)
Vector compute_cov(reb_simulation *const n_body_sim, int i_low,
		int i_high) {
	// Get simulation information and initialize sums
	reb_particle *particles = n_body_sim->particles;
	Vector v_times_m(zero_vec);
	double m_tot = 0.0;

	// Sum masses and mass-weighted velocities
	for (int i = i_low; i < i_high; i++) {
		v_times_m += particles[i].m * Vector( { particles[i].vx,
				particles[i].vy, particles[i].vz });
		m_tot += particles[i].m;
	}

	return v_times_m / m_tot;
}

// Recenter resolved body on its center of mass
// Caution: only affects particles in set [i_low, i_high)
void subtract_com(reb_simulation *const n_body_sim, int i_low, int i_high) {

	// Find center of mass of resolved body
	Vector CoM = compute_com(n_body_sim, i_low, i_high);

	// Move body
	move_resolved(n_body_sim, -CoM, zero_vector, i_low, i_high);
}

// Recenter frame of resolved body on its center of velocity
// Caution: only affects particles in set [i_low, i_high)
void subtract_cov(reb_simulation *const n_body_sim, int i_low,
		int i_high) {
	// Find center of velocity of resolved body
	Vector CoV = compute_cov(n_body_sim, i_low, i_high);

	// Move body
	move_resolved(n_body_sim, zero_vector, -CoV, i_low, i_high);
}

// Move simulation to coordinate frame of body defined by particles [i_low,i_high)
// Shifts all particles
void center_sim(reb_simulation *const n_body_sim, int i_low,
		int i_high) {
	// Find CoM of request particles
	Vector CoM = compute_com(n_body_sim, i_low, i_high);

	// Shift coordinates of all particles
	for (int i = 0; i < (n_body_sim->N); i++) {
		n_body_sim->particles[i].x -= CoM.getX();
		n_body_sim->particles[i].y -= CoM.getY();
		n_body_sim->particles[i].z -= CoM.getZ();
	}
}

/***********************/
/* Rotational routines */
/***********************/

// Compute the moment of inertia tensor of particles in range [i_low,i_high) with respect to center of mass
Matrix mom_inertia(reb_simulation *const n_body_sim, int i_low,
		int i_high) {
	// Get particle information and initialize inertia matrix
	reb_particle *particles = n_body_sim->particles;
	Matrix inertia(zero_mat);

	// Calculate center of mass of specified body
	Vector CoM = compute_com(n_body_sim, i_low, i_high);

	// Calculate moment of inertia of particles
	for (int i = i_low; i < i_high; i++) {
		Vector dx = { particles[i].x - CoM.getX(), particles[i].y - CoM.getY(),
				particles[i].z - CoM.getZ() };
		inertia += particles[i].m * pow(dx.len(), 2.0) * I;
		inertia -= particles[i].m * outer(dx, dx);
	}

	// Return inertia tensor
	return inertia;
}

// Compute (spin) angular momentum vector of particles in range [i_low, i_high) with respect to their center of mass position and velocity
// Caution: Can measure angular momentum of the entire system, but will be with respect to center of mass of entire system
Vector measure_L(reb_simulation *const n_body_sim, int i_low,
		int i_high) {
	// Get particle info
	reb_particle *particles = n_body_sim->particles;

	// Find center of mass and center of velocity
	Vector CoM = compute_com(n_body_sim, i_low, i_high);
	Vector CoV = compute_cov(n_body_sim, i_low, i_high);

	// Initialize angular momentum
	Vector L(zero_vec);

	// Calculate angular momentum of each particle with respect to centers of mass and velocity
	for (int i = i_low; i < i_high; i++) {
		Vector dx = { particles[i].x, particles[i].y, particles[i].z };
		dx -= CoM;
		Vector dv = { particles[i].vx, particles[i].vy, particles[i].vz };
		dv -= CoV;
		Vector dL = cross(dx, dv);
		L += particles[i].m * dL; // angular momentum vector
	}

	// Return total angular momentum
	return L;

}

// Compute spin vector of body with particles in range [i_low,i_high) using inverse of moment of inertia matrix
// Also returns eigenvalues of moment of inertia matrix, eigs[0] >= eigs[1] >= eigs[2]
Vector body_spin(reb_simulation *const n_body_sim, int i_low, int i_high,
		double eigs[3]) {

	// Compute moment of inertia matrix
	Matrix inertia_mat = mom_inertia(n_body_sim, i_low, i_high);

	// Compute inverse of moment of inertia matrix
	Matrix inv_mom_inert = inverse(inertia_mat);

	// Compute angular momentum vector with respect to its center of mass position and velocity
	Vector L = measure_L(n_body_sim, i_low, i_high);

	// Compute eigenvalues and spin vector
	eigenvalues(inertia_mat, eigs);
	return inv_mom_inert * L;
}

// Start particles in range [i_low,i_high) spinning with vector omega about their center of mass
void spin_body(reb_simulation *const n_body_sim, int i_low, int i_high,
		Vector omega) {
	// Get particle info
	reb_particle *particles = n_body_sim->particles;

	// Compute center of mass and velocity
	Vector CoM = compute_com(n_body_sim, i_low, i_high);
	Vector CoV = compute_cov(n_body_sim, i_low, i_high);

	// Operate only if spin is large enough
	if (omega.len() > 1e-5) {
		// Add necessary velocity to each particle
		for (int i = i_low; i < i_high; i++) {

			Vector dx = { particles[i].x, particles[i].y, particles[i].z };
			dx -= CoM;
			Vector r_cross_omega = cross(dx, omega);

			// Set it spinning with respect to center of mass
			particles[i].vx = CoV.getX() + r_cross_omega.getX();
			particles[i].vy = CoV.getY() + r_cross_omega.getY();
			particles[i].vz = CoV.getZ() + r_cross_omega.getZ();
		}
	}

}

// Compute the orbital angular momentum vector of particles in range [i_low,i_high)
// About central mass if number of perturbers is 1, otherwise about the center of mass of all the perturbers
Vector compute_Lorb(reb_simulation *const n_body_sim, int i_low,
		int i_high) {
	reb_particle *particles = n_body_sim->particles;
	static int first = 0;
	static double tm = 0.0;  // its mass

	if (first == 0) { // only calculate once
		for (int i = i_low; i < i_high; i++)
			tm += n_body_sim->particles[i].m; // mass of resolved body
		first = 1;
	}

	// Find center of mass and center of velocity of requested particles
	Vector CoM = compute_com(n_body_sim, i_low, i_high);
	Vector CoV = compute_cov(n_body_sim, i_low, i_high);

	// Find center of mass and center of velocity of perturbing particles
	Vector CoM0;
	Vector CoV0;
	if (num_perts == 1) {
		// Get index, location and velocity of primary perturber
		int im1 = n_body_sim->N - 1;
		CoM0 = { particles[im1].x, particles[im1].y, particles[im1].z };
		CoV0 = { particles[im1].vx, particles[im1].vy, particles[im1].vz };
	} else {
		// Get index range for perturbing masses
		int iml = n_body_sim->N - num_perts;
		int imh = n_body_sim->N;

		// Get center of mass and velocity of perturbing masses
		CoM0 = compute_com(n_body_sim, iml, imh);
		CoV0 = compute_cov(n_body_sim, iml, imh);
	}

	// Calculate and return orbital angular momentum
	Vector dx = CoM0 - CoM;
	Vector dv = CoV0 - CoV;
	return cross(dx, dv);

}

// Compute non-translational kinetic energy of particles in range [i_low,i_high)
// Caution: Includes kinetic energy in vibration as well as rotation
double compute_rot_kin(reb_simulation *const n_body_sim, int i_low,
		int i_high) {
	// Get particle info
	reb_particle *particles = n_body_sim->particles;

	// Find center of velocity of requested particles
	Vector CoV = compute_cov(n_body_sim, i_low, i_high);

	// Sum all non-body kinetic energy (rotational, vibrational)
	double KE = 0.0;
	for (int i = i_low; i < i_high; i++) {
		// Get velocity difference from bulk
		Vector v = { particles[i].vx, particles[i].vy, particles[i].vz };
		v -= CoV;

		// Add kinetic energy
		KE += 0.5 * particles[i].m * pow(v.len(), 2.0);
	}
	return KE;
}

// Compute (spin) angular momentum vector of all particles in range [i_low, i_high) with respect to origin
Vector measure_L_origin(reb_simulation *const n_body_sim, int i_low,
		int i_high) {
	// Get particle info
	reb_particle *particles = n_body_sim->particles;

	// Calculate angular momentum
	Vector L;
	for (int i = i_low; i < i_high; i++) {
		Vector dx = { particles[i].x, particles[i].y, particles[i].z };
		Vector dv = { particles[i].vx, particles[i].vy, particles[i].vz };
		L += particles[i].m * cross(dx, dv);
	}

	// Return angular momentum vector
	return L;
}

// Rotate a body with particle indices [i_low, i_high) about center of mass using Euler angles
// Rotates both position and velocity
// Center of mass position and velocity is not changed
void rotate_body(reb_simulation *const n_body_sim, int i_low, int i_high,
		double alpha, double beta, double gamma) {
	reb_particle *particles = n_body_sim->particles;

	// Compute center of mass and velocity
	Vector CoM = compute_com(n_body_sim, i_low, i_high);
	Vector CoV = compute_cov(n_body_sim, i_low, i_high);

	// Set rotation matrices
	Matrix Rz1 = getRotMatZ(alpha); // Rotate about z
	Matrix Rx = getRotMatX(beta); // Rotate about x'
	Matrix Rz2 = getRotMatZ(gamma); // Rotate about z''

	// Calculate new position and velocity of each particle
	for (int i = i_low; i < i_high; i++) {

		// Get position relative to CoM
		Vector x0 = { particles[i].x, particles[i].y, particles[i].z };
		x0 -= CoM;

		// Rotate
		Vector x_rot = Rz2 * Rx * Rz1 * x0;

		// Update relative to CoM
		particles[i].x = x_rot.getX() + CoM.getX();
		particles[i].y = x_rot.getY() + CoM.getY();
		particles[i].z = x_rot.getZ() + CoM.getZ();

		// Get velocity relative to CoV
		Vector v0 = { particles[i].vx, particles[i].vy, particles[i].vz };
		v0 -= CoV;

		// Rotate
		Vector v_rot = Rz2 * Rx * Rz1 * v0;

		// Update relative to CoV
		particles[i].vx = v_rot.getX() + CoV.getX();
		particles[i].vy = v_rot.getY() + CoV.getY();
		particles[i].vz = v_rot.getZ() + CoV.getZ();
	}
}

// Rotate all particles in set [i_low, i_high) about origin
// Rotates both position and velocity
void rotate_origin(reb_simulation *const n_body_sim, int i_low,
		int i_high, double alpha, double beta, double gamma) {
	// Get particle info
	reb_particle *particles = n_body_sim->particles;

	// Set rotation matrices
	Matrix Rz1 = getRotMatZ(alpha); // Rotate about z
	Matrix Rx = getRotMatX(beta); // Rotate about x'
	Matrix Rz2 = getRotMatZ(gamma); // Rotate about z''

	// Calculate new position and velocity of each particle
	for (int i = i_low; i < i_high; i++) {

		// Get current position
		Vector x0 = { particles[i].x, particles[i].y, particles[i].z };

		// Rotate
		Vector x_rot = Rz2 * Rx * Rz1 * x0;

		// Update particle location
		particles[i].x = x_rot.getX();
		particles[i].y = x_rot.getY();
		particles[i].z = x_rot.getZ();

		// Get current velocity
		Vector v0 = { particles[i].vx, particles[i].vy, particles[i].vz };

		// Rotate
		Vector v_rot = Rz2 * Rx * Rz1 * v0;

		// Update particle velocity
		particles[i].vx = v_rot.getX();
		particles[i].vy = v_rot.getY();
		particles[i].vz = v_rot.getZ();
	}
}

// Rotate particles so that z is along biggest eigenvalue, x along smallest
void rotate_to_principal(reb_simulation *const n_body_sim, int i_low,
		int i_high) {
	// Get particle info
	reb_particle *particles = n_body_sim->particles;

	// Calculate moment of inertia matrix
	Matrix inertia_mat = mom_inertia(n_body_sim, i_low, i_high);
	std::cout.precision(3);
	std::cout << "rotate_to_principal initial moment of inertia: "
			<< inertia_mat << "\n";

	// Get eigenvalues of moment of inertia
	double eigs[3];
	eigenvalues(inertia_mat, eigs);
	std::cout << "Moment of inertia eigenvalues: " << eigs[0] << ", " << eigs[1]
			<< ", " << eigs[2] << "\n";

	// Get eigenvectors of moment of inertia matrix
	// Note order
	Vector eigvec1 = eigenvector(inertia_mat, eigs[2]);
	std::cout << "MoI eigvec 1: << " << eigvec1 << "\n";
	Vector eigvec2 = eigenvector(inertia_mat, eigs[1]);
	std::cout << "MoI eigvec 2: << " << eigvec2 << "\n";
	Vector eigvec3 = eigenvector(inertia_mat, eigs[0]);
	std::cout << "MoI eigvec 3: << " << eigvec3 << "\n";

	// Get eigenvector matrix for rotations
	Matrix eigvec_mat = { eigvec1, eigvec2, eigvec3 };

	// Rotate using eigenvector matrix
	for (int i = i_low; i < i_high; i++) {

		// Get current position
		Vector x0 = { particles[i].x, particles[i].y, particles[i].z };

		// Rotate
		Vector x_rot = eigvec_mat * x0;

		// Update position
		particles[i].x = x_rot.getX();
		particles[i].y = x_rot.getY();
		particles[i].z = x_rot.getZ();

		// Get current velocity
		Vector v0 = { particles[i].vx, particles[i].vy, particles[i].vz };

		// Rotate
		Vector v_rot = eigvec_mat * v0;

		// Update velocities
		particles[i].vx = v_rot.getX();
		particles[i].vy = v_rot.getY();
		particles[i].vz = v_rot.getZ();
	}

	// Check if rotation worked
	inertia_mat = mom_inertia(n_body_sim, i_low, i_high);
	std::cout << "rotate_to_principal final moment of inertia: " << inertia_mat
			<< "\n";
	eigenvalues(inertia_mat, eigs);
	std::cout << "Moment of inertia eigenvalues: " << eigs[0] << ", " << eigs[1]
			<< ", " << eigs[2] << std::endl;
}

/*******************/
/* Linear routines */
/*******************/

// Sum total momentum of particles
Vector total_mom(reb_simulation *const n_body_sim, int i_low,
		int i_high) {
	// Get particle info, initialize momentum vector
	reb_particle *particles = n_body_sim->particles;
	Vector p(zero_vec);

	// Sum momentum
	for (int i = i_low; i < i_high; i++) {
		p += particles[i].m * Vector( { particles[i].vx, particles[i].vy,
				particles[i].vz });
	}

	// Return momentum
	return p;
}

// Shift the position and velocity of a resolved body with particles in set [i_low, i_high) by dx, dv
void move_resolved(reb_simulation *const n_body_sim, Vector dx,
		Vector dv, int i_low, int i_high) {
	// Get particle info
	reb_particle *particles = n_body_sim->particles;

	// Update requested particles
	for (int i = i_low; i < i_high; i++) {
		particles[i].x += dx.getX();
		particles[i].y += dx.getY();
		particles[i].z += dx.getZ();
		particles[i].vx += dv.getX();
		particles[i].vy += dv.getY();
		particles[i].vz += dv.getZ();
	}
}

/********************/
/* Energy and power */
/********************/

// Compute power lost in damping from a specific spring
double dEdt(reb_simulation *const n_body_sim, spring spr) {
	// Get damping coefficient and particle info
	reb_particle *particles = n_body_sim->particles;
	double gamma = spr.gamma;

	// If no damping, return immediately
	if (gamma == 0.0)
		return 0.0;

	// Get spring length vector
	Vector dx = spring_r(n_body_sim, spr);
	double len = dx.len() + L_EPS;
	Vector len_hat = dx / len;

	// Get particle masses
	int i = spr.particle_1;
	int j = spr.particle_2;
	double m_i = particles[i].m;
	double m_j = particles[j].m;

	// Get speed of spread of endpoints
	Vector v_i = { particles[i].vx, particles[i].vy, particles[i].vz };
	Vector v_j = { particles[j].vx, particles[j].vy, particles[j].vz };
	Vector dv = v_i - v_j;

	// Strain rate
	double strain_rate = dot(dv, len_hat);

	// Reduced mass
	double m_red = m_i * m_j / (m_i + m_j);

	// Return power (dE/dt)
	return gamma * m_red * pow(strain_rate, 2.0);
}

// Compute power lost to damping over all springs
double dEdt_total(reb_simulation *const n_body_sim) {
	// Initialize sum
	double power_loss = 0.0;

	// Sum power lost for each spring
	for (int i = 0; i < num_springs; i++) {
		power_loss += dEdt(n_body_sim, springs[i]);
	}

	// Return total power lost
	return power_loss;
}

// Compute total potential energy in spring network
double spring_potential_energy(reb_simulation *const n_body_sim) {
	// Init total potential energy
	double SPE = 0.0;

	// Calculate potential energy in each spring
	for (int i = 0; i < num_springs; i++) {
		// Get length and spring constant
		double len0 = springs[i].rs0;
		double len = spring_r(n_body_sim, springs[i]).len() + L_EPS;
		double k = springs[i].k;

		// Add to potential energy
		SPE += k * pow(len - len0, 2.0) / 2.0;
	}

	// Return total potential energy
	return SPE;
}

// Compute total gravitational potential energy between particles in range [i_low, i_high)
double grav_potential_energy(reb_simulation *const n_body_sim, int i_low,
		int i_high) {
	// Get particle info
	reb_particle *particles = n_body_sim->particles;

	// Init total potential energy
	double GPE = 0.0;

	// Calculate potential energy between each pair of particles
	for (int i = i_low; i < i_high; i++) {
		for (int j = i + 1; j < i_high; j++) {
			// Get distance
			Vector x_1 = { particles[i].x, particles[i].y, particles[i].z };
			Vector x_2 = { particles[j].x, particles[j].y, particles[j].z };
			Vector dx = x_1 - x_2;

			GPE -= n_body_sim->G * particles[i].m * particles[j].m / dx.len();
		}
	}

	// Return total potential energy
	return GPE;
}

/********/
/* Misc */
/********/

// Zero all particle accelerations
void zero_accel(reb_simulation *const n_body_sim) {
	for (int i = 0; i < (n_body_sim->N); i++) {
		n_body_sim->particles[i].ax = 0.0;
		n_body_sim->particles[i].ay = 0.0;
		n_body_sim->particles[i].az = 0.0;
	}
}

// Multiply all masses to the right of x_min by factor m_fac, then rescale
void adjust_mass_side(reb_simulation *const n_body_sim, double m_fac,
		double x_min) {
	// Get particle info
	reb_particle *particles = n_body_sim->particles;
	int i_low = 0;
	int i_high = n_body_sim->N - num_perts;

	// If particle if right of x_min, adjust mass and count as modified
	double total_mass = 0.0;
	int num_mod = 0;
	for (int i = i_low; i < i_high; i++) {
		if (particles[i].x > x_min) {
			particles[i].m *= m_fac;
			num_mod++;
		}
	}

	// Renormalize particle masses
	for (int i = i_low; i < i_high; i++)
		total_mass += particles[i].m;
	for (int i = i_low; i < i_high; i++)
		particles[i].m /= total_mass;

	std::cout << "adjust_mass_side: Modified masses of " << num_mod
			<< " particles." << std::endl;
}
