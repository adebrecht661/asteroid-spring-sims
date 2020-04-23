#ifdef __cplusplus
# 	ifdef __GNUC__
#		define restrict __restrict__
#	else
#		define restrict
#	endif
#endif

#include <cmath>
extern "C" {
#include "rebound.h"
}
#include "matrix_math.h"
#include "springs.h"
#include "kepcart.h"
#include "physics.h"
#include "orb.h"

extern int num_perts;

// Add binary perturbers separated by distance sep in a circular orbit
// Center of mass of new particles set to be in orbit with given orbital elements
// r_prim is radius for display of primary
// Orbital elements refer to orbit of resolved body (particles [0,N)) about binary
// The resolved body is assumed to be alone and at origin
// Returns mean motion of resolved body around binary
double add_bin_kep(reb_simulation *const n_body_sim, double m_prim,
		double r_prim, double m_ratio, double sep, OrbitalElements orb_el) {

	// Initialize particle
	reb_particle pt;
	pt.ax = 0.0;
	pt.ay = 0.0;
	pt.az = 0.0;
	pt.vx = 0.0;
	pt.vy = 0.0;
	pt.vz = 0.0;
	pt.x = 0.0;
	pt.y = 0.0;
	pt.z = 0.0;

	// Get new particle masses
	double m_1 = m_prim;
	double m_2 = m_ratio * m_1;
	double mtot = m_1 + m_2;

	// Get particle locations relative to center of mass of new binary
	double x_1 = sep * m_2 / mtot;
	double x_2 = -sep * m_1 / mtot;

	// Get particle velocities relative to center of velocity of new binary (?????)
	double vc_bin = sqrt(n_body_sim->G * mtot / sep);
	double vy_1 = vc_bin * m_2 / mtot;
	double vy_2 = -vc_bin * m_1 / mtot;

	// Get display radius of particles
	double r_1 = r_prim;
	double r_2 = r_prim * pow(m_2 / m_1, 0.333333); // display radius
	if (r_2 < 1.0)
		r_2 = 1.0;

	// Get phase space state
	// +1 here for extended body
	double GM = n_body_sim->G * (mtot + 1.0);
	PhaseState state = kep_to_cart(GM, orb_el);

	// Add secondary
	pt.x = state.x.getX() + x_2;
	pt.y = state.x.getY();
	pt.z = state.x.getZ();
	pt.vx = state.v.getX();
	pt.vy = state.v.getY() + vy_2;
	pt.vz = state.v.getZ();
	pt.m = m_2;
	pt.r = r_2;
	reb_add(n_body_sim, pt);

	// Add primary
	pt.x = state.x.getX() + x_1;
	pt.y = state.x.getY();
	pt.z = state.x.getZ();
	pt.vx = state.v.getX();
	pt.vy = state.v.getY() + vy_1;
	pt.vz = state.v.getZ();
	pt.m = m_1;
	pt.r = r_1;
	reb_add(n_body_sim, pt);

	// Get orbital angular velocity (mean motion) of resolved body
	double a = orb_el.a;
	double omega_orb = sqrt(GM / a) / a;

	// Ratio of reduced mass to binary mass
	double muBratio = m_ratio / pow(1.0 + m_ratio, 2.0);

	// Correct mean motion by quadratic moment of binary
	double quadratic_moment = 1.0 + (3.0 / 8.0) * muBratio * pow(sep / a, 2.0);
	omega_orb *= quadratic_moment;

	// Return mean motion of resolved body around binary
	return omega_orb;
}

// Add a point mass perturber of mass m, display radius r
// If i_p >= 0, put point mass in orbit about another point mass with index i_p (i_low, i_high are not used)
// If i_p < 0
//    The extended mass has indices [i_low, i_high)
//    The point mass is added at origin
//    If extended body is rotating then it still rotates, but if orbit is tilted then obliquity will not be the same
// Returns mean motion of new orbit
double add_pt_mass_kep(reb_simulation *const n_body_sim, int i_low,
		int i_high, int i_p, double mass, double radius,
		OrbitalElements orb_el) {
	// Get particle info
	reb_particle *particles = n_body_sim->particles;

	// Initialize particle
	reb_particle pt;
	pt.m = mass;
	pt.r = radius;
	pt.ax = 0;
	pt.ay = 0;
	pt.az = 0;

	double m0, GM;
	PhaseState state;
	// If i_p < 0, move resolved body
	if (i_p < 0) {
		// Get mass of resolved body
		m0 = sum_mass(n_body_sim, i_low, i_high);

		// Zero the center of mass position of resolved body
		subtract_com(n_body_sim, i_low, i_high);

		// Zero the center of mass velocity of resolved body
		subtract_cov(n_body_sim, i_low, i_high);

		// Get phase space state corresponding to requested orbit
		GM = n_body_sim->G * (mass + m0);
		state = kep_to_cart(GM, orb_el);

		// New point mass is at origin
		pt.x = 0.0;
		pt.y = 0.0;
		pt.z = 0.0;
		pt.vx = 0.0;
		pt.vy = 0.0;
		pt.vz = 0.0;

		// Move resolved body to orbit given by orb_el
		move_resolved(n_body_sim, state.x, state.v, i_low, i_high);

		// Otherwise, new particle has motion w.r.t to particle at i_p
	} else {
		// Get mass of particle it orbits
		m0 = particles[i_p].m;

		// Get position and velocity of particle it orbits
		Vector x0 = { particles[i_p].x, particles[i_p].y, particles[i_p].z };
		Vector v0 = { particles[i_p].vx, particles[i_p].vy, particles[i_p].vz };

		// Get phase space state corresponding to requested orbit
		GM = n_body_sim->G * (mass + m0);
		state = kep_to_cart(GM, orb_el);

		// New particle orbits particle i_p with orbit given by orb_el
		pt.x = state.x.getX() + x0.getX();
		pt.y = state.x.getY() + x0.getY();
		pt.z = state.x.getZ() + x0.getZ();
		pt.vx = state.v.getX() + v0.getX();
		pt.vy = state.v.getY() + v0.getY();
		pt.vz = state.v.getZ() + v0.getZ();
	}

	// Add new particle
	reb_add(n_body_sim, pt);

	// Return mean motion of new orbit
	double a = orb_el.a;
	return sqrt(GM / a) / a;
}

// Drift the orbits of two point masses part_1 and part_2 using Beauge et al. 06's formulae
// (Beauge, Michtchenko, Ferraz-Mello 2006, MNRAS 365, 1160)
// Maintains center of mass of binary
// Note: If inv_tau_a > 0, they are drifting apart
void drift_bin(reb_simulation *const n_body_sim, double timestep,
		double inv_tau_a, double inv_tau_e, int part_1, int part_2) {
	// Get particle info
	reb_particle *particles = n_body_sim->particles;

	// Get masses
	double m1 = particles[part_1].m;
	double m2 = particles[part_2].m;

	// Get position and velocity vectors
	Vector x_2 =
			{ particles[part_2].x, particles[part_2].y, particles[part_2].z };
	Vector x_1 =
			{ particles[part_1].x, particles[part_1].y, particles[part_1].z };
	Vector x = x_2 - x_1;

	Vector v_2 = { particles[part_2].vx, particles[part_2].vy,
			particles[part_2].vz };
	Vector v_1 = { particles[part_1].vx, particles[part_1].vy,
			particles[part_1].vz };
	Vector v = v_2 - v_1;

	// Get centripetal velocity (theoretical)
	double GM = n_body_sim->G * (m1 + m2);
	double v_cent = sqrt(GM / x.len());

	// Get unit vector in direction of orbit (r x r x v = r x L)
	Vector r_cross_L = cross(x, cross(x, v));
	r_cross_L /= r_cross_L.len();

	// Compare actual velocity and centripetal velocity
	Vector v_cent_vec = v_cent * r_cross_L;
	Vector delta_v = v - v_cent_vec;

	// Compute changes in velocity
	Vector dv = timestep * (v * inv_tau_a / 2.0 + delta_v * inv_tau_e);

	// Update velocities in such a way as to conserve momentum of binary
	// Lower mass is affected much more than higher mass
	Vector dv_1 = m1 * dv / (m1 + m2);
	Vector dv_2 = m2 * dv / (m1 + m2);

	// Particle 1
	particles[part_1].vx -= dv_1.getX();
	particles[part_1].vy -= dv_1.getY();
	particles[part_1].vz -= dv_1.getZ();

	// Particle 2
	particles[part_2].vx -= dv_2.getX();
	particles[part_2].vy -= dv_2.getY();
	particles[part_2].vz -= dv_2.getZ();
}

// Drift orbits of a particle and a resolved body
// See drift_bin for description of inv_tau_a, inv_tau_e
void drift_resolved(reb_simulation *const n_body_sim, double timestep,
		double inv_tau_a, double inv_tau_e, int i_part, int i_low, int i_high) {
	// Get particle info
	reb_particle *particles = n_body_sim->particles;

	// Get mass of particle and resolved body
	double m1 = particles[i_part].m;
	double m2 = sum_mass(n_body_sim, i_low, i_high);

	// Get center of mass and velocity of resolved body
	Vector CoM = compute_com(n_body_sim, i_low, i_high);
	Vector CoV = compute_cov(n_body_sim, i_low, i_high);

	// Find phase space state of resolved body relative to particle
	Vector x_0 =
			{ particles[i_part].x, particles[i_part].y, particles[i_part].z };
	Vector x = CoM - x_0;
	Vector v_0 = { particles[i_part].vx, particles[i_part].vy,
			particles[i_part].vz };
	Vector v = CoV - v_0;

	// Centripetal velocity (theoretical)
	double GM = n_body_sim->G * (m1 + m2);
	double rad = x.len();
	double v_cent = sqrt(GM / rad);

	// Get unit vector in direction of orbit (r x r x v = r x L)
	Vector r_cross_L = cross(x, cross(x, v));
	r_cross_L /= r_cross_L.len();

	// Compare actual velocity and centripetal velocity
	Vector v_cent_vec = v_cent * r_cross_L;
	Vector delta_v = v - v_cent_vec;

	// Compute changes in velocity
	Vector dv = timestep * (v * inv_tau_a / 2.0 + delta_v * inv_tau_e);

	// Update velocities in such a way as to conserve momentum of binary
	// Lower mass is affected much more than higher mass
	Vector dv_1 = m1 * dv / (m1 + m2);
	Vector dv_2 = m2 * dv / (m1 + m2);

	// Particle
	particles[i_part].vx -= dv_2.getX();
	particles[i_part].vy -= dv_2.getY();
	particles[i_part].vz -= dv_2.getZ();

	// Resolved body -- update velocities only
	move_resolved(n_body_sim, zero_vector, dv_1, i_low, i_high);
}

// Apply quadrupole force of particle with index i_p onto all particles
// J2_p is unitless
// R_p is radius of planet
// i_p is the array index of the planet with the quadrupole moment
// See https://en.wikipedia.org/wiki/Geopotential_model
// Takes into account orientation of planet's north pole
void quadrupole_accel(reb_simulation *const n_body_sim, double J2_p,
		double R_p, double phi_p, double theta_p, int i_p) {
	// Return if particle index is out of range
	if (i_p >= n_body_sim->N)
		return;

	// Return if quadrupole moment is null
	if (J2_p == 0.0)
		return;

	// Get particle info
	reb_particle *particles = n_body_sim->particles;

	// Get constant values (not dependent on interacting particle)
	double Gmc = n_body_sim->G * particles[i_p].m;
	double constant = 1.5 * Gmc * J2_p * pow(R_p, 2.0);
	Vector xhat = { sin(theta_p) * cos(phi_p), sin(theta_p) * sin(phi_p), cos(
			theta_p) };
	Vector diff = { 1.0, 1.0, 3.0 };

	// Get interaction terms and accelerations for each particle
	for (int i = 0; i < n_body_sim->N; i++) {
		if (i != i_p) {
			// Get displacement
			Vector x_1 = { particles[i].x, particles[i].y, particles[i].z };
			Vector x_2 =
					{ particles[i_p].x, particles[i_p].y, particles[i_p].z };
			Vector r = x_1 - x_2;

			// Arbitrary softening of r^5
			double r5 = pow(r.len(), 5.0) + 0.5;

			// cos^2 -- r_vec dot x_hat / r
			double cos2 = pow(dot(r, xhat) / r.len(), 2.0);

			// Get acceleration of particle i
			Vector a = constant * r / r5 * (5.0 * cos2 - diff);

			// Note: Forces can also be written as
			// F = 3/2 G M J2_p R_p^2 r_vec / r^5 (5 z^2/r^2 - <1, 1, 3>)
			// Instead of z/r, use dot product of pole direction with displacement

			// Update acceleration of particle i
			particles[i].ax += a.getX();
			particles[i].ay += a.getY();
			particles[i].az += a.getZ();

			// Update acceleration of particle i_p
			double mfac = particles[i].m / particles[i_p].m;
			particles[i_p].ax -= a.getX() * mfac;
			particles[i_p].ay -= a.getY() * mfac;
			particles[i_p].az -= a.getZ() * mfac;
		}
	}
}

/***************************************/
/* Orbital properties of resolved body */
/***************************************/

// Compute orbital properties of resolved body (particles in range [i_low, i_high)) about all perturbing particles
// L is orbital angular momentum per unit mass
void compute_orb(reb_simulation *const n_body_sim, int i_low,
		int i_high, double *a, double *mean_motion, double *e, double *i,
		double *L) {
	// Get particle info
	reb_particle *particles = n_body_sim->particles;

	// Get total mass of resolved body
	double m_tot = 0.0;
	for (int i = i_low; i < i_high; i++) {
		m_tot += n_body_sim->particles[i].m; // mass of resolved body
	}

	// Get center of mass and velocity of resolved body
	Vector CoM = compute_com(n_body_sim, i_low, i_high);
	Vector CoV = compute_cov(n_body_sim, i_low, i_high);

	// Get low and high index for perturbing particles
	int i_pert_low = n_body_sim->N - num_perts;
	int i_pert_high = n_body_sim->N;

	// Get center of mass and velocity of perturbers
	Vector CoM0 = compute_com(n_body_sim, i_pert_low, i_pert_high);
	Vector CoV0 = compute_cov(n_body_sim, i_pert_low, i_pert_high);

	// Get total perturbing mass
	double m = 0.0;
	for (int i = i_pert_low; i < i_pert_high; i++)
		m += particles[i].m;

	// Displacement and relative velocity of resolved body
	Vector dx = CoM0 - CoM;
	Vector dv = CoV0 - CoV;

	// Total (resolved and perturbing) mass
	double M = m + m_tot;

	// Kinetic energy per unit mass
	double KE = 0.5 * pow(dv.len(), 2.0);

	// Gravitational potential energy per unit mass
	double GM = n_body_sim->G * M;
	double GPE = -GM / dx.len();

	// Semi-major axis
	double E = KE + GPE;
	*a = -0.5 * GM / E;

	// Mean motion
	*mean_motion = sqrt(GM / (*a * *a * *a));

	// Orbital angular momentum per unit mass
	Vector r_cross_v = cross(dx, dv);
	*L = r_cross_v.len();

	// Eccentricity
	double e2 = 1.0 - pow(*L, 2.0) / (*a * GM);
	if (e2 > 0.0) {
		*e = sqrt(e2);
	} else {
		*e = 0.0;
	}

	// Inclination
	*i = acos(r_cross_v.getZ() / *L);
}

/********/
/* Misc */
/********/

// Compute total mass of particles in particle range [i_low, i_high)
double sum_mass(reb_simulation *const n_body_sim, int i_low,
		int i_high) {
	// Get particle info
	reb_particle *particles = n_body_sim->particles;

	// Sum masses of each particle
	double m_tot = 0.0;
	for (int i = i_low; i < i_high; i++) {
		m_tot += particles[i].m;
	}
	return m_tot;
}
