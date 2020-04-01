#ifdef __cplusplus
# 	ifdef __GNUC__
#		define restrict __restrict__
#	else
#		define restrict
#	endif
#endif

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <cmath>
#include <time.h>
#include <sys/time.h>
#include <stdbool.h>
#include <iostream>
extern "C" {
#include "rebound.h"
}
#include "spring.h"
#include "tools.h"
#include "output.h"
#include "kepcart.h"
#include "matrix_math.h"

extern int NS; // Current number of springs
int NSmax = 0; // Max number of springs (size of array)

/********************/
/* Spring functions */
/********************/

// Delete spring at spring_index
void del_spring(struct reb_simulation *const n_body_sim, int spring_index) {
	// Overwrite spring at index with last spring, and reduce count (last spring is still present at original location)
	if (NS > 0) {
		springs[spring_index] = springs[NS - 1];
		NS--;
	}
}

// Add a spring connecting particle_1 and particle_2
// Check that a spring doesn't already exist between these two particles
// Set the natural distance of the spring to the current interparticle distance
// Spring constant is not scaled
// Return -1 if no spring because you're trying to connect a particle to itself
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
	for (int spring_index = 0; spring_index < NS; spring_index++) {
		if ((springs[spring_index].i == particle_low)
				&& (springs[spring_index].j == particle_high))
			return spring_index;
	}

	// No spring connects these two indices. Create one.
	spr.i = particle_low;
	spr.j = particle_high;
	spr.rs0 = spring_length(n_body_sim, spr); // rest spring length
	add_spring_helper(spr);
	return NS - 1; // index of new spring!
}

// Helper function that adds a spring to array
// Caution: no sanity checking
void add_spring_helper(spring spr) {
	// If max size is smaller than current size, increase
	while (NSmax <= NS) {
		NSmax += 128;
		springs = (spring*) realloc(springs, sizeof(spring) * NSmax);
	}

	// Add spring to end of array, update count
	springs[NS] = spr;
	NS++;
}

// Return spring length
double spring_length(struct reb_simulation *const n_body_sim, spring spr) {
	// Get simulation and spring information
	struct reb_particle *particles = n_body_sim->particles;
	int i = spr.i;
	int j = spr.j;

	// Get locations of endpoints of spring
	Vector x_i = { particles[i].x, particles[i].y, particles[i].z };
	Vector x_j = { particles[j].x, particles[j].y, particles[j].z };

	// Return length of spring
	return (x_i - x_j).len();
}

// Compute spring midpoint location from arbitrary center (Cartesian coordinates)
Vector spr_mid(struct reb_simulation *const n_body_sim, spring spr,
		Vector center) {
	// Get simulation and spring information
	struct reb_particle *particles = n_body_sim->particles;
	int i = spr.i;
	int j = spr.j;

	// Get locations of endpoints of spring
	Vector x_i = { particles[i].x, particles[i].y, particles[i].z };
	Vector x_j = { particles[j].x, particles[j].y, particles[j].z };

	// Calculate location of midpoint of spring from arbitrary center position
	return 0.5 * (x_i + x_j) - center;
}

// Return spring strain
double strain(struct reb_simulation *const n_body_sim, spring spr) {
	double len = spring_length(n_body_sim, spr); // Spring length
	return (len - spr.rs0) / spr.rs0; // Positive under extension
}

// Connect springs to all particles with interparticle distances less than max_dist apart for particle index range [i_low,i_high-1]
// Spring is added with rest length at current length
void connect_springs_dist(struct reb_simulation *const n_body_sim,
		double max_dist, int i_low, int i_high, spring spr) {
	// Get particle info from sim
	reb_particle* particles = n_body_sim->particles;

	// Return if low index is above high index
	if (i_high <= i_low)
		return;

	// Create all springs for near neighbors
	// Scan through particles from i_low to i_high-1 (i_high can't pair with itself)
	for (int i = i_low; i < i_high - 1; i++) {
		Vector x_i = {particles[i].x, particles[i].y, particles[i].z};

		// Scan through all possible pairs, without repetition
		for (int j = i + 1; j < i_high; j++) {
			Vector x_j = {particles[j].x, particles[j].y, particles[j].z};

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
	for (int i = 0; i < NS; i++) {
		springs[i].gamma = new_gamma;
	}
	std::cout.precision(2);
	std::cout << "\n gamma set to " << new_gamma << std::endl;
}

/*****************************************/
/* !!!!!!!!! Unsorted functions !!!!!!!! */
/*****************************************/

// compute the closest distance between any pair of particles
// with index in range [imin,imax-1]
double mindist(struct reb_simulation *const r, int imin, int imax) {
	double dist = 1e10;
	struct reb_particle *particles = r->particles;
	for (int i = imin; i < imax - 1; i++) {
		for (int j = i + 1; j < imax; j++) {
			double dx = particles[i].x - particles[j].x;
			double dy = particles[i].y - particles[j].y;
			double dz = particles[i].z - particles[j].z;
			double dr = sqrt(dx * dx + dy * dy + dz * dz);
			if (dr < dist)
				dist = dr;
		}
	}
	return dist;
}

// compute center of mass coordinates in particle range [il,ih)
Vector compute_com(struct reb_simulation *const n_body_sim, int il, int ih) {
	// Get simulation information and initialize sums
	struct reb_particle *particles = n_body_sim->particles;
	double m_tot = 0.0;
	Vector x_times_m(zero_vec);

	// Sum masses and location-weighted masses
	for (int i = il; i < ih; i++) {
		x_times_m += particles[i].m * Vector({ particles[i].x,
						particles[i].y, particles[i].z });
		m_tot += particles[i].m;
	}

	// Calculate center of mass
	return x_times_m / m_tot;
}

// compute total mass of particles in particle range [il,ih)
double sum_mass(struct reb_simulation *const r, int il, int ih) {
	struct reb_particle *particles = r->particles;
	double msum = 0.0;
	for (int i = il; i < ih; i++) {
		msum += particles[i].m;
	}
	return msum;
}

// compute center of velocity particles in particle range [il,ih)
// values returned in vxc, vyc, vzc
void compute_cov(struct reb_simulation *const r, int il, int ih, double *vxc, double *vyc, double *vzc) {
	double vxsum = 0.0;
	double vysum = 0.0;
	double vzsum = 0.0;
	double msum = 0.0;
	struct reb_particle *particles = r->particles;
	for (int i = il; i < ih; i++) {
		vxsum += particles[i].vx * particles[i].m;
		vysum += particles[i].vy * particles[i].m;
		vzsum += particles[i].vz * particles[i].m;
		msum += particles[i].m;
	}
	*vxc = vxsum / msum;
	*vyc = vysum / msum;
	*vzc = vzsum / msum;
}

// go to coordinate frame of body defined by vertices/particles [il,ih)
// mass weighted center of mass
// only coordinates changed,  particle velocities not changed
// all particles are shifted, not just the extended body
void centerbody(struct reb_simulation *const r, int il, int ih) {
	Vector CoM = compute_com(r, il, ih);
	for (int i = 0; i < (r->N); i++) { // all particles shifted
		r->particles[i].x -= CoM.getX();
		r->particles[i].y -= CoM.getY();
		r->particles[i].z -= CoM.getZ();
	}
}

// subtract center of velocity from the resolved body
// only changing particles in the resolved body [il,ih)
void subtractcov(struct reb_simulation *const r, int il, int ih) {
	double vxc = 0.0;
	double vyc = 0.0;
	double vzc = 0.0;
	compute_cov(r, il, ih, &vxc, &vyc, &vzc); // center of velocity of resolved body
	move_resolved(r, 0.0, 0.0, 0.0, -vxc, -vyc, -vzc, il, ih);
}

// subtract center of mass position from the resolved body
// only changing particles in the resolved body [il,ih)
void subtractcom(struct reb_simulation *const r, int il, int ih) {
	Vector CoM = compute_com(r, il, ih); // center of mass of resolved body
	move_resolved(r, -CoM.getX(), -CoM.getY(), -CoM.getZ(), 0.0, 0.0, 0.0, il, ih);
}

// calculate spring forces
// viscoelastic model, Clavet et al.  05
// "Particle-based Viscoelastic Fluid Simulation"
// Eurographics/ACM SIGGRAPH Symposium on Computer Animation (2005)
// K. Anjyo, P. Faloutsos (Editors)
// by Simon Clavet, Philippe Beaudoin, and Pierre Poulin 
#define L_EPS 0e-16; // softening for spring length
void spring_forces(struct reb_simulation *const r) {
	for (int i = 0; i < NS; i++) {  // spring forces
		double L = spring_length(r, springs[i]) + L_EPS
		; // spring length
		double rs0 = springs[i].rs0; // rest length
		int ii = springs[i].i;
		int jj = springs[i].j;
		double dx = r->particles[ii].x - r->particles[jj].x;
		double dy = r->particles[ii].y - r->particles[jj].y;
		double dz = r->particles[ii].z - r->particles[jj].z;
		double mii = r->particles[ii].m;
		double mjj = r->particles[jj].m;
		double ks = springs[i].ks;
		double fac = -ks * (L - rs0) / L; // L here to normalize direction
		// apply elastic forces
		// accelerations are force divided by mass
		r->particles[ii].ax += fac * dx / mii;
		r->particles[jj].ax -= fac * dx / mjj;
		r->particles[ii].ay += fac * dy / mii;
		r->particles[jj].ay -= fac * dy / mjj;
		r->particles[ii].az += fac * dz / mii;
		r->particles[jj].az -= fac * dz / mjj;

		// apply damping, depends on strain rate
		double ggamma = springs[i].gamma;
		if (ggamma > 0.0) {
			double dvx = r->particles[ii].vx - r->particles[jj].vx;
			double dvy = r->particles[ii].vy - r->particles[jj].vy;
			double dvz = r->particles[ii].vz - r->particles[jj].vz;
			double dLdt = (dx * dvx + dy * dvy + dz * dvz) / L;
			// divide dL/dt by L to get strain rate
			double mbar = mii * mjj / (mii + mjj); // reduced mass
			double dampfac = ggamma * mbar * dLdt / L;
			// factor L here to normalize dx,dy,dz
			r->particles[ii].ax -= dampfac * dx / mii;
			r->particles[jj].ax += dampfac * dx / mjj;
			r->particles[ii].ay -= dampfac * dy / mii;
			r->particles[jj].ay += dampfac * dy / mjj;
			r->particles[ii].az -= dampfac * dz / mii;
			r->particles[jj].az += dampfac * dz / mjj;
			// ggamma is in units of 1/time
			// force is gamma*dL/dt*mbar * dx/L = gamma*mbar*deps/dt*L*dx/L
		}

	}
}

// compute total potential energy in spring network
double spring_potential_energy(struct reb_simulation *const r) {
	double pe = 0.0;
	for (int i = 0; i < NS; i++) {
		double rs0 = springs[i].rs0; // rest length
		double L = spring_length(r, springs[i]) + L_EPS
		; // spring length
		double ks = springs[i].ks; // spring constant
		pe += ks * pow(L - rs0, 2.0) / 2.0;
	}
	return pe;
}

// compute total gravitational potential energy in resolved body 
double grav_potential_energy(struct reb_simulation *const r, int il, int ih) {
	double pesum = 0.0;
	for (int i = il; i < ih; i++) {
		for (int j = i + 1; j < ih; j++) {
			double dx = r->particles[i].x - r->particles[j].x;
			double dy = r->particles[i].y - r->particles[j].y;
			double dz = r->particles[i].z - r->particles[j].z;
			double invr = 1.0 / sqrt(dx * dx + dy * dy + dz * dz);
			double pe2 = r->particles[i].m * r->particles[j].m * invr * r->G;
			pesum += pe2;
		}
	}
	return -pesum; // note sign
}

//////////////////////////////////
// compute power lost in damping from a specific spring
// due to damping (viscoelasticity)
double dEdt(struct reb_simulation *const r, struct spring spr) {
	double ggamma = spr.gamma;
	if (ggamma == 0.0)
		return 0.0;

	double L = spring_length(r, spr) + L_EPS
	; // spring length
	int ii = spr.i;
	int jj = spr.j;
	double dx = r->particles[ii].x - r->particles[jj].x;
	double dy = r->particles[ii].y - r->particles[jj].y;
	double dz = r->particles[ii].z - r->particles[jj].z;
	// double dr = sqrt(dx*dx + dy*dy + dz*dz);
	double mii = r->particles[ii].m;
	double mjj = r->particles[jj].m;
	double mbar = mii * mjj / (mii + mjj); // reduced mass
	double dvx = r->particles[ii].vx - r->particles[jj].vx;
	double dvy = r->particles[ii].vy - r->particles[jj].vy;
	double dvz = r->particles[ii].vz - r->particles[jj].vz;
	double dLdt = (dx * dvx + dy * dvy + dz * dvz) / L;
	// divide dL/dt by L to get strain rate
	double de = ggamma * mbar * dLdt * dLdt;
	return de;  // units power de/dt as expected
	// we do not need to multiply by timestep to get de/dt
}

// compute power summing over all springs
double dEdt_total(struct reb_simulation *const r) {
	double sum = 0.0;
	for (int i = 0; i < NS; i++) {
		sum += dEdt(r, springs[i]);
	}
	return sum;
}

// zero all particle accelerations
void zero_accel(struct reb_simulation *const r) {
	for (int i = 0; i < (r->N); i++) {
		r->particles[i].ax = 0.0;
		r->particles[i].ay = 0.0;
		r->particles[i].az = 0.0;
	}
}

// create a ~uniform random particle distribution with total mass  = total_mass
// fill particles within a football shape given by semi- axes ax, by, cz
// spacing set by parameter dist: no closer particles allowed
// to make this uniform first create a sphere that is large
// enough to hold entire football
// this should reduce non-uniformity from surface effects
void rand_football_from_sphere(struct reb_simulation *r, double dist, double ax,
		double by, double cz, double total_mass) {
	double rhold = ax + 2 * dist; // size of sphere that we need to hold everything
	struct reb_particle pt;
	int npart = (int) (40.0 * pow(2.0 * rhold / dist, 3.0));
	// guess for number of random particles we need to generate
	pt.ax = 0.0;
	pt.ay = 0.0;
	pt.az = 0.0;
	pt.vx = 0.0;
	pt.vy = 0.0;
	pt.vz = 0.0;
	pt.m = 1.0;
	double particle_radius = dist / 2.0;
	pt.r = particle_radius;
	int N = r->N;
	int i0 = N;
	for (int i = 0; i < npart; i++) {
		pt.x = reb_random_uniform(-rhold, rhold);
		pt.y = reb_random_uniform(-rhold, rhold);
		pt.z = reb_random_uniform(-rhold, rhold);
		double x2 = pow(pt.x, 2.0);
		double y2 = pow(pt.y, 2.0);
		double z2 = pow(pt.z, 2.0);
		double rval = sqrt(x2 + y2 + z2);
		if (rval < rhold) { // within rhold
			int toonear = 0;  // is there a particle too nearby?
			int j = i0;
			N = r->N;
			while ((toonear == 0) && (j < N)) {
				double dx = pt.x - r->particles[j].x;
				double dy = pt.y - r->particles[j].y;
				double dz = pt.z - r->particles[j].z;
				double dr = sqrt(dx * dx + dy * dy + dz * dz);
				if (dr < dist)
					toonear = 1;
				j++;
			}
			if (toonear == 0)
				reb_add(r, pt);
			// only add particle if not near any other
		}
	}
	N = r->N;
	// now remove all particles outside our football
	int imax = N;
	for (int i = i0; i < imax; i++) {
		double x = r->particles[i].x;
		double y = r->particles[i].y;
		double z = r->particles[i].z;
		double xa2 = pow(x / ax, 2.0);
		double ya2 = pow(y / by, 2.0);
		double za2 = pow(z / cz, 2.0);
		double rval = sqrt(xa2 + ya2 + za2);
		if (rval > 1.0) { // outside ellipsoid
			reb_remove(r, i, 0);  // remove particle
			i--; // we copy in a particle and we need to look at it
			imax = r->N;
		}
	}
	N = r->N;

	// adjust mass of each particle so that sums to desired total mass
	double particle_mass = total_mass / (N - i0);
	// fix masses
	for (int ii = i0; ii < N; ii++) { // all particles!
		r->particles[ii].m = particle_mass;
	}
	double md = mindist(r, i0, N);
	printf("rand_football_from_sphere: Nparticles=%d min_d=%.2f\n", N - i0, md);
}

// create a ~uniform random particle distribution with total mass  = total_mass
// fill particles within a football shape given by semi- axes ax, by, cz
// spacing set by parameter dist: no closer particles allowed
void rand_football(struct reb_simulation *const r, double dist, double ax,
		double by, double cz, double total_mass) {
	struct reb_particle pt;
	int npart = 40 * pow(2.0 * ax / dist, 3.0);
	// guess for number of random particles we need to generate
	printf("rand_football: npart=%d\n", npart);
	pt.ax = 0.0;
	pt.ay = 0.0;
	pt.az = 0.0;
	pt.vx = 0.0;
	pt.vy = 0.0;
	pt.vz = 0.0;
	pt.m = 1.0;
	double particle_radius = dist / 2.0;
	pt.r = particle_radius / 2.0;  // XXXXxxxxx temp
	int N = r->N;
	int i0 = N;
	for (int i = 0; i < npart; i++) {
		pt.x = reb_random_uniform(-ax, ax);
		pt.y = reb_random_uniform(-by, by);
		pt.z = reb_random_uniform(-cz, cz);
		double xa2 = pow(pt.x / ax, 2.0);
		double yb2 = pow(pt.y / by, 2.0);
		double zc2 = pow(pt.z / cz, 2.0);
		double rval = sqrt(xa2 + yb2 + zc2);
		if (rval < 1.0) { // within football
			int toonear = 0;  // is there a particle too nearby?
			int j = i0;
			N = r->N;
			while ((toonear == 0) && (j < N)) {
				double dx = pt.x - r->particles[j].x;
				double dy = pt.y - r->particles[j].y;
				double dz = pt.z - r->particles[j].z;
				double dr = sqrt(dx * dx + dy * dy + dz * dz);
				if (dr < dist)
					toonear = 1;
				j++;
			}
			if (toonear == 0)
				reb_add(r, pt);
			// only add particle if not near any other
		}
	}
	N = r->N;

	// adjust mass of each particle so that sums to desired total mass
	double particle_mass = total_mass / (N - i0);
	// fix masses
	for (int ii = i0; ii < N; ii++) { // all particles!
		r->particles[ii].m = particle_mass;
	}
	double md = mindist(r, i0, N);
	printf("rand_football: Nparticles=%d min_d=%.2f\n", N - i0, md);
}

// create a ~uniform random particle distribution with total mass  = total_mass
// fill particles within a rectangle shape given by total lengths ax, by, cz
// spacing set by parameter dist: no closer particles allowed
void rand_rectangle(struct reb_simulation *const r, double dist, double ax,
		double by, double cz, double total_mass) {
	struct reb_particle pt;
	int npart = 40.0 * ax * by * cz / pow(dist, 3.0);
	// guess for number of random particles we need to generate
	printf("rand_rectangle: npart=%d\n", npart);
	pt.ax = 0.0;
	pt.ay = 0.0;
	pt.az = 0.0;
	pt.vx = 0.0;
	pt.vy = 0.0;
	pt.vz = 0.0;
	pt.m = 1.0;
	double particle_radius = dist / 2.0;
	pt.r = particle_radius / 3.0;  // XXXXxxxxx temp
	int N = r->N;
	int i0 = N;
	for (int i = 0; i < npart; i++) {
		pt.x = reb_random_uniform(-ax / 2, ax / 2);
		pt.y = reb_random_uniform(-by / 2, by / 2);
		pt.z = reb_random_uniform(-cz / 2, cz / 2);
		{    // within rectangle
			int toonear = 0;  // is there a particle too nearby?
			int j = i0;
			N = r->N;
			while ((toonear == 0) && (j < N)) {
				double dx = pt.x - r->particles[j].x;
				double dy = pt.y - r->particles[j].y;
				double dz = pt.z - r->particles[j].z;
				double dr = sqrt(dx * dx + dy * dy + dz * dz);
				if (dr < dist)
					toonear = 1;
				j++;
			}
			if (toonear == 0)
				reb_add(r, pt);
			// only add particle if not near any other
		}
	}
	N = r->N;

	// adjust mass of each particle so that sums to desired total mass
	double particle_mass = total_mass / (N - i0);
	// fix masses
	for (int ii = i0; ii < N; ii++) { // all particles!
		r->particles[ii].m = particle_mass;
	}
	double md = mindist(r, i0, N);
	printf("rand_rectangle: Nparticles=%d min_d=%.2f\n", N - i0, md);
}

// create a ~uniform random particle distribution with total mass  = total_mass
// fill particles within a rectangle shape given by total lengths ax, by in 2d
// spacing set by parameter dist: no closer particles allowed
// set z=0 everywhere
void rand_rectangle_2d(struct reb_simulation *const r, double dist, double ax,
		double by, double total_mass) {
	struct reb_particle pt;
	int npart = 40.0 * ax * by / pow(dist, 2.0);
	printf("rand_rectangle_2d: npart=%d\n", npart);
	pt.ax = 0.0;
	pt.ay = 0.0;
	pt.az = 0.0;
	pt.vx = 0.0;
	pt.vy = 0.0;
	pt.vz = 0.0;
	pt.m = 1.0;
	pt.x = 0.0;
	pt.y = 0.0;
	pt.z = 0.0;
	double particle_radius = dist / 2.0;
	pt.r = particle_radius / 3.0;  // XXXXxxxxx temp
	int N = r->N;
	int i0 = N;
	for (int i = 0; i < npart; i++) {
		pt.x = reb_random_uniform(-ax / 2, ax / 2);
		pt.y = reb_random_uniform(-by / 2, by / 2);
		{    // within rectangle
			int toonear = 0;  // is there a particle too nearby?
			int j = i0;
			N = r->N;
			while ((toonear == 0) && (j < N)) {
				double dx = pt.x - r->particles[j].x;
				double dy = pt.y - r->particles[j].y;
				double dz = pt.z - r->particles[j].z;
				double dr = sqrt(dx * dx + dy * dy + dz * dz);
				if (dr < dist)
					toonear = 1;
				j++;
			}
			if (toonear == 0)
				reb_add(r, pt);
			// only add particle if not near any other
		}
	}
	N = r->N;

	// adjust mass of each particle so that sums to desired total mass
	double particle_mass = total_mass / (N - i0);
	// fix masses
	for (int ii = i0; ii < N; ii++) { // all particles!
		r->particles[ii].m = particle_mass;
	}
	double md = mindist(r, i0, N);
	printf("rand_rectangle: Nparticles=%d min_d=%.2f\n", N - i0, md);
}

// create a ~uniform random particle distribution with total mass  = total_mass
// fill particles within a cone shape given by base radius rb, and slope h/rb
// spacing set by parameter dist: no closer particles allowed
void rand_cone(struct reb_simulation *const r, double dist, double rb,
		double slope, double total_mass) {
	struct reb_particle pt;
	int npart = 40 * pow(2.0 * rb / dist, 3.0);
	// guess for number of random particles we need to generate
	pt.ax = 0.0;
	pt.ay = 0.0;
	pt.az = 0.0;
	pt.vx = 0.0;
	pt.vy = 0.0;
	pt.vz = 0.0;
	pt.m = 1.0;
	double particle_radius = dist / 2.0;
	pt.r = particle_radius / 2.0;  // XXXXxxxxx temp
	int N = r->N;
	int i0 = N;
	for (int i = 0; i < npart; i++) {
		pt.x = reb_random_uniform(-rb, rb);
		pt.y = reb_random_uniform(-rb, rb);
		pt.z = reb_random_uniform(-rb, rb);
		double radius = sqrt(pt.x * pt.x + pt.y * pt.y);
		double zval = rb * slope - slope * radius; // h - slope*r
		double zratio = fabs(pt.z / zval);
		if ((radius < rb) && (zratio < 1.0)) { // within cone
			int toonear = 0;  // is there a particle too nearby?
			int j = i0;
			N = r->N;
			while ((toonear == 0) && (j < N)) {
				double dx = pt.x - r->particles[j].x;
				double dy = pt.y - r->particles[j].y;
				double dz = pt.z - r->particles[j].z;
				double dr = sqrt(dx * dx + dy * dy + dz * dz);
				if (dr < dist)
					toonear = 1;
				j++;
			}
			if (toonear == 0)
				reb_add(r, pt);
			// only add particle if not near any other
		}
	}
	N = r->N;

	// adjust mass of each particle so that sums to desired total mass
	double particle_mass = total_mass / (N - i0);
	// fix masses
	for (int ii = i0; ii < N; ii++) { // all particles!
		r->particles[ii].m = particle_mass;
	}
	double md = mindist(r, i0, N);
	printf("rand_cone: Nparticles=%d min_d=%.2f\n", N - i0, md);
}

// return mean spring constant of all springs
double mean_ks(struct reb_simulation *const r, int type) {
	double sum = 0.0;
	int nn = 0;
	for (int i = 0; i < NS; i++) {
		sum += springs[i].ks;
		nn++;
	}
	return sum / nn;
}

// compute distance of particle with index ii from coordinate (xc,yc,zc)
double rad_com(struct reb_simulation *const r, int ii, double xc, double yc,
		double zc) {
	double dx = r->particles[ii].x - xc;
	double dy = r->particles[ii].y - yc;
	double dz = r->particles[ii].z - zc;
	double rad = sqrt(dx * dx + dy * dy + dz * dz);
	return rad;
}

// compute Young's modulus of springs 
// using midpoints in radial range 
//    from center of mass [rmin,rmax]
// using equation 20 by Kot et al. 2014  sum_i k_iL_i^2/(6V)
// uses rest lengths
// only computes center of mass using particles index range [il,ih)
double Young_mesh(struct reb_simulation *const n_body_sim, int il, int ih, double rmin,
		double rmax) {
	double sum = 0.0;
	Vector CoM = compute_com(n_body_sim, il, ih); // center of mass coords for particles in range
	for (int i = 0; i < NS; i++) {
		Vector x_mid = spr_mid(n_body_sim, springs[i], CoM);

		double r_c = x_mid.len(); // center of spring
		if ((r_c < rmax) && (r_c > rmin)) {
			double ks = springs[i].ks;
			double Li = springs[i].rs0;
			sum += ks * Li * Li;
		}
	}
	double volume = (4.0 * M_PI / 3.0) * (pow(rmax, 3.0) - pow(rmin, 3.0)); // in shell
	double E = sum / (6.0 * volume); // equation 20 by Kot et al. 2014
	return E; // return Young's modulus
}

// alternate routine using every spring
double Young_mesh_big(struct reb_simulation *const r, int il, int ih) {
	double sum = 0.0;
	for (int i = 0; i < NS; i++) {
		double ks = springs[i].ks;
		double Li = springs[i].rs0;
		sum += ks * Li * Li;
	}
	double volume = 4.0 * M_PI / 3.0; // sphere of radius 1 assumed?
	double E = sum / (6.0 * volume); // equation 20 by Kot et al. 2014
	return E; // return Young's modulus
}

// compute mean rest length of springs 
double mean_L(struct reb_simulation *const r) {
	double sum = 0.0;
	for (int i = 0; i < NS; i++) {
		sum += springs[i].rs0;
	}
	return sum / NS;
}

// start the body spinning particles indexes in [il,ih)
// with spin value omegax,omegay,omegaz,  about center of mass
// spin vector is omegax,omegay,omegaz
// does not change center of mass coordinates or velocity
void spin(struct reb_simulation *const r, int il, int ih, Vector omega) {
	Vector CoM = compute_com(r, il, ih); // compute center of mass
	double vxc, vyc, vzc;
	compute_cov(r, il, ih, &vxc, &vyc, &vzc);

	// Operate only if spin is large enough
	if (omega.len() > 1e-5) {
		for (int i = il; i < ih; i++) {

			Vector dx = { r->particles[i].x - CoM.getX(), r->particles[i].y - CoM.getY(), r->particles[i].z - CoM.getZ() };
			Vector r_cross_omega = cross(dx, omega);

			// set it spinning with respect to center of mass
			r->particles[i].vx = vxc + r_cross_omega.getX();
			r->particles[i].vy = vyc + r_cross_omega.getY();
			r->particles[i].vz = vzc + r_cross_omega.getZ();
		}
	}

}

// make a binary with two masses m1,m2 spinning with vector omega
// masses are separated by distance sep
// connect two masses with spring with values given by spring_vals
// center of mass is set to origin
void make_binary_spring(struct reb_simulation *const r, double m1, double m2,
		double sep, Vector omega, struct spring spring_vals) {
	const int i_low = r->N;
	struct reb_particle pt;
	pt.ax = 0.0;
	pt.ay = 0.0;
	pt.az = 0.0;
	pt.vx = 0.0;
	pt.vy = 0.0;
	pt.vz = 0.0;
	pt.y = 0.0;
	pt.z = 0.0;
	pt.m = m1;
	pt.x = sep * m2 / (m1 + m2);
	pt.r = sep * 0.3;
	reb_add(r, pt);
	pt.m = m2;
	pt.x = -sep * m1 / (m1 + m2);
	pt.r *= pow(m1 / m2, 0.33333);
	reb_add(r, pt);
	int i_high = i_low + 2;
	spin(r, i_low, i_high, omega); // spin it
	connect_springs_dist(r, sep * 1.1, i_low, i_high, spring_vals); // add spring

}

// return angle between body i and body j in xy plane
double get_angle(struct reb_simulation *const r, int i, int j) {
	double dx = r->particles[i].x - r->particles[j].x;
	double dy = r->particles[i].y - r->particles[j].y;
	double theta = atan2(dy, dx);
	return theta;
}

// compute angular momentum vector of a body with particles in range [il,ih)
// with respect to its center of mass position and velocity 
// can measure angular momentum of the entire system if il=0 and ih=N
// but with respect to center of mass of entire system
Vector measure_L(struct reb_simulation *const r, int il, int ih) {
	struct reb_particle *particles = r->particles;
	Vector CoM = compute_com(r, il, ih);
	double vxc = 0.0;
	double vyc = 0.0;
	double vzc = 0.0;
	compute_cov(r, il, ih, &vxc, &vyc, &vzc);

	Vector L(zero_vec);
	for (int i = il; i < ih; i++) {
		Vector dx = { particles[i].x - CoM.getX(), particles[i].y
								- CoM.getY(), particles[i].z - CoM.getZ() };
		Vector dv = { particles[i].vx - vxc, particles[i].vy
								- vyc, particles[i].vz - vzc };
		Vector dL = cross(dx, dv);
		L += particles[i].m * dL; // angular momentum vector
	}

	return Vector(L);

}

// compute angular momentum vector of all particles in range [il,ih)
// with respect to origin
void measure_L_origin(struct reb_simulation *const r, int il, int ih,
		double *llx, double *lly, double *llz) {
	struct reb_particle *particles = r->particles;

	double lx = 0.0;
	double ly = 0.0;
	double lz = 0.0;
	for (int i = il; i < ih; i++) {
		double dx = (particles[i].x);
		double dy = (particles[i].y);
		double dz = (particles[i].z);
		double dvx = (particles[i].vx);
		double dvy = (particles[i].vy);
		double dvz = (particles[i].vz);
		lx += particles[i].m * (dy * dvz - dz * dvy); // angular momentum vector
		ly += particles[i].m * (dz * dvx - dx * dvz);
		lz += particles[i].m * (dx * dvy - dy * dvx);
	}
	*llx = lx;
	*lly = ly;
	*llz = lz;
}

// compute the moment of inertia tensor of a body with particle indices [il,ih)
// with respect to center of mass
Matrix mom_inertia(struct reb_simulation *const n_body_sim, int i_low,
		int i_high) {

	reb_particle *particles = n_body_sim->particles;
	Matrix inertia(zero_mat);

	// Calculate center of mass of specified body
	Vector CoM = compute_com(n_body_sim, i_low, i_high);

	// Calculate moment of inertia of particles
	for (int i = i_low; i < i_high; i++) {
		Vector dx ={ particles[i].x - CoM.getX(), particles[i].y
								- CoM.getY(), particles[i].z - CoM.getZ() };
		inertia += particles[i].m * pow(dx.len(), 2.0) * I;
		inertia -= particles[i].m * outer(dx, dx);
	}
	return inertia;
}

// compute orbital properties of body
// resolved body indexes [il,ih)
// primary mass is at index im1  
// and compute: 
//   mean-motion:nn, semi-major axis:aa 
//   eccentricity:ee and inclination:ii
//   LL orbital angular momentum per unit mass
void compute_semi(struct reb_simulation *const r, int il, int ih, int im1,
		double *aa, double *meanmo, double *ee, double *ii, double *LL) {
	struct reb_particle *particles = r->particles;
	static int first = 0;
	static double tm = 0.0;  // its mass
	if (first == 0) { // only calculate once
		for (int i = il; i < ih; i++)
			tm += r->particles[i].m; // mass of resolved body
		first = 1;
	}
	Vector CoM = compute_com(r, il, ih); // center of mass of resolved body
	double vxc = 0.0;
	double vyc = 0.0;
	double vzc = 0.0;
	compute_cov(r, il, ih, &vxc, &vyc, &vzc); // center of velocity of resolved body

	// int im1 = r->N -1; // index for primary perturber
	double x0 = particles[im1].x;
	double vx0 = particles[im1].vx;
	double y0 = particles[im1].y;
	double vy0 = particles[im1].vy;
	double z0 = particles[im1].z;
	double vz0 = particles[im1].vz;
	double dv2 = pow(vx0 - vxc, 2.0) + // square of velocity difference
			pow(vy0 - vyc, 2.0) + pow(vz0 - vzc, 2.0);
	double m1 = particles[im1].m;
	double MM = m1 + tm; // total mass
	// double mu = tm*m1/MM;  // reduced mass
	double GMM = r->G * MM;
	double ke = 0.5 * dv2; // kinetic energy /mu  (per unit mass)
	double dr = sqrt(pow(x0 - CoM.getX(), 2.0) +  // distance between
			pow(y0 - CoM.getY(), 2.0) + pow(z0 - CoM.getZ(), 2.0));
	double pe = -GMM / dr; // potential energy/mu, interaction term  tm*m1 = GM*mu
	double E = ke + pe;  // total energy per unit mass
	double a = -0.5 * GMM / E; // semi-major axis
	*aa = a;
	*meanmo = sqrt(GMM / (a * a * a)); // mean motion
	// printf("dr=%.2f dv2=%.2f\n",dr,dv2);
	// compute orbital angular momentum
	double dx = x0 - CoM.getX();
	double dy = y0 - CoM.getY();
	double dz = z0 - CoM.getZ();
	double dvx = vx0 - vxc;
	double dvy = vy0 - vyc;
	double dvz = vz0 - vzc;
	double lx = dy * dvz - dz * dvy;
	double ly = dz * dvx - dx * dvz;
	double lz = dx * dvy - dy * dvx;
	double ltot = sqrt(lx * lx + ly * ly + lz * lz);
	*LL = ltot; // orbital angular momentum per unit mass
	double e2 = 1.0 - ltot * ltot / (a * GMM);
	*ee = 0.0;
	if (e2 > 0.0)
		*ee = sqrt(e2); // eccentricity
	*ii = acos(lz / ltot); // inclination ? if lz==ltot then is zero

}

// compute the orbital angular momentum vector of a resolved body
// resolved body at indexes [il,ih)
// returns orbital angular momentum vector
// computes orbital angular momentum about central mass if there is one mass
// otherwise computes it about the center of mass of all the perturbers
void compute_Lorb(struct reb_simulation *const r, int il, int ih, int npert,
		double *llx, double *lly, double *llz) {
	struct reb_particle *particles = r->particles;
	static int first = 0;
	static double tm = 0.0;  // its mass
	if (first == 0) { // only calculate once
		for (int i = il; i < ih; i++)
			tm += r->particles[i].m; // mass of resolved body
		first = 1;
	}
	Vector CoM = compute_com(r, il, ih); // center of mass of resolved body
	double vxc = 0.0;
	double vyc = 0.0;
	double vzc = 0.0;
	compute_cov(r, il, ih, &vxc, &vyc, &vzc); // center of velocity of resolved body

	Vector CoM0;
	double vx0 = 0.0;
	double vy0 = 0.0;
	double vz0 = 0.0;
	if (npert == 1) {
		int im1 = r->N - 1; // index for primary perturber
		CoM0 = { particles[im1].x, particles[im1].y, particles[im1].z };
		vx0 = particles[im1].vx;
		vy0 = particles[im1].vy;
		vz0 = particles[im1].vz;
	} else {
		int iml = r->N - npert; // index range for perturbing masses
		int imh = r->N;
		CoM0 = compute_com(r, iml, imh); // center of mass of perturbing bodies
		compute_cov(r, iml, imh, &vx0, &vy0, &vz0); // center of velocity of perturbing bodies
	}

	Vector dx = CoM0 - CoM;
	Vector dv = {vx0 - vxc, vy0 - vyc, vz0 - vzc};
	Vector L = cross(dx, dv);
	*llx = L.getX();
	*lly = L.getY();
	*llz = L.getZ();
}

// compute rotational kinetic energy
// resolved body indexes [il,ih)
// this includes kinetic energy in vibrations (so not just rotational)
double compute_rot_kin(struct reb_simulation *const r, int il, int ih) {
	struct reb_particle *particles = r->particles;
	double vxc = 0.0;
	double vyc = 0.0;
	double vzc = 0.0;
	compute_cov(r, il, ih, &vxc, &vyc, &vzc); // center of velocity of resolved body
	double kesum = 0.0;
	for (int i = il; i < ih; i++) {
		double vx = particles[i].vx - vxc;
		double vy = particles[i].vy - vyc;
		double vz = particles[i].vz - vzc;
		double v2 = vx * vx + vy * vy + vz * vz;
		double ke = 0.5 * v2 * particles[i].m; // kinetic energy of particle
		kesum += ke;
	}
	return kesum;
}

// compute orbital properties of body, with respect
// to center of mass of a binary with two masses at index [N-1] and N-2 (or up to npert)
// resolved body indexes [il,ih)
// and compute: 
//   mean-motion:nn, semi-major axis:aa 
//   eccentricity:ee and inclination:ii
//   LL orbital angular momentum per unit mass
void compute_semi_bin(struct reb_simulation *const r, int il, int ih, int npert,
		double *aa, double *meanmo, double *ee, double *ii, double *LL) {
	struct reb_particle *particles = r->particles;
	static int first = 0;
	static double tm = 0.0;  // its mass
	if (first == 0) { // only calculate once
		for (int i = il; i < ih; i++)
			tm += r->particles[i].m; // mass of resolved body
		first = 1;
	}
	Vector CoM = compute_com(r, il, ih); // center of mass of resolved body
	double vxc = 0.0;
	double vyc = 0.0;
	double vzc = 0.0;
	compute_cov(r, il, ih, &vxc, &vyc, &vzc); // center of velocity of resolved body

	Vector CoM0;
	double vxc0 = 0.0;
	double vyc0 = 0.0;
	double vzc0 = 0.0;
	double m1 = 0.0;
	if (npert == 1) {
		int im1 = r->N - 1; // index for primary perturber
		CoM0 = {particles[im1].x, particles[im1].y, particles[im1].z};
		vxc0 = particles[im1].vx;
		vyc0 = particles[im1].vy;
		vzc0 = particles[im1].vz;
		m1 = particles[im1].m;
	} else {
		int iml = r->N - npert; // index range for perturbing masses
		int imh = r->N;
		CoM0 = compute_com(r, iml, imh); // center of mass of perturbing bodies
		compute_cov(r, iml, imh, &vxc0, &vyc0, &vzc0); // center of velocity of perturbing body
		for (int i = iml; i < imh; i++)
			m1 += particles[i].m; // total perturbing mass
	}

	double dv2 = pow(vxc0 - vxc, 2.0) + // square of velocity difference
			pow(vyc0 - vyc, 2.0) + pow(vzc0 - vzc, 2.0);
	double MM = m1 + tm; // total mass
	// double mu = tm*m1/MM;  // reduced mass
	double GMM = r->G * MM;
	double ke = 0.5 * dv2; // kinetic energy /mu  (per unit mass)
	double dr = (CoM0 - CoM).len();  // distance between
	double pe = -GMM / dr; // potential energy/mu, interaction term  tm*m1 = GM*mu
	double E = ke + pe;  // total energy
	double a = -0.5 * GMM / E; // semi-major axis
	*aa = a;
	*meanmo = sqrt(GMM / (a * a * a)); // mean motion

	// compute orbital angular momentum
	Vector dx = CoM0 - CoM;
	Vector dv = {vxc0 - vxc, vyc0 - vyc, vzc0 - vzc};
	Vector L = cross(dx, dv);
	double ltot = L.len();
	*LL = ltot; // orbital angular momentum per unit mass
	double e2 = fabs(1.0 - ltot * ltot / (a * GMM));
	*ee = sqrt(e2 + 1e-16); // eccentricity
	*ii = acos(L.getZ() / ltot); // inclination ? if lz==ltot then is zero

}

// sum total momentum of particles
void total_mom(struct reb_simulation *const r, int il, int ih, double *ppx,
		double *ppy, double *ppz) {
	double px = 0.0;
	double py = 0.0;
	double pz = 0.0;
	for (int i = il; i < ih; i++) {
		px += r->particles[i].m * r->particles[i].vx;
		py += r->particles[i].m * r->particles[i].vy;
		pz += r->particles[i].m * r->particles[i].vz;
	}
	*ppx = px;
	*ppy = py;
	*ppz = pz;

}

// using Euler angles rotate all particles [il,ih) about origin
// rotate both position and velocities
void rotate_origin(struct reb_simulation *const r, int il, int ih, double alpha,
		double beta, double ggamma) {
	struct reb_particle *particles = r->particles;
	for (int i = il; i < ih; i++) {
		double x0 = particles[i].x;
		double y0 = particles[i].y;
		double z0 = particles[i].z;
		double x1 = x0 * cos(alpha) - y0 * sin(alpha); // rotate about z axis in xy plane
		double y1 = x0 * sin(alpha) + y0 * cos(alpha);
		double z1 = z0;
		double x2 = x1;                      // rotate about x' axis in yz plane
		double y2 = y1 * cos(beta) - z1 * sin(beta);
		double z2 = y1 * sin(beta) + z1 * cos(beta);
		double x3 = x2 * cos(ggamma) - y2 * sin(ggamma); // rotate about z'' axis in xy plane
		double y3 = x2 * sin(ggamma) + y2 * cos(ggamma);
		double z3 = z2;
		particles[i].x = x3;
		particles[i].y = y3;
		particles[i].z = z3;

		double vx0 = particles[i].vx;
		double vy0 = particles[i].vy;
		double vz0 = particles[i].vz;
		double vx1 = vx0 * cos(alpha) - vy0 * sin(alpha); // rotate about z axis in xy plane
		double vy1 = vx0 * sin(alpha) + vy0 * cos(alpha);
		double vz1 = vz0;
		double vx2 = vx1;                    // rotate about x' axis in yz plane
		double vy2 = vy1 * cos(beta) - vz1 * sin(beta);
		double vz2 = vy1 * sin(beta) + vz1 * cos(beta);
		double vx3 = vx2 * cos(ggamma) - vy2 * sin(ggamma); // rotate about z'' axis in xy plane
		double vy3 = vx2 * sin(ggamma) + vy2 * cos(ggamma);
		double vz3 = vz2;
		particles[i].vx = vx3;
		particles[i].vy = vy3;
		particles[i].vz = vz3;
	}
}

// using Euler angles rotate a body with particle indices [il,ih)
// about center of mass
// rotate both position and velocities
// center of mass position and velocity is not changed
void rotate_body(struct reb_simulation *const r, int il, int ih, double alpha,
		double beta, double gamma) {
	struct reb_particle *particles = r->particles;

	// Set rotation matrices
	Matrix Rz1 = {{cos(alpha), -sin(alpha), 0}, {sin(alpha), cos(alpha), 0}, {0, 0, 1}}; // Rotate about z
	Matrix Rx = {{1, 0, 0}, {0, cos(beta), -sin(beta)}, {0, sin(beta), cos(beta)}}; // Rotate about x'
	Matrix Rz2 = {{cos(gamma), -sin(gamma), 0}, {sin(gamma), cos(gamma), 0}, {0, 0, 1}}; // Rotate about z''

	Vector CoM = compute_com(r, il, ih);
	double vxc, vyc, vzc;
	compute_cov(r, il, ih, &vxc, &vyc, &vzc);

	for (int i = il; i < ih; i++) {

		Vector x0 = {particles[i].x, particles[i].y, particles[i].z};
		x0 -= CoM;

		Vector x_rot = Rz1*Rx*Rz2*x0;

		particles[i].x = x_rot.getX() + CoM.getX();
		particles[i].y = x_rot.getY() + CoM.getY();
		particles[i].z = x_rot.getZ() + CoM.getZ();

		Vector v0 = {particles[i].vx - vxc, particles[i].vy - vyc, particles[i].vz - vzc};
		Vector v_rot = Rz1*Rx*Rz2*v0;

		particles[i].vx = v_rot.getX() + vxc;
		particles[i].vy = v_rot.getY() + vyc;
		particles[i].vz = v_rot.getZ() + vzc;
	}
}

// for testing above rotations
void rotate_vector(double x, double y, double z, double *xr, double *yr,
		double *zr, double alpha, double beta, double ggamma) {
	double x0 = x;
	double y0 = y;
	double z0 = z;
	double x1 = x0 * cos(alpha) - y0 * sin(alpha); // rotate about z axis in xy plane
	double y1 = x0 * sin(alpha) + y0 * cos(alpha);
	double z1 = z0;
	double x2 = x1;                          // rotate about x' axis in yz plane
	double y2 = y1 * cos(beta) - z1 * sin(beta);
	double z2 = y1 * sin(beta) + z1 * cos(beta);
	double x3 = x2 * cos(ggamma) - y2 * sin(ggamma); // rotate about z'' axis in xy plane
	double y3 = x2 * sin(ggamma) + y2 * cos(ggamma);
	double z3 = z2;
	*xr = x3;
	*yr = y3;
	*zr = z3;
}

// compute spin vector of a body with indices [il,ih)
//  spin vector  is omx, omy, omz 
//    computed using inverse of moment of inertia matrix
// also return eigenvalues of moment of inertia  matrix
// order big>=middle>=small (
Vector body_spin(struct reb_simulation *const r, int il, int ih, double eigs[3]) {

	//compute moment of inertia matrix
	Matrix inertia_mat = mom_inertia(r, il, ih);

	// compute inverse of moment of inertia matrix
	Matrix inv_mom_inert = inverse(inertia_mat);

	// compute angular momentum vector body with respect to its center
	// of mass position and velocity
	Vector L = measure_L(r, il, ih);

	// compute spin vector
	eigenvalues(inertia_mat, eigs);
	return inv_mom_inert*L;
}

// adjust ks and gamma  and k_heat
// for springs with midpoints between rmin and rmax of center of mass
void adjust_ks(struct reb_simulation *const r, int npert, double ksnew,
		double gammanew, double kheatnew, double rmin, double rmax) {
	// struct reb_particle* particles = r->particles;
	int il = 0;
	int ih = r->N - npert;
	Vector CoM = compute_com(r, il, ih);
	int NC = 0;
	for (int i = 0; i < NS; i++) {
		// compute spring mid point from central position
		Vector x_mid = spr_mid(r, springs[i], CoM);
		double rmid = x_mid.len();
		if ((rmid >= rmin) && (rmid <= rmax)) {
			springs[i].ks = ksnew;
			springs[i].gamma = gammanew;
			springs[i].k_heat = kheatnew;
			NC++;
		}
	}
	printf("adjust_ks: number of springs changed %d\n", NC);
	printf("adjust_ks: faction of springs changed %.3f\n",
			(double) NC / (double) NS);

}

// adjust ks and gamma and kheat
// for springs with midpoints within or outside of ellipsoid set 
// by (x-x0)^2/a^2 + (y-y0)^2/b^2 + (z-z0)^2/c^2 = 1
// inside==1 then inside ellipsoid
// inside==0 then outside
void adjust_ks_abc(struct reb_simulation *const r, int npert, double ksnew,
		double gammanew, double kheatnew, double a, double b, double c,
		Vector x0, int inside) {

	int il = 0;
	int ih = r->N - npert;
	Vector CoM = compute_com(r, il, ih);
	int NC = 0;
	for (int i = 0; i < NS; i++) {
		// compute spring mid point from central position
		Vector x_mid = spr_mid(r, springs[i], CoM);
		Vector x = x_mid - x0;
		double rmid2 = pow(x.getX() / a, 2.0)
					 + pow(x.getY() / b, 2.0)
					 + pow(x.getZ() / c, 2.0);
		if ((rmid2 <= 1.0) && (inside == 1)) {
			springs[i].ks = ksnew;
			springs[i].gamma = gammanew;
			springs[i].k_heat = kheatnew;
			NC++;
		}
		if ((rmid2 >= 1.0) && (inside == 0)) {
			springs[i].ks = ksnew;
			springs[i].gamma = gammanew;
			springs[i].k_heat = kheatnew;
			NC++;
		}
	}
	printf("adjust_ks_abc: number of springs changed %d\n", NC);
}

// adjust ks and gamma and kheat by a factor *= 1.0+fac*cos(mm*(phi-phi0))
// here phi is the angle in the xy plane
// for springs with midpoints within or outside of ellipsoid set 
// by (x-x0)^2/a^2 + (y-y0)^2/b^2 + (z-z0)^2/c^2 = 1
// inside==1 then inside ellipsoid
// inside==0 then outside
void adjust_ks_abc_fac(struct reb_simulation *const r, int npert, double ks_fac,
		double gamma_fac, double kheat_fac, int mm, double phi0, double a,
		double b, double c, Vector x0, int inside) {
	int il = 0;
	int ih = r->N - npert;
	Vector CoM = compute_com(r, il, ih);
	int NC = 0;
	for (int i = 0; i < NS; i++) {
		// compute spring mid point from central position
		Vector x_mid = spr_mid(r, springs[i], CoM);
		Vector x = x_mid - x0;
		double rmid2 = pow(x.getX() / a, 2.0)
					 + pow(x.getY() / b, 2.0)
					 + pow(x.getZ() / c, 2.0);
		double phi = atan2(x.getY(), x.getX());

		double angfac = cos(mm * (phi - phi0));
		if ((rmid2 <= 1.0) && (inside == 1)) {
			springs[i].ks *= 1.0 + ks_fac * angfac;
			springs[i].gamma *= 1.0 + gamma_fac * angfac;
			springs[i].k_heat *= 1.0 + kheat_fac * angfac;
			NC++;
		}
		if ((rmid2 >= 1.0) && (inside == 0)) {
			springs[i].ks *= 1.0 + ks_fac * angfac;
			springs[i].gamma *= 1.0 + gamma_fac * angfac;
			springs[i].k_heat *= 1.0 + kheat_fac * angfac;
			NC++;
		}
	}
	printf("adjust_ks_abc_fac: number of springs changed %d\n", NC);
}

// change all mass nodes with r within or without ellipsoid 
// with (x-x0)^2/a^2 + (y-y0)^2/b^2 + (z-z0)^2/c^2 = 1
// by factor mfac (multiplied)
// done with respect to center of mass
// then normalize so that total mass sum to 1 again
// mfac is the density ratio if nodes are evenly distributed in space  
// inside==1 then inside ellipsoid
// inside==0 then outside
void adjust_mass_abc(struct reb_simulation *const r, int npert, double mfac,
		double a, double b, double c, Vector x0, int inside) {
	struct reb_particle *particles = r->particles;
	int il = 0;
	int ih = r->N - npert;
	int ncore = 0;
	Vector CoM = compute_com(r, il, ih);
	double rfac = pow(mfac, 1.0 / 2.0); // interesting choice!!!!
	for (int i = il; i < ih; i++) {
		Vector x = {particles[i].x, particles[i].y, particles[i].z};
		x -= CoM;
		x -= x0;
		double rmid2 = pow(x.getX() / a, 2.0)
					 + pow(x.getY() / b, 2.0)
					 + pow(x.getZ() / c, 2.0);
		if (inside == 1) {
			if (rmid2 < 1.0) {
				ncore++; // numbers of nodes in core
				particles[i].m *= mfac;
				particles[i].r *= rfac;
			} else {
				particles[i].r /= rfac;
			}
		}
		if (inside == 0) {
			if (rmid2 > 1.0) {
				particles[i].m *= mfac;
				particles[i].r *= rfac;
			} else {
				particles[i].r /= rfac;
				ncore++; // numbers of nodes in core
			}
		}
	}
	// if nodes are evenly distributed in volume then mass of node scales with density

	printf("adjust_mass_abc: ncore=%d nshell=%d\n", ncore,
			r->N - npert - ncore);
	double tm = 0.0; // total mass
	for (int i = il; i < ih; i++)
		tm += particles[i].m; // total mass!
	for (int i = il; i < ih; i++)
		particles[i].m /= tm;
	// now normalized mass to 1
}

// change all masses with x>xmin by factor mfac, then rescale so sum is still 1
// done with respect to origin
void adjust_mass_side(struct reb_simulation *const r, int npert, double mfac,
		double xmin) {
	struct reb_particle *particles = r->particles;
	int il = 0;
	int ih = r->N - npert;
	double tm = 0.0;
	int npart = 0;
	for (int i = il; i < ih; i++) {
		if (particles[i].x > xmin) {
			particles[i].m *= mfac;
			npart++;
		}
	}
	for (int i = il; i < ih; i++)
		tm += particles[i].m;
	for (int i = il; i < ih; i++)
		particles[i].m /= tm;
	printf("npart=%d \n", npart);

}

// returns 0 if even 1 if odd
int mym(int k) {
	int z = floor(k / 2);
	return abs(k - 2 * z); // returns 0 or 1
}

// make a football of particles hcp lattice
// and total mass = total_mass
// hexagonal close packed setup
// dd is minimum distance between particles
// ax,by,cz are semi-axes
double fill_hcp(struct reb_simulation *r, double dd, double ax, double by,
		double cz, double total_mass) {
	// struct reb_particle* particles = r->particles;
	struct reb_particle pt;
	int i0 = r->N; // store initial particle number
	double dx = 1.08 * dd;
	int nx = (int) (1.2 * ax / dx);
	// double dx = 2.0*ax/nx;
	double zfac = dx * sqrt(2.0 / 3.0); // factors for hcp lattice
	double yfac = dx * sin(M_PI / 3.0);
	double xfac_s = dx / 2.0; // shifting
	// printf("dx =%.3f zfac=%.3f yfac=%.3f\n",dx,zfac,yfac);
	// fill a cube
	double yy = 1.2 * by / yfac;
	int ny = (int) yy;
	double zz = 1.2 * cz / zfac;
	int nz = (int) zz;
	double midvalx = 0.0 * nx * dx; // center of ball
	double midvaly = 0.0 * ny * yfac;
	double midvalz = 0.0 * nz * zfac;
	double particle_radius = dx / 2.0; // temporary set
	// make an hcp grid
	for (int k = -nz; k <= nz; k++) {
		double z = zfac * k;
		for (int j = -ny; j <= ny; j++) {
			double y = yfac * (j + 0.5 * mym(k));
			for (int i = -nx; i <= nx; i++) {
				double x = dx * i + xfac_s * mym(j) + xfac_s * mym(k);
				// printf("%.2f %.2f %.2f\n",x-midval,y-midval,z-midval);
				pt.m = 1.0;
				pt.x = x - midvalx;
				pt.y = y - midvaly;
				pt.z = z - midvalz;
				pt.vx = 0;
				pt.ax = 0;
				pt.vy = 0;
				pt.ay = 0;
				pt.vz = 0;
				pt.az = 0;
				pt.r = particle_radius;
				double xa2 = pow(pt.x / ax, 2.0);
				double yb2 = pow(pt.y / by, 2.0);
				double zc2 = pow(pt.z / cz, 2.0);
				double rval = sqrt(xa2 + yb2 + zc2);
				if (rval <= 1.0) {
					reb_add(r, pt);
				}
			}
		}
	}
	int N = r->N;
	// printf("i0=%d N=%d\n",i0,N);
	double particle_mass = total_mass / (N - i0);
	double min_d = mindist(r, i0, N);
	// correct the radii and mass of the particles

	for (int i = i0; i < N; i++) {
		r->particles[i].m = particle_mass;
		r->particles[i].r = min_d / 4.0;
	}

	printf("fill_hcp: Nparticles =%d dx=%.2e dr=%.2e rad=%.2e m=%.2e\n", N - i0,
			dx, min_d, min_d / 2.0, particle_mass);
	return min_d;
}

// make a football of particles cubic lattice
// and total mass = total_mass
// dd is minimum distance between particles (cube side)
// ax,by,cz are semi-axes
double fill_cubic(struct reb_simulation *r, double dd, double ax, double by,
		double cz, double total_mass) {
	struct reb_particle pt;
	int i0 = r->N; // store initial particle number
	int nz = (int) (1.2 * ax / dd);
	int ny = nz;
	int nx = nz;
	double zfac = dd;
	for (int k = -nz; k <= nz; k++) {
		double z = zfac * k;
		for (int j = -ny; j <= ny; j++) {
			double y = zfac * j;
			for (int i = -nx; i <= nx; i++) {
				double x = zfac * i;
				pt.m = 1.0;
				pt.x = x;
				pt.y = y;
				pt.z = z;
				pt.vx = 0;
				pt.ax = 0;
				pt.vy = 0;
				pt.ay = 0;
				pt.vz = 0;
				pt.az = 0;
				pt.r = 1.0;
				double xa2 = pow(pt.x / ax, 2.0);
				double yb2 = pow(pt.y / by, 2.0);
				double zc2 = pow(pt.z / cz, 2.0);
				double rval = sqrt(xa2 + yb2 + zc2);
				if (rval <= 1.0) {
					reb_add(r, pt);
				}
			}
		}
	}
	int N = r->N;
	double particle_mass = total_mass / (N - i0);
	double min_d = mindist(r, i0, N);
	// correct the radii and mass of the particles
	for (int i = i0; i < N; i++) {
		r->particles[i].m = particle_mass;
		r->particles[i].r = min_d / 4.0;
	}
	printf("fill_cubic: Nparticles =%d dx=%.2e rad=%.2e m=%.2e\n", N - i0, dd,
			min_d / 2.0, particle_mass);
	return min_d;
}

// shift the position of a resolved body by dx,dy,dz,dvx,dvy,dvz
// shift all particles positions and velocities in the body [il,ih)
void move_resolved(struct reb_simulation *r, double dx, double dy, double dz,
		double dvx, double dvy, double dvz, int il, int ih) {
	// struct reb_particle* particles = r->particles;
	// struct reb_particle pt;
	for (int i = il; i < ih; i++) {
		r->particles[i].x += dx;
		r->particles[i].y += dy;
		r->particles[i].z += dz;
		r->particles[i].vx += dvx;
		r->particles[i].vy += dvy;
		r->particles[i].vz += dvz;
	}
}

// read in a vertex (shape) file, list of vertices
// the file is filename
// file is ascii in form "v x y z"  were v is ascii v if a vertex
// and xyz are vertex positions units km
void read_vertex_file(struct reb_simulation *r, char *filename) {
	// struct reb_particle* particles = r->particles;
	int il = r->N;
	struct reb_particle pt;
	FILE *fp;
	fp = fopen(filename, "r");
	char ss[300];
	pt.ax = 0.0;
	pt.ay = 0.0;
	pt.az = 0.0;
	pt.vx = 0.0;
	pt.vy = 0.0;
	pt.vz = 0.0;
	pt.m = 1.0; // arbitrary
	pt.r = 1.0; // arbitrary
	while (fgets(ss, 300, fp) != NULL) {
		double x, y, z;
		char vf;
		sscanf(ss, "%c %lf %lf %lf", &vf, &x, &y, &z);
		if (vf == 'v') {  // is a vertex not a face
			// printf("%s\n",ss);
			pt.x = x;
			pt.y = y;
			pt.z = z;
			reb_add(r, pt);
		}
	}
	fclose(fp);
	int ih = r->N;
	printf("read_vertex_file: number of vertices read in %d\n", ih - il);
	double dd = mindist(r, il, ih);
	for (int i = il; i < ih; i++) { // adjust radius of each vertex particle node
		r->particles[i].r = dd / 4.0;
	}
	// return dd;
}

// return maximum radius of a particle distribution
double max_radius(struct reb_simulation *r, int il, int ih) {
	double max_r2 = 0.0;
	for (int i = il; i < ih; i++) {
		double r2 = pow(r->particles[i].x, 2.0);
		r2 += pow(r->particles[i].y, 2.0);
		r2 += pow(r->particles[i].z, 2.0);
		if (r2 > max_r2)
			max_r2 = r2;
	}
	return sqrt(max_r2);
}

// return minimum radius of a particle distribution
double min_radius(struct reb_simulation *r, int il, int ih) {
	double min_r2 = 0.0;
	for (int i = il; i < ih; i++) {
		double r2 = pow(r->particles[i].x, 2.0);
		r2 += pow(r->particles[i].y, 2.0);
		r2 += pow(r->particles[i].z, 2.0);
		if (r2 < min_r2)
			min_r2 = r2;
	}
	return sqrt(min_r2);
}

// return index of closest particle to position x,y,z, within [il,ih)
// before doing this routine check to see if radius is within largest to smallest
// of shape model
// this just finds nearest particle, can be used outside of shape model stuff
// note, we have not necessarily centered the body prior to calling this routine
// origin might not be center of body
int nearest_to_shape(struct reb_simulation *r, int il, int ih, double x,
		double y, double z) {
	double dist2 = 1e30;
	int i0 = 0;
	// printf("%.3f %.3f %.3f \n",x,y,z);
	for (int i = il; i < ih; i++) {
		double r2 = pow(r->particles[i].x - x, 2.0);
		r2 += pow(r->particles[i].y - y, 2.0);
		r2 += pow(r->particles[i].z - z, 2.0);
		if (r2 < dist2) {
			i0 = i;
			dist2 = r2;
			// printf("%d %.3f\n",i0,dist2);
		}
		// printf("%.3f %.3f %.3f %.3e\n",r->particles[i].x, r->particles[i].y,r->particles[i].z,r2);
	}
	// printf("final %d %.3f\n",i0, r->particles[i0].z);
	// exit(0);
	return i0; // ACQ fix xxxxxx
}

// is position x,y,z within the shape given by vertices [il,ih-1]?
// minr, maxr are the maximum and mininum radii of shape vertices, precomputed
int within_shape(struct reb_simulation *r, int il, int ih, double minr,
		double maxr, double x, double y, double z) {
	double r2 = pow(x, 2.0) + pow(y, 2.0) + pow(z, 2.0);
	double rad = sqrt(r2);
	if (rad < minr)
		return 1;  // is within shape
	if (rad > maxr)
		return 0;  // is outside shape
	// find index of nearest particle in shape model
	int j0 = nearest_to_shape(r, il, ih, x, y, z);
	// compute radius of this particle in shape model
	r2 = pow(r->particles[j0].x, 2.0);
	r2 += pow(r->particles[j0].y, 2.0);
	r2 += pow(r->particles[j0].z, 2.0);
	double rs = sqrt(r2);
	if (rad > rs)
		return 0; // outside shape locally
	else
		return 1;  // inside shape locally

}

// remove the shape vertices
void rmshape_vertices(struct reb_simulation *const r, int N_bennu) {
	for (int i = 0; i < N_bennu; i++) {
		reb_remove(r, 0, 1); // last 0 is without keeping order, last 1 shifts particles
		// here we keep removing first particle in list
	}
}

// create a ~uniform random particle distribution with total mass  = total_mass
// fill particles within a bennu's shape given by a set of particles 
//  already read in 
// spacing set by parameter dist: no closer particles allowed
void rand_bennu(struct reb_simulation *const r, double dist,
		double total_mass) {
	struct reb_particle pt;
	// we assume that bennu shape model vertices have already been read in
	// and they are currently [0,N-1]
	int N_bennu = r->N;
	double min_r_shape = min_radius(r, 0, N_bennu); // min and max radii of shape model
	double max_r_shape = max_radius(r, 0, N_bennu); // particles
	int npart = 40 * pow(2.0 * max_r_shape / dist, 3.0);
	// guess for number of random particles we need to generate
	pt.ax = 0.0;
	pt.ay = 0.0;
	pt.az = 0.0;
	pt.vx = 0.0;
	pt.vy = 0.0;
	pt.vz = 0.0;
	pt.m = 1.0;
	pt.r = dist / 4.0;
	int N = r->N;
	int i0 = N; // newly added particles start at i0
	for (int i = 0; i < npart; i++) { // generate particles
		pt.x = reb_random_uniform(-max_r_shape, max_r_shape);
		pt.y = reb_random_uniform(-max_r_shape, max_r_shape);
		pt.z = reb_random_uniform(-max_r_shape, max_r_shape);
		int ws = within_shape(r, 0, N_bennu, min_r_shape, max_r_shape, pt.x,
				pt.y, pt.z);
		if (ws == 1.0) { // within shape
			int toonear = 0;  // is there a particle too nearby?
			int j = i0;
			N = r->N;
			while ((toonear == 0) && (j < N)) {
				double dx = pt.x - r->particles[j].x;
				double dy = pt.y - r->particles[j].y;
				double dz = pt.z - r->particles[j].z;
				double dr = sqrt(dx * dx + dy * dy + dz * dz);
				if (dr < dist)
					toonear = 1;
				j++;
			}
			if (toonear == 0)
				reb_add(r, pt);
			// only add particle if not near any other
		}
	}
	N = r->N;

	// adjust mass of each particle so that sums to desired total mass
	double particle_mass = total_mass / (N - i0);
	// fix masses
	for (int ii = i0; ii < N; ii++) { // all particles!
		r->particles[ii].m = particle_mass;
	}
	double md = mindist(r, i0, N);
	printf("rand_bennu: Nparticles=%d min_d=%.2f\n", N - i0, md);

	// remove the shape model vertices
	// rmshape_vertices(r,  N_bennu);
	// printf("rand_bennu: shape removed \n");
}

void rescale_xyz(struct reb_simulation *const r, int il, int ih,
		double scalefac) {
	for (int i = il; i < ih; i++) {
		r->particles[i].x *= scalefac;
		r->particles[i].y *= scalefac;
		r->particles[i].z *= scalefac;
	}
}

// using shape model vertices
// create and fill an array that marks particles that are near the surface
// array has a 0 if particle is not near the surface
// otherwise is 1
// index starts from 0 first particle not in shape model
// mind is minimum distance from a surface particle
int* marksurface(struct reb_simulation *const r, int N_bennu, double mind) {
	int npsurf = 0; // number of particles in surface
	int *surfp;
	surfp = (int*) malloc(sizeof(int) * r->N);
	for (int j = N_bennu; j < r->N; j++) {
		int i0 = j - N_bennu;
		surfp[i0] = 0; // index of particle once shape model removed
		double x = r->particles[j].x; // position of particle not in shape model
		double y = r->particles[j].y;
		double z = r->particles[j].z;
		int j0 = nearest_to_shape(r, 0, N_bennu, x, y, z); // j0 tells index of
		// nearest particle in shape model
		double xj = r->particles[j0].x;
		double yj = r->particles[j0].y;
		double zj = r->particles[j0].z;
		double r2 = pow(x - xj, 2.0) + pow(y - yj, 2.0) + pow(z - zj, 2.0);
		double rdist = sqrt(r2); // distance to nearest point in shape model

		if (rdist < mind) {
			surfp[i0] = 1;  // it is near the surface
			npsurf++;
		} else {
			r->particles[j].r = 0.001; // shrink so hard to see!
			// leaving visible only surface particles
		}
	}
	printf("npsurf_bennu = %d\n", npsurf);
	return surfp;
}

// mark particles near the surface of a rand_football
int* marksurface_football(struct reb_simulation *const r, double mind,
		double ax, double by, double cz) {
	int npsurf = 0; // number of particles in surface
	int *surfp;
	surfp = (int*) malloc(sizeof(int) * r->N);
	for (int j = 0; j < r->N; j++)
		surfp[j] = 0.0;
	for (int j = 0; j < r->N; j++) {
		double x = r->particles[j].x;  //
		double y = r->particles[j].y;
		double z = r->particles[j].z;
		double xa2 = pow(x / ax, 2.0);
		double ya2 = pow(y / by, 2.0);
		double za2 = pow(z / cz, 2.0);
		double rval = sqrt(xa2 + ya2 + za2);
		double rdist = 1.0 - rval;
		if ((rdist < mind) && (rdist > 0)) {
			surfp[j] = 1;  // it is near the surface
			npsurf++;
		} else {
			// r->particles[j].r = 0.001; // shrink so hard to see!
			// leaving visible only surface particles
		}
	}
	printf("npsurf_football = %d\n", npsurf);
	return surfp;
}

// mark particles near the surface of a rand_cone    
int* marksurface_cone(struct reb_simulation *const r, double mind, double rb,
		double slope) {
	int npsurf = 0; // number of particles in surface
	int *surfp;
	surfp = (int*) malloc(sizeof(int) * r->N);
	double ifac = sqrt(1.0 + slope * slope);
	double hcone = rb * slope;
	for (int j = 0; j < r->N; j++)
		surfp[j] = 0.0;
	for (int j = 0; j < r->N; j++) {
		double x = r->particles[j].x;  //
		double y = r->particles[j].y;
		double z = r->particles[j].z;
		double radius = sqrt(x * x + y * y);
		double dplus = fabs(z + slope * radius - hcone) / ifac; // distances to lines
		double dminus = fabs(z - slope * radius + hcone) / ifac;
		if ((dplus < mind) || (dminus < mind)) {
			surfp[j] = 1;  // it is near the surface
			npsurf++;
		} else {
			r->particles[j].r = 0.001; // shrink so hard to see!
			// leaving visible only surface particles
		}
	}
	printf("npsurf_cone = %d\n", npsurf);
	return surfp;
}

// ZYH vector to longitude and latitude
// input v[3] = (x,y,z)  output a[2] = (long, lat) in radians
// long in [0,2pi]
// lat in [0,pi] note! could change
void potoang(double v[3], double a[2]) {
	double x = v[0];
	double y = v[1];
	double z = v[2];
	double sinlat = z / sqrt(x * x + y * y + z * z); // z/r  range [-1,1]
	a[1] = asin(sinlat);  // range [-pi/2 to pi/2]?  latitude
	// double w = sqrt(x*x+y*y);
	a[0] = atan2(y, x);  // range [-pi/pi]  longitude
	if (a[0] < 0.0)
		a[0] += 2.0 * M_PI;  // range [0,2pi] longitude
	a[1] += M_PI / 2;  // returns [0,pi] note offset!!!!  range [0,pi]
}

// connect springs to all particles with interparticle
// distances less than h_dist_max apart
// distances greater than h_dist_min apart
// for particle index range [i0,imax-1]
// spring added with rest length at current length
// nodemax is the maximum number of connections a single node can have
// don't allow springs to be connected if the node they connect to has more than nodemax 
// previously connected springs
void connect_springs_dist_nodemax(struct reb_simulation *const r,
		double h_dist_min, double h_dist_max, int i0, int imax,
		struct spring spring_vals, int nodemax) {
	if (imax <= i0)
		return;
	int *n_nodes;
	n_nodes = (int*) malloc(r->N * sizeof(int));
	for (int ii = 0; ii < r->N; ii++)
		n_nodes[ii] = 0;
	for (int k = 0; k < NS; k++) {
		int ii = springs[k].i;
		int jj = springs[k].j;
		n_nodes[ii]++;
		n_nodes[jj]++;
	}
	// find all the springs for near neighbors
	for (int ii = i0; ii < imax - 1; ii++) {
		// printf("hello %d \n",ii);
		double xi = r->particles[ii].x;
		double yi = r->particles[ii].y;
		double zi = r->particles[ii].z;
		for (int jj = ii + 1; jj < imax; jj++) { // all pairs
			double xj = r->particles[jj].x;
			double yj = r->particles[jj].y;
			double zj = r->particles[jj].z;
			double dr = sqrt(
					(xi - xj) * (xi - xj) + (yi - yj) * (yi - yj)
							+ (zi - zj) * (zi - zj));
			if ((dr < h_dist_max) && (dr > h_dist_min)) { // try to add the spring
				if ((n_nodes[ii] < nodemax) && (n_nodes[jj] < nodemax)) {
					int ns_temp = NS;
					add_spring(r, ii, jj, spring_vals); // will not be added if there is already
					// one there
					// spring added at rest distance
					if (NS > ns_temp) {  // a spring was added
						n_nodes[ii]++;
						n_nodes[jj]++;
					}
				}

			}
		}
	}
	printf("add springs nodemax: NS=%d\n", NS);
	free(n_nodes);
}

// return the spring force vector for spring i on particle springs[i].i
// also return length vector between particle i and j as Lx,Ly,Lz
// Lx,Ly,Lz = r_ij
// ignore damping
void spring_force_k_one(struct reb_simulation *const r, int i_s, double *Fx,
		double *Fy, double *Fz, double *Lx, double *Ly, double *Lz) {

	double L = spring_length(r, springs[i_s]) + L_EPS
	; // spring length
	double rs0 = springs[i_s].rs0; // rest length
	int ii = springs[i_s].i;
	int jj = springs[i_s].j;
	double dx = r->particles[ii].x - r->particles[jj].x;
	double dy = r->particles[ii].y - r->particles[jj].y;
	double dz = r->particles[ii].z - r->particles[jj].z;
	// double mii = r->particles[ii].m;
	// double mjj = r->particles[jj].m;
	double ks = springs[i_s].ks;
	double fac = -ks * (L - rs0) / L; // L here to normalize direction
	// accelerations are force divided by mass
	*Fx = fac * dx;
	*Fy = fac * dy;
	*Fz = fac * dz;
	*Lx = dx;
	*Ly = dy;
	*Lz = dz;
}

void zero_symtensor(struct symtensor *S) {
	S->xx = 0.0;
	S->yy = 0.0;
	S->zz = 0.0;
	S->xy = 0.0;
	S->yz = 0.0;
	S->xz = 0.0;
}

void stretch(struct reb_simulation *const r, int il, int ih, double scalex,
		double scaley, double scalez) {
	for (int i = il; i < ih; i++) {
		r->particles[i].x *= scalex;
		r->particles[i].y *= scaley;
		r->particles[i].z *= scalez;
	}
}

// rotate the body so that z is along biggest eigenvalue, x along smallest!
// seems to work!
void rotate_to_principal(struct reb_simulation *const r, int il, int ih) {
	Matrix inertia_mat = mom_inertia(r, il, ih); // get moment of inertial matrix
	std::cout.precision(3);
	std::cout << "Moment of inertia: " << inertia_mat << std::endl;
	double eigs[3];
	eigenvalues(inertia_mat, eigs);
	printf("I eigenvalues %.3f %.3f %.3f\n", eigs[0], eigs[1], eigs[2]);

	Vector eigvec1 = eigenvector(inertia_mat, eigs[2]); // note order eig3 here
	printf("evec1 %.3f %.3f %.3f\n", eigvec1.getX(), eigvec1.getY(), eigvec1.getZ());
	Vector eigvec2 = eigenvector(inertia_mat, eigs[1]);
	printf("evec2 %.3f %.3f %.3f\n", eigvec2.getX(), eigvec2.getY(), eigvec2.getZ());
	Vector eigvec3 = eigenvector(inertia_mat, eigs[0]); // eig1 here
	printf("evec3 %.3f %.3f %.3f\n", eigvec3.getX(), eigvec3.getY(), eigvec3.getZ());

	double ex1 = eigvec1.getX(), ex2 = eigvec2.getX(), ex3 = eigvec3.getX(), ey1 =
			eigvec1.getY(), ey2 = eigvec2.getY(), ey3 = eigvec3.getY(), ez1 = eigvec1.getZ(),
			ez2 = eigvec2.getZ(), ez3 = eigvec3.getZ();

	struct reb_particle *particles = r->particles;
	for (int i = il; i < ih; i++) {
		double xr, yr, zr, vxr, vyr, vzr;
		double x0 = particles[i].x;
		double y0 = particles[i].y;
		double z0 = particles[i].z;
		vec_mul(ex1, ey1, ez1, // rotate coordinates
				ex2, ey2, ez2, ex3, ey3, ez3, x0, y0, z0, &xr, &yr, &zr);
		particles[i].x = xr;
		particles[i].y = yr;
		particles[i].z = zr;
		double vx0 = particles[i].vx;
		double vy0 = particles[i].vy;
		double vz0 = particles[i].vz;
		vec_mul(ex1, ey1, ez1, // rotate velocities too
				ex2, ey2, ez2, ex3, ey3, ez3, vx0, vy0, vz0, &vxr, &vyr, &vzr);
		particles[i].vx = vxr;
		particles[i].vy = vyr;
		particles[i].vz = vzr;
	}

	inertia_mat = mom_inertia(r, il, ih); // get moment of inertial matrix
	std::cout.precision(3);
	std::cout << "I matrix: \n" << inertia_mat << std::endl;
	eigenvalues(inertia_mat, eigs);
	std::cout << "I eigenvalues: " << eigs[0] << ", " << eigs[1] << ", " << eigs[2] << std::endl;
}

// multiply a vector by a matrix, helper routine for above
void vec_mul(double Axx, double Axy, double Axz, double Ayx, double Ayy,
		double Ayz, double Azx, double Azy, double Azz, double bx, double by,
		double bz, double *cx, double *cy, double *cz) {
	*cx = Axx * bx + Axy * by + Axz * bz;
	*cy = Ayx * bx + Ayy * by + Ayz * bz;
	*cz = Azx * bx + Azy * by + Azz * bz;
}

