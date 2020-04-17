/*
 * shapes.cpp
 *
 * Create and modify distributions of particles
 *
 *  Created on: Apr 3, 2020
 *      Author: alex
 */

#include <cfloat>
#include <cmath>
extern "C" {
#include "rebound.h"
}
#include "springs.h"
#include "shapes.h"

/**************/
/* Generators */
/**************/

// Create an approximately uniform random particle distribution with total_mass within rectangular prism given by lengths x, y, z
// No particles closer than min_dist created
void rand_rectangle(struct reb_simulation *const n_body_sim, double min_dist,
		double x, double y, double z, double total_mass) {
	// Get particle info
	struct reb_particle *particles = n_body_sim->particles;

	// Guess at number of random particles we need to generate
	int n_part = 40.0 * x * y * z / pow(min_dist, 3.0);
	std::cout << "rand_rectangle: Trying to create " << n_part << "particles."
			<< std::endl;

	// Set up particle defaults
	struct reb_particle pt;
	pt.ax = 0.0;
	pt.ay = 0.0;
	pt.az = 0.0;
	pt.vx = 0.0;
	pt.vy = 0.0;
	pt.vz = 0.0;
	pt.m = 1.0;
	double particle_radius = min_dist / 2.0;
	pt.r = particle_radius / 3.0;  // temporary

	// Get current number of particles
	int N = n_body_sim->N;
	int i_0 = N;

	// Get uniform distribution inside rectangle
	for (int i = 0; i < n_part; i++) {
		pt.x = reb_random_uniform(-x / 2, x / 2);
		pt.y = reb_random_uniform(-y / 2, y / 2);
		pt.z = reb_random_uniform(-z / 2, z / 2);

		// Is there a particle too nearby?
		bool too_near = false;
		N = n_body_sim->N;
		for (int j = i_0; not too_near && j < N; j++) {
			Vector dx = { pt.x - particles[j].x, pt.y - particles[j].y, pt.z
					- particles[j].z };
			if (dx.len() < min_dist)
				too_near = true;
		}

		// Only add particle if not near any other
		if (not too_near)
			reb_add(n_body_sim, pt);
	}

	// Get current number of particles
	N = n_body_sim->N;

	// Adjust mass of each particle so that they sum to desired total mass
	double particle_mass = total_mass / (N - i_0);
	for (int i = i_0; i < N; i++) {
		particles[i].m = particle_mass;
	}

	// Double-check min_dist in created particles
	double min_d = mindist(n_body_sim, i_0, N);
	std::cout << "rand_rectangle: Created " << N - i_0
			<< " particles separated by minimum distance " << min_d
			<< std::endl;
}

// Create an approximately uniform random particle distribution with total_mass within rectangle (2D) given by lengths x, y
// No particles closer than min_dist created
void rand_rectangle_2d(struct reb_simulation *const n_body_sim, double min_dist,
		double x, double y, double total_mass) {
	// Get particle info
	reb_particle *particles = n_body_sim->particles;

	// Guess at number of particles to create
	int n_part = 40.0 * x * y / pow(min_dist, 2.0);
	printf("rand_rectangle_2d: n_part=%d\n", n_part);

	// Set particle defaults
	struct reb_particle pt;
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
	double particle_radius = min_dist / 2.0;
	pt.r = particle_radius / 3.0;  // temporary

	// Get current number of particles
	int N = n_body_sim->N;
	int i_0 = N;

	// Get uniform distribution inside 2D rectangle
	for (int i = 0; i < n_part; i++) {
		pt.x = reb_random_uniform(-x / 2, x / 2);
		pt.y = reb_random_uniform(-y / 2, y / 2);

		// Is there a particle too nearby?
		bool too_near = false;
		N = n_body_sim->N;
		for (int j = i_0; not too_near && j < N; j++) {
			Vector dx = { pt.x - particles[j].x, pt.y - particles[j].y, pt.z
					- particles[j].z };
			if (dx.len() < min_dist)
				too_near = true;
		}

		// Only add particle if not near any other
		if (not too_near)
			reb_add(n_body_sim, pt);
	}

	// Get current number of particles
	N = n_body_sim->N;

	// Adjust mass of each particle so that they sum to desired total mass
	double particle_mass = total_mass / (N - i_0);
	for (int i = i_0; i < N; i++) {
		particles[i].m = particle_mass;
	}

	// Double check min_dist
	double min_d = mindist(n_body_sim, i_0, N);
	std::cout << "rand_rectangle: Created " << N - i_0
			<< " particles separated by minimum distance " << min_d
			<< std::endl;
}

// Create an approximately uniform random particle distribution with total_mass within a cone given by base radius, height
// No particles closer than min_dist created
void rand_cone(struct reb_simulation *const n_body_sim, double min_dist,
		double radius, double height, double total_mass) {
	// Get particle info
	reb_particle *particles = n_body_sim->particles;

	// Calculate slope of sides of cone
	double slope = height / radius;

	// Guess at number of random particles we need to generate
	int n_part = 40 * pow(2.0 * radius / min_dist, 3.0);

	// Set defaults for particles
	struct reb_particle pt;
	pt.ax = 0.0;
	pt.ay = 0.0;
	pt.az = 0.0;
	pt.vx = 0.0;
	pt.vy = 0.0;
	pt.vz = 0.0;
	pt.m = 1.0;
	double particle_radius = min_dist / 2.0;
	pt.r = particle_radius / 2.0;  // temporary

	// Get current number of particles
	int N = n_body_sim->N;
	int i_0 = N;

	// Get uniform distribution of particles in cone
	for (int i = 0; i < n_part; i++) {
		pt.x = reb_random_uniform(-radius, radius);
		pt.y = reb_random_uniform(-radius, radius);
		pt.z = reb_random_uniform(0, height);

		// Get current location, maximum height
		double cur_rad = sqrt(pt.x * pt.x + pt.y * pt.y);
		double max_height = height - cur_rad * slope; // h - slope*r

		// Check if created particle is inside cone
		if ((cur_rad < radius) && (pt.z < max_height)) {
			// Is there a particle too nearby?
			bool too_near = false;
			N = n_body_sim->N;
			for (int j = i_0; not too_near && j < N; j++) {
				Vector dx = { pt.x - particles[j].x, pt.y - particles[j].y, pt.z
						- particles[j].z };
				if (dx.len() < min_dist)
					too_near = true;
			}

			// Only add particle if not near any other
			if (not too_near)
				reb_add(n_body_sim, pt);
		}
	}

	// Get current number of particles
	N = n_body_sim->N;

	// Adjust mass of each particle so that sums to desired total mass
	double particle_mass = total_mass / (N - i_0);
	for (int i = i_0; i < N; i++) {
		particles[i].m = particle_mass;
	}

	// Double check min_dist
	double min_d = mindist(n_body_sim, i_0, N);
	std::cout << "rand_cone: Created " << N - i_0
			<< " particles separated by minimum distance " << min_d
			<< std::endl;
}

// Create an approximately uniform random particle distribution with total_mass within an ellipsoid shape given by semi-axes x, y, z
// No particles closer than min_dist created
void rand_ellipsoid(struct reb_simulation *const n_body_sim, double min_dist,
		double x, double y, double z, double total_mass) {
	// Get particle info
	reb_particle *particles = n_body_sim->particles;

	// Guess at number of random particles we need to generate
	int n_part = 40 * pow(2.0 * x / min_dist, 3.0);
	std::cout << "rand_ellipsoid: Creating " << n_part << " particles."
			<< std::endl;

	// Set particle defaults
	struct reb_particle pt;
	pt.ax = 0.0;
	pt.ay = 0.0;
	pt.az = 0.0;
	pt.vx = 0.0;
	pt.vy = 0.0;
	pt.vz = 0.0;
	pt.m = 1.0;
	double particle_radius = min_dist / 2.0;
	pt.r = particle_radius / 2.0;  // temporary

	// Get current number of particles
	int N = n_body_sim->N;
	int i_0 = N;

	// Get uniform distribution inside ellipsoid
	for (int i = 0; i < n_part; i++) {
		pt.x = reb_random_uniform(-x, x);
		pt.y = reb_random_uniform(-y, y);
		pt.z = reb_random_uniform(-z, z);
		Vector pos = { pt.x / x, pt.y / y, pt.z / z };

		// Check if particle is inside ellipsoid
		if (pos.len() < 1.0) {

			// Is there a particle too nearby?
			bool too_near = false;
			N = n_body_sim->N;
			for (int j = i_0; not too_near && j < N; j++) {
				Vector dx = { pt.x - particles[j].x, pt.y - particles[j].y, pt.z
						- particles[j].z };
				if (dx.len() < min_dist)
					too_near = true;
			}

			// Only add particle if not near any other
			if (not too_near)
				reb_add(n_body_sim, pt);
		}
	}

	// Get current number of particles
	N = n_body_sim->N;

	// Adjust mass of each particle so that they sum to desired total mass
	double particle_mass = total_mass / (N - i_0);
	for (int i = i_0; i < N; i++) {
		n_body_sim->particles[i].m = particle_mass;
	}

	// Double check min_dist
	double min_d = mindist(n_body_sim, i_0, N);
	std::cout << "rand_ellipsoid: Created " << N - i_0
			<< " particles separated by minimum distance " << min_d
			<< std::endl;
}

// Create an ellipsoid of particles using hexagonal close-packed lattice with total_mass, semi-axes x, y, z
// No particles closer than min_dist created
void hcp_ellipsoid(struct reb_simulation *n_body_sim, double min_dist, double x,
		double y, double z, double total_mass) {
	// Get particle info
	struct reb_particle *particles = n_body_sim->particles;

	// Set particle defaults
	struct reb_particle pt;
	pt.vx = 0;
	pt.ax = 0;
	pt.vy = 0;
	pt.ay = 0;
	pt.vz = 0;
	pt.az = 0;
	pt.m = 1.0;

	// Get current number of particles
	int N = n_body_sim->N;
	int i_0 = N;

	// Lattice spacing
	double dx = 1.08 * min_dist;
	pt.r = dx / 2.0; // temporary

	// Lattice spacing factors (cf hexagons)
	double xfac = dx / 2.0;
	double yfac = dx * sin(M_PI / 3.0);
	double zfac = dx * sqrt(2.0 / 3.0);

	// Number of points in each direction
	int nx = (int) (1.2 * x / dx);
	int ny = (int) (1.2 * y / yfac);
	int nz = (int) (1.2 * z / zfac);

	// Make an HCP grid
	for (int k = -nz; k <= nz; k++) {

		// z location
		double cur_z = zfac * k;

		for (int j = -ny; j <= ny; j++) {

			// y location
			double cur_y = yfac * (j + 0.5 * (k % 2));

			for (int i = -nx; i <= nx; i++) {

				// x location
				double cur_x = dx * i + xfac * (j % 2) + xfac * (k % 2);

				// Set particle location
				pt.x = cur_x;
				pt.y = cur_y;
				pt.z = cur_z;

				// Add if particle is inside ellipsoid
				Vector pos = { pt.x / x, pt.y / y, pt.z / z };
				if (pos.len() <= 1.0)
					reb_add(n_body_sim, pt);
			}
		}
	}

	// Get current number of particles
	N = n_body_sim->N;

	// Adjust mass of each particle so that they sum to desired total mass
	// Adjust radius of each particle so they have radius of 1/4 smallest distance between any particles
	double particle_mass = total_mass / (N - i_0);
	double min_d = mindist(n_body_sim, i_0, N);
	for (int i = i_0; i < N; i++) {
		particles[i].m = particle_mass;
		particles[i].r = min_d / 4.0;
	}

	std::cout << "hcp_ellipsoid: Created " << N - i_0
			<< " particles with lattice spacing " << dx
			<< ", separated by minimum distance " << min_d << std::endl;
}

// Create an ellipsoid of particles using a cubic lattice with total_mass, semi-axes x, y, z
// No particles closer than min_dist created
void cubic_ellipsoid(struct reb_simulation *n_body_sim, double min_dist,
		double x, double y, double z, double total_mass) {
	// Get particle info
	reb_particle *particles = n_body_sim->particles;

	// Set particle defaults
	struct reb_particle pt;
	pt.vx = 0;
	pt.vy = 0;
	pt.vz = 0;
	pt.ax = 0;
	pt.ay = 0;
	pt.az = 0;
	pt.r = 1.0;
	pt.m = 1.0;

	// Get current number of particles
	int N = n_body_sim->N;
	int i_0 = N;

	// Lattice size
	int nx, ny, nz;
	nx = ny = nz = (int) (1.2 * x / min_dist);

	// Lattice spacing factor
	double fac = min_dist;

	// Get particles in ellipsoid
	for (int k = -nz; k <= nz; k++) {

		// z location
		double cur_z = fac * k;

		for (int j = -ny; j <= ny; j++) {

			// y location
			double cur_y = fac * j;

			for (int i = -nx; i <= nx; i++) {

				// x location
				double cur_x = fac * i;

				// Set particle location
				pt.x = cur_x;
				pt.y = cur_y;
				pt.z = cur_z;

				// Add particle if inside ellipsoid
				Vector pos = { pt.x / x, pt.y / y, pt.z / z };
				if (pos.len() <= 1.0)
					reb_add(n_body_sim, pt);
			}
		}
	}

	// Get current number of particles
	N = n_body_sim->N;

	// Adjust mass of each particle so that they sum to desired total mass
	// Adjust radius of each particle so they have radius of 1/4 smallest distance between any particles
	double particle_mass = total_mass / (N - i_0);
	double min_d = mindist(n_body_sim, i_0, N);
	for (int i = i_0; i < N; i++) {
		particles[i].m = particle_mass;
		particles[i].r = min_d / 4.0;
	}

	std::cout << "hcp_ellipsoid: Created " << N - i_0
			<< " particles, separated by minimum distance " << min_d
			<< std::endl;
}

// Create an approximately uniform random particle distribution with total_mass within shape given by a set of vertices already read in (stored in particles)
// No particles closer than min_dist created
void rand_shape(struct reb_simulation *const n_body_sim, double min_dist,
		double total_mass) {
	// Get particle info
	reb_particle *particles = n_body_sim->particles;

	// Shape model vertices should already have been read in
	// Should be particles [0, N-1]
	int N_shape = n_body_sim->N;

	// Find minimum and maximum radius of arbitrary shape
	double min_r_shape = min_radius(n_body_sim, 0, N_shape);
	double max_r_shape = max_radius(n_body_sim, 0, N_shape);

	// Guess at number of random particles we need to generate
	int n_part = 40 * pow(2.0 * max_r_shape / min_dist, 3.0);

	// Set particle defaults
	struct reb_particle pt;
	pt.ax = 0.0;
	pt.ay = 0.0;
	pt.az = 0.0;
	pt.vx = 0.0;
	pt.vy = 0.0;
	pt.vz = 0.0;
	pt.m = 1.0;
	pt.r = min_dist / 4.0;

	// Get current number of particles
	int N = n_body_sim->N;
	int i_0 = N;

	// Get uniform distribution inside arbitrary shape
	for (int i = 0; i < n_part; i++) {
		pt.x = reb_random_uniform(-max_r_shape, max_r_shape);
		pt.y = reb_random_uniform(-max_r_shape, max_r_shape);
		pt.z = reb_random_uniform(-max_r_shape, max_r_shape);
		bool inside = within_shape(n_body_sim, 0, N_shape, min_r_shape,
				max_r_shape, Vector( { pt.x, pt.y, pt.z }));

		// Check if particle is inside shape
		if (inside) {

			// Is there a particle too nearby?
			bool too_near = false;
			N = n_body_sim->N;
			for (int j = i_0; not too_near && j < N; j++) {
				Vector dx = { pt.x - particles[j].x, pt.y - particles[j].y, pt.z
						- particles[j].z };
				if (dx.len() < min_dist)
					too_near = true;
			}

			// Only add particle if not near any other
			if (not too_near)
				reb_add(n_body_sim, pt);
		}
	}

	// Get current number of particles
	N = n_body_sim->N;

	// Adjust mass of each particle so that they sum to desired total mass
	double particle_mass = total_mass / (N - i_0);
	for (int i = i_0; i < N; i++) {
		particles[i].m = particle_mass;
	}

	// Double-check min_dist
	double min_d = mindist(n_body_sim, i_0, N);
	std::cout << "rand_shape: Created " << N - i_0
			<< "particle separated by a minimum of " << min_d << std::endl;
}

/*************/
/* Modifiers */
/*************/

// Rescale by same factor in each direction
void stretch(struct reb_simulation *const n_body_sim, int i_low, int i_high,
		double scale) {
	for (int i = i_low; i < i_high; i++) {
		n_body_sim->particles[i].x *= scale;
		n_body_sim->particles[i].y *= scale;
		n_body_sim->particles[i].z *= scale;
	}
}

// Rescale by different factor in each direction
void stretch_xyz(struct reb_simulation *const n_body_sim, int i_low, int i_high,
		double x_scale, double y_scale, double z_scale) {
	for (int i = i_low; i < i_high; i++) {
		n_body_sim->particles[i].x *= x_scale;
		n_body_sim->particles[i].y *= y_scale;
		n_body_sim->particles[i].z *= z_scale;
	}
}

/*****************/
/* Mark surfaces */
/*****************/

// Create and fill array that marks particles that are within surf_dist of the surface of arbitrary shape
bool* mark_surf_shrink_int_shape(struct reb_simulation *const n_body_sim,
		int i_low, int i_high, double surf_dist) {
	// Get particle info
	reb_particle *particles = n_body_sim->particles;

	// Initialize vars
	int num_at_surf = 0;
	bool *isSurf = (bool*) malloc(sizeof(bool) * n_body_sim->N);

	// Find particles at surface
	// 0 to i_low
	for (int j = 0; j < i_low; j++) {

		// Get index after removal of shape particles and initialize
		int i = j;
		isSurf[i] = false;

		// Get current position
		Vector x = { particles[j].x, particles[j].y, particles[j].z };

		// Determine closest particle in shape vertices
		int nearest_particle = nearest_to_point(n_body_sim, i_low, i_high, x);

		// Get distance to nearest point
		Vector x_shape = { particles[nearest_particle].x,
				particles[nearest_particle].y, particles[nearest_particle].z };
		Vector dx = x - x_shape;

		// If near surface, set true and add to number of surface particles
		if (dx.len() < surf_dist) {
			isSurf[i] = true;
			num_at_surf++;
			// Otherwise, shrink particle so only surface is visible
		} else {
			particles[j].r = 0.001;
		}
	}

	// i_high to end
	for (int j = i_high; j < n_body_sim->N; j++) {

		// Get index after removal of shape particles and initialize
		int i = j - i_high + i_low;

		// Get current position
		Vector x = { particles[j].x, particles[j].y, particles[j].z };

		// Determine closest particle in shape vertices
		int nearest_particle = nearest_to_point(n_body_sim, 0, i_high, x);

		// Get distance to nearest point
		Vector x_shape = { particles[nearest_particle].x,
				particles[nearest_particle].y, particles[nearest_particle].z };
		Vector dx = x - x_shape;

		// If near surface, set true and add to number of surface particles
		if (dx.len() < surf_dist) {
			isSurf[i] = true;
			num_at_surf++;
			// Otherwise, shrink particle so only surface is visible
		} else {
			isSurf[i] = false;
			particles[j].r = 0.001;
		}
	}

	// Return boolean array of which particles are on surface
	std::cout << "Number of vertices on surface of shape: " << num_at_surf
			<< std::endl;
	return isSurf;
}

// Mark particles near the surface of a cone
bool* mark_surf_shrink_int_cone(struct reb_simulation *const n_body_sim,
		double surf_dist, double radius, double height) {
	// Get particle info
	reb_particle *particles = n_body_sim->particles;

	// Calculate slope
	double slope = height / radius;

	// Initialize vars
	int num_at_surf = 0;
	bool *isSurf = (bool*) malloc(sizeof(bool) * n_body_sim->N);

	// Check if particles are close to surface
	double ifac = sqrt(1.0 + slope * slope);
	for (int j = 0; j < n_body_sim->N; j++) {

		// Get particle location and cone location
		double x = particles[j].x;
		double y = particles[j].y;
		double z = particles[j].z;
		double cur_rad = sqrt(x * x + y * y);
		double max_height = cur_rad * slope;

		// Distances to surface??????
		double dplus = abs(z + max_height - height) / ifac;
		double dminus = abs(z - max_height + height) / ifac;

		// If near surface, mark and increment
		if ((dplus < surf_dist) || (dminus < surf_dist)) {
			isSurf[j] = true;
			num_at_surf++;
			// Otherwise, mark and shrink so only surface particles show
		} else {
			isSurf[j] = false;
			particles[j].r = 0.001;
		}
	}

	// Return boolean array of which particles are on surface
	std::cout << "Number of vertices on surface of cone: " << num_at_surf
			<< std::endl;
	return isSurf;
}

// Mark particles near the surface of an ellipsoid defined by semi-axes x, y, z
bool* mark_surf_shrink_int_ellipsoid(struct reb_simulation *const n_body_sim,
		double surf_dist, double x, double y, double z) {
	// Get particle info
	reb_particle *particles = n_body_sim->particles;

	// Initialize vars
	int num_at_surf = 0;
	bool *isSurf = (bool*) malloc(sizeof(bool) * n_body_sim->N);

	// Check if particle is close to surface
	for (int i = 0; i < n_body_sim->N; i++) {

		// Get particle position
		double cur_x = particles[i].x;
		double cur_y = particles[i].y;
		double cur_z = particles[i].z;
		Vector pos = { cur_x / x, cur_y / y, cur_z / z };

		// Distance from surface
		double dist = 1.0 - pos.len();

		// If at surface, mark and increment number of surface particles
		if ((dist < surf_dist) && (dist > 0)) {
			isSurf[i] = true;
			num_at_surf++;
			// Otherwise, shrink so that only surface is visible
		} else {
			isSurf[i] = false;
			n_body_sim->particles[i].r = 0.001;
		}
	}

	// Return boolean array of which particles are on surface
	std::cout << "Number of vertices on surface of ellipsoid: " << num_at_surf
			<< std::endl;
	return isSurf;
}

/***********/
/* Helpers */
/***********/

// Remove particles [i_low, i_high)
void rm_particles(struct reb_simulation *const n_body_sim, int i_low,
		int i_high) {
	// Throw error if low index is above high index
	if (i_low > i_high)
		throw;

	// Remove i_low repeatedly
	// Particles are shifted downward after each removal
	const int to_remove = i_low;
	const int retain_order = 1;
	for (int i = i_low; i < i_high; i++) {
		reb_remove(n_body_sim, to_remove, retain_order);
	}
}

// Is position x,y,z within the shape given by vertices [i_low, i_high-1]?
// r_min, r_max are the maximum and minimum radii of shape vertices, to short-circuit test
bool within_shape(struct reb_simulation *n_body_sim, int i_low, int i_high,
		double r_min, double r_max, Vector pos) {
	// Get particle info
	reb_particle *particles = n_body_sim->particles;

	// Get radius of point
	double rad = pos.len();

	// If radius is inside smallest or outside largest, we can return
	if (rad < r_min) {
		return true;	// is within shape
	} else if (rad > r_max) {
		return false;	// is outside shape
	}

	// Otherwise, find index of nearest particle in shape model
	int particle = nearest_to_point(n_body_sim, i_low, i_high, pos);

	// Compute radius of this particle in shape model
	Vector part_pos = { particles[particle].x, particles[particle].y,
			particles[particle].z };

	// Determine if inside or outside shape
	if (rad > part_pos.len()) {
		return false; // outside shape locally
	} else {
		return true;  // inside shape locally
	}
}

// Return index of closest particle to position x,y,z within [i_low, i_high)
// Caution: Origin might not be center of body
int nearest_to_point(struct reb_simulation *n_body_sim, int i_low, int i_high,
		Vector pos) {
	// Get particle info
	reb_particle *particles = n_body_sim->particles;

	// Initialize squared distance
	double dist = DBL_MAX;

	// Find nearest particle
	int i_0;
	for (int i = i_low; i < i_high; i++) {
		Vector part_pos = { particles[i].x, particles[i].y, particles[i].z };
		Vector dx = part_pos - pos;
		double r = dx.len();
		if (r < dist) {
			i_0 = i;
			dist = r;
		}
	}

	// Return index of nearest particle
	return i_0;
}

// Return minimum radius of particles in index range [i_low, i_high)
double min_radius(struct reb_simulation *const n_body_sim, int i_low,
		int i_high) {
	// Get particle info
	reb_particle *particles = n_body_sim->particles;

	// Find smallest radius of particles
	double min_r = DBL_MAX;
	for (int i = i_low; i < i_high; i++) {
		Vector x = { particles[i].x, particles[i].y, particles[i].z };
		double r = x.len();
		if (r < min_r)
			min_r = r;
	}

	// Return smallest radius
	return min_r;
}

// Return maximum radius of particles in index range [i_low, i_high)
double max_radius(struct reb_simulation *const n_body_sim, int i_low,
		int i_high) {
	// Get particle info
	reb_particle *particles = n_body_sim->particles;

	// Find largest radius of particles
	double max_r = 0.0;
	for (int i = i_low; i < i_high; i++) {
		Vector x = { particles[i].x, particles[i].y, particles[i].z };
		double r = x.len();
		if (r > max_r)
			max_r = r;
	}

	// Return largest radius
	return max_r;
}

// Compute the shortest distance between any pair of particles in range [i_low, i_high)
double mindist(struct reb_simulation *const n_body_sim, int i_low, int i_high) {
	// Get particle info
	struct reb_particle *particles = n_body_sim->particles;

	// Find shortest distance
	double min_dist = DBL_MAX;
	for (int i = i_low; i < i_high - 1; i++) {
		Vector x1 = { particles[i].x, particles[i].y, particles[i].z };
		for (int j = i + 1; j < i_high; j++) {
			Vector x2 = { particles[j].x, particles[j].y, particles[j].z };
			Vector dx = x1 - x2;
			double dist = dx.len();
			if (dist < min_dist)
				min_dist = dist;
		}
	}

	// Return shortest distance
	return min_dist;
}
