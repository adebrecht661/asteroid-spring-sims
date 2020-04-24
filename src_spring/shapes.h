/*
 * shapes.h
 *
 *	Create and modify distributions of particles
 *
 *  Created on: Apr 3, 2020
 *      Author: alex
 */

#ifndef SRC_SPRING_SHAPES_H_
#define SRC_SPRING_SHAPES_H_

/**************/
/* Generators */
/**************/

// Create particles uniformly spaced along a line of length 1 in the y direction, with an initial offset of y_off
void uniform_line(reb_simulation* const n_body_sim, int num_parts, double y_off);
// Create particles approximately evenly distributed inside rectangular prism defined by sides x, y, z, no closer than min_dist, with total_mass
void rand_rectangle(reb_simulation *const n_body_sim, double min_dist,
		double x, double y, double z, double total_mass);
// Create particles approximately evenly distributed inside rectangle (2D) defined by sides x, y, no closer than min_dist, with total_mass
void rand_rectangle_2d(reb_simulation *const n_body_sim, double min_dist,
		double x, double y, double total_mass);
// Create particles approximately evenly distributed inside cone of radius, height, no closer than min_dist, with total_mass
void rand_cone(reb_simulation *const n_body_sim, double min_dist,
		double radius, double height, double total_mass);
// Create particles approximately evenly distributed inside ellipsoid of semi-axes x, y, z, no closer than min_dist, with total_mass
void rand_ellipsoid(reb_simulation *const n_body_sim, double min_dist,
		double x, double y, double z, double total_mass);
// Create particles in hexagonal close-packed lattice inside ellipsoid of semi-axes x, y, z, no closer than min_dist, with total_mass
void hcp_ellipsoid(reb_simulation *n_body_sim, double min_dist,
		double x, double y, double z, double total_mass);
// Create particles in cubic lattice inside ellipsoid of semi-axes x, y, z, no closer than min_dist, with total_mass
void cubic_ellipsoid(reb_simulation *n_body_sim, double min_dist,
		double x, double y, double z, double total_mass);
// Create particles approximately evenly distributed inside shape from vertices already read in
// CAUTION: Make sure to remove particles that store shape vertices
void rand_shape(reb_simulation *const n_body_sim, double min_dist,
		double total_mass);

/*************/
/* Modifiers */
/*************/

// Stretch particles in [i_low, i_high) by same amount in each direction
void stretch(reb_simulation *const n_body_sim, int i_low, int i_high,
		double scale);
// Stretch particles in [i_low, i_high) by different amount in each direction
void stretch_xyz(reb_simulation *const n_body_sim, int i_low, int i_high,
		double x_scale, double y_scale, double z_scale);

/*****************/
/* Mark surfaces */
/*****************/

// Mark surface particles and shrink interior particles of arbitrary shape defined by particles in set [i_low, i_high)
void mark_surf_shrink_int_shape(reb_simulation *const n_body_sim,
		int i_low, int i_high, double surf_dist);
// Mark surface particles and shrink interior particles of cone of radius and height
void mark_surf_shrink_int_cone(reb_simulation *const n_body_sim,
		double surf_dist, double radius, double height);
// Mark surface particles and shrink interior particles of ellipsoid with semi-axes x, y, z
void mark_surf_shrink_int_ellipsoid(reb_simulation *const n_body_sim,
		double surf_dist, double x, double y, double z);

/***********/
/* Helpers */
/***********/

// Remove particles in set [i_low, i_high) from simulation
void rm_particles(reb_simulation *const n_body_sim, int i_low,
		int i_high);
// Determine if pos is inside shape defined by [i_low, i_high)
bool within_shape(reb_simulation *n_body_sim, int i_low, int i_high,
		double r_min, double r_max, Vector pos);
// Find nearest particle to pos
int nearest_to_point(reb_simulation *n_body_sim, int i_low, int i_high,
		Vector pos);
// Return minimum distance from <0,0,0> in particle set [i_low, i_high)
double min_radius(reb_simulation *n_body_sim, int i_low, int i_high);
// Return maximum distance from <0,0,0> in particle set [i_low, i_high)
double max_radius(reb_simulation *n_body_sim, int i_low, int i_high);
// Return minimum distance between any two particles in set [i_low, i_high)
double mindist(reb_simulation *const n_body_sim, int i_low, int i_high);

#endif /* SRC_SPRING_SHAPES_H_ */
