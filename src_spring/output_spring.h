/*
 * output.h
 *
 *  Created on: Apr 15, 2020
 *      Author: alex
 */

#ifndef SRC_SPRING_OUTPUT_SPRING_H_
#define SRC_SPRING_OUTPUT_SPRING_H_

#include <string>
using std::string;

/********************/
/* Output functions */
/********************/

// Write information about surface particles to specified filename
// t x y z vx vy vz ax ay az r m
void write_surf_part(reb_simulation *const n_body_sim, int i_low, int i_high,
		string filename);
// Write information about resolved body to specified filename
// t x y z vx vy vz omx omy omz llx lly llz Ixx Iyy Izz Ixy Iyz Ixz KErot PEspr PEgrav Etot dEdtnow dEdtave
void write_resolved_with_E(reb_simulation *const n_body_sim, int i_low,
		int i_high, string filename, double dEdt_ave);
// Write information about resolved body to specified filename
// t x y z vx vy vz omx omy omz llx lly llz Ixx Iyy Izz Ixy Iyz Ixz
void write_resolved_no_E(reb_simulation *const n_body_sim, int i_low,
		int i_high, string filename);
// Write information about particle with largest x value and particle with largest z value at first print
// t x y z vx vy vz x y z vx vy vz
void write_resolved_2nodes(reb_simulation *const n_body_sim, int i_low,
		int i_high, string filename);
// Write information about point mass perturber
// i_pert is index in particles array, pert_num is number out of total point masses
// t x y z vx vy vz m
void write_pt_mass(reb_simulation *const n_body_sim, int i_pert, int pert_num,
		string filename);
// Write information about resolved body orbiting one or more perturbing particles
// t a n e i omx omy omz eigs[0] eigs[1] eigs[2] E lx ly lz ang lox loy loz Ixx Iyy Izz Ixy Iyz Ixz
void write_resolved_orb(reb_simulation *const n_body_sim, string filename);
// Write information about resolved body and point mass with binary info, if present
void write_resolved_bin(reb_simulation *const n_body_sim, string filename);
// Write spring info to file fileroot_%06d_springs.txt
void write_springs(reb_simulation *const n_body_sim, string fileroot,
		int index);
// Write particle info to file fileroot_%06d_particles.txt
void write_particles(reb_simulation *const n_body_sim, string fileroot,
		int index);
// Write all node info to file
void write_nodes(reb_simulation *const n_body_sim, string filename);
// Write all stress info to file
void write_stresses(reb_simulation *const n_body_sim, string filename);
// Write all heat info to file
void write_heat(reb_simulation *const n_body_sim, string filename,
		int num_timesteps, double power_fac);
// Print doubles to file and standard out
void print_run_double(double quantity, string label, std::ofstream *outfile);

/****************/
/* Heat helpers */
/****************/

// Normalize total power by number of timesteps for output
void normalize_tot_power(double ndt);

/********************/
/* Filename helpers */
/********************/

// Get filename for heat file
string heat_filename(reb_simulation *const n_body_sim, string root,
		double print_interval);
// Get filename for node file
string node_filename(reb_simulation *const n_body_sim, string root,
		double print_interval);
// Get filename for stress file
string stress_filename(reb_simulation *const n_body_sim, string root,
		double print_interval);

#endif /* SRC_SPRING_OUTPUT_SPRING_H_ */
