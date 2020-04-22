/*
 * output.h
 *
 *  Created on: Apr 15, 2020
 *      Author: alex
 */

#ifndef SRC_SPRING_OUTPUT_H_
#define SRC_SPRING_OUTPUT_H_

/********************/
/* Output functions */
/********************/

// Write information about surface particles to specified filename
// t x y z vx vy vz ax ay az r m
void write_surf_part(struct reb_simulation *const n_body_sim, int i_low, int i_high, bool is_surf,
		string filename);
// Write information about resolved body to specified filename
// t x y z vx vy vz omx omy omz llx lly llz Ixx Iyy Izz Ixy Iyz Ixz KErot PEspr PEgrav Etot dEdtnow dEdtave
void write_resolved_with_E(struct reb_simulation *const n_body_sim, int i_low, int i_high,
		string filename, double dEdt_ave);
// Write information about resolved body to specified filename
// t x y z vx vy vz omx omy omz llx lly llz Ixx Iyy Izz Ixy Iyz Ixz
void write_resolved_no_E(struct reb_simulation *const n_body_sim, int i_low, int i_high,
		string filename, double dEdt_ave);
// Write information about particle with largest x value and particle with largest z value at first print
// t x y z vx vy vz x y z vx vy vz
void write_resolved_2nodes(struct reb_simulation *const n_body_sim, int i_low,
		int i_high, string filename);
// Write information about point mass perturber
// i_pert is index in particles array, pert_num is number out of total point masses
// t x y z vx vy vz m
void write_pt_mass(struct reb_simulation *const n_body_sim, int i_pert, int pert_num,
		string filename);
// Write information about resolved body orbiting one or more perturbing particles
// t a n e i omx omy omz eigs[0] eigs[1] eigs[2] E lx ly lz ang lox loy loz Ixx Iyy Izz Ixy Iyz Ixz
void write_resolved_orb(struct reb_simulation *const n_body_sim, string filename);
// Write information about resolved body and point mass with binary info, if present
void write_resolved_bin(struct reb_simulation *const n_body_sim, string filename);
// Write spring info to file fileroot_%06d_springs.txt
void write_springs(struct reb_simulation *const n_body_sim, string fileroot, int index);
// Write particle info to file fileroot_%06d_particles.txt
void write_particles(struct reb_simulation *const n_body_sim, string fileroot, int index);
// Write all node info to file
void write_nodes(struct reb_simulation *const n_body_sim, string filename);
// Write all stress info to file
void write_stresses(struct reb_simulation *const n_body_sim, string filename);

/********************/
/* Filename helpers */
/********************/

// Get filename for node file
void node_filename(struct reb_simulation *const n_body_sim, string root, double num_to_print);
// Get filename for stress file
void stress_filename(struct reb_simulation *const n_body_sim, string root, double num_to_print);

#endif /* SRC_SPRING_OUTPUT_H_ */
