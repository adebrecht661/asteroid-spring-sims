/**
 * @file        spring.h
 * @brief       springs 
 * @author      Alice Quillen
 */

#ifndef _SPRING_H
#define _SPRING_H

#include <stdbool.h>
#include "particle.h"
#include "matrix_math.h"

// Spring structure
typedef struct spring {
	double ks;   // spring constant
	double rs0;  // distance of no force
	double gamma; // damping coefficient
	double k_heat;  //  heat diffusion coefficient!
	int i;    // vertex 1 referring to a particle
	int j;    // vertex 2 referring to a particle
} spring;

// Coordination structure for each node
typedef struct coord {
	int ncoord;  // number of springs connected to this node
	int coordlist[50]; // a list of spring indices ** note number here is hard wired, giving a max!
	double lx;
	double ly;
	double lz;
} coord;

typedef struct node {
	int surf;  // is 1 if is a surface node
	double temp;  // temperature
	double cv;    // specific heat, probably integrated over mass of node?
} node;

typedef struct symtensor {
	double xx;
	double yy;
	double zz;
	double xy;
	double yz;
	double xz;
} symtensor;

extern struct spring *springs;
extern int NS; // numbers of springs
extern int NPERT; // number of point masses
extern double b_distance; // mush formation
extern double mush_distance; // mush spring connection 
extern double t_reform; // formation springs timescale
extern double gamma_all; // for gamma  of all springs

/*********************/
/* Spring operations */
/*********************/

// Delete spring i
void del_spring(struct reb_simulation *const n_body_sim, int i);
// Add a spring with properties of spr between particle 1 and particle 2. Returns index of spring connecting these particles.
int add_spring(struct reb_simulation *const n_body_sim, int particle_1,
		int particle_2, spring spr);
// Helper function to add a spring at end of current array, expand if required (no sanity checking)
void add_spring_helper(spring spr);
// Get length of passed spring
double spring_length(struct reb_simulation *const n_body_sim, spring spr);
// Compute spring midpoint location from arbitrary center in Cartesian coordinates
Vector spr_mid(struct reb_simulation *const n_body_sim, spring spr,
		Vector center);
// Returns strain on given spring
double strain(struct reb_simulation *const n_body_sim, spring spr);
// Connects springs between all particles closer than dist in index range i_low -> i_high-1
void connect_springs_dist(struct reb_simulation *const n_body_sim,
		double max_dist, int i_low, int i_high, spring spr);
// Set damping coefficient of all springs
void set_gamma(double new_gamma);
void connect_springs_dist_nodemax();

/****************************************/
/* !!!!!! Unsorted operations !!!!!!!!! */
/****************************************/

void spring_forces(struct reb_simulation* n_body_sim);
double spring_potential_energy(struct reb_simulation *const r);
double grav_potential_energy(struct reb_simulation *const r, int il, int ih);

// list of springs related subroutines
struct node* mknodevec();
void surface_nodes();
void nfilename();
void print_node();
void transport_heat();
void heat_nodes_tidal();
void heat_nodes_radiogenic();
double Kthermal_mush();
void adjust_spring_temp_lin();
void adjust_spring_temp_ts();

void spring_force_k_one();
void zero_symtensor();

void print_stress();
void update_stresstensor();
int markfailure();
void sfilename();
void killsprings();

int nearest_to_shape(struct reb_simulation *r, int il, int ih, double x,
		double y, double z);
void print_surf();
int surface_shape(); //ZYH
void surfaceparticle_display(); //ZYH
void potoang(); //ZYH
int* marksurface();
int* marksurface_football();
int* marksurface_cone();
void rescale_xyz();
double min_radius();
double max_radius();
void rand_bennu();
void rmvertices();
void rmshape_vertices();
void read_vertex_file();
void adjust_ks_abc();
void adjust_ks_abc_fac();
void adjust_mass_abc();
void subtractcov(struct reb_simulation *const r, int il, int ih);
void subtractcom(struct reb_simulation *const r, int il, int ih);
void print_extended();
void print_extended_simp();
void print_extended_2nodes();
void print_pm();
void hfilename();
void quadJ2pole();
double add_pt_mass_kep();
void move_resolved(struct reb_simulation *r, double dx, double dy, double dz,
		double dvx, double dvy, double dvz, int il, int ih);
double sum_mass(struct reb_simulation *const r, int il, int ih);
double dEdt(struct reb_simulation *const r, struct spring spr);
double dEdt_total(struct reb_simulation *const r);
Vector compute_com(struct reb_simulation *const n_body_sim, int il, int ih);
void compute_cov(struct reb_simulation *const r, int il, int ih, double *vxc, double *vyc, double *vzc);
void compute_Lorb(struct reb_simulation *const r, int il, int ih, int npert, double *llx, double *lly, double *llz);
double mindist();
void centerbody(struct reb_simulation *const r, int il, int ih);
void rand_football();
void rand_rectangle();
void rand_rectangle_2d();
void rand_cone();
void rand_football_from_sphere();
double Young_mesh(struct reb_simulation *const r, int il, int ih, double rmin,
		double rmax);
double Young_mesh_big(struct reb_simulation *const r, int il, int ih);
void spin(struct reb_simulation *const r, int il, int ih, Vector omega);
void make_binary_spring();
Matrix mom_inertia(struct reb_simulation *const n_body_sim, int i_low, int i_high);
Vector measure_L(struct reb_simulation *const r, int il, int ih);
void measure_L_origin();
void compute_semi(struct reb_simulation *const r, int il, int ih, int im1,
		double *aa, double *meanmo, double *ee, double *ii, double *LL);
void compute_semi_bin(struct reb_simulation *const r, int il, int ih, int npert,
		double *aa, double *meanmo, double *ee, double *ii, double *LL);
void total_mom();
double mean_L(struct reb_simulation *const r);
void spring_init(struct reb_simulation *r);
void output_png();
void output_png_single();
Vector body_spin(struct reb_simulation *const r, int il, int ih, double eigs[3]);
void print_tab();
void print_bin();
void print_heat();
void write_particles(struct reb_simulation *const r, char *fileroot, int index);
void write_springs(struct reb_simulation *const r, char *fileroot, int index);
void toistring(char *istring, int i);
void zero_accel(struct reb_simulation* n_body_sim);
void rotate_to_principal();
void vec_mul(double Axx, double Axy, double Axz, double Ayx, double Ayy,
		double Ayz, double Azx, double Azy, double Azz, double bx, double by,
		double bz, double *cx, double *cy, double *cz);
void adjust_ks();
void adjust_mass_side();
void rotate_body();
void rotate_origin();
void rotate_vector();
double fill_hcp();
double fill_cubic();
double add_pluto_charon();
double add_pluto_charon_kep();
double add_one_mass_kep();
void add_one_mass_cartesian();
void dodrift_bin();
void dodrift_res();
void addto_heatvec();
void norm_heatvec(int ndt);
double compute_rot_kin(struct reb_simulation *const r, int il, int ih);
void adjust_nodes_cp();
void stretch();

#endif // _SPRING_H

