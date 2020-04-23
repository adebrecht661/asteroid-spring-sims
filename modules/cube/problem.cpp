#ifdef __cplusplus
# 	ifdef __GNUC__
#		define restrict __restrict__
#	else
#		define restrict
#	endif
#endif

/**
 * resolved mass spring model
 * using the leap frog integrator.
 */
#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
extern "C" {
#include "rebound.h"
}
#include "matrix_math.h"
#include "input_spring.h"
#include "output_spring.h"
#include "shapes.h"
#include "springs.h"
#include "physics.h"
#include "stress.h"

using std::string;
using std::vector;

int num_springs;  // number of springs
vector<spring> springs;  // springs structure
double surfdist;       // for identifying surface particles
double Pressure; // vertical pressure

double def_gamma;
double gamma_fac; // for adjustment of gamma of all springs
double t_damp;    // end faster damping, relaxation
double print_interval;   // for printouts

string fileroot;   // output files

int num_perts = 0;

// Forward declarations
void reb_springs(reb_simulation *const n_body_sim); // to pass springs to display
void table_top(reb_simulation *const n_body_sim);
void top_push(reb_simulation *const n_body_sim);
void heartbeat(reb_simulation *const n_body_sim);
void additional_forces(reb_simulation *n_body_sim);

int main(int argc, char *argv[]) {
	struct reb_simulation *const n_body_sim = reb_create_simulation();
	// Setup constants
	n_body_sim->integrator = reb_simulation::REB_INTEGRATOR_LEAPFROG;
	n_body_sim->gravity = reb_simulation::REB_GRAVITY_NONE;
	n_body_sim->boundary = reb_simulation::REB_BOUNDARY_NONE;
	n_body_sim->G = 0;
	n_body_sim->additional_forces = additional_forces; // setup callback function for additional forces
	double mcube = 1.0;          // total mass of all particles
	double rcube = 1.0;          // length of cube
	double tmax = 0.0;           // if 0 integrate forever

// things to set! ////////////////////// could be read in with parameter file
	double dt;             // timestep
	double b_distance;     // mininum interparticle spacing
	double mush_fac1;       // mush_fac*b_distance is maximum spring length
	double mush_fac2;       //
	double mush_fac3;       //
	def_gamma = 0.01;  // damping parm
	double ks1 = 0.0;      // spring constant
	double ks2 = 0.0;      // spring constants
	double ks3 = 0.0;      // spring constants

	if (argc == 1) {
		fileroot = "t1";   // to make output files
		dt = 1e-3;    // Timestep
		b_distance = 0.15; // for creating random sphere, min separation between particles
		mush_fac1 = 2.3; // ratio of smallest spring distance to minimum interparticle dist
		mush_fac2 = 2.3; // ratio of smallest spring distance to minimum interparticle dist
		mush_fac3 = 2.3; // ratio of smallest spring distance to minimum interparticle dist
		// spring damping
		def_gamma = 1.0;    // initial factor for initial damping value for springs
		t_damp = 1.0;    // gamma to final values for all springs at this time
		ks1 = 8e-2;   // spring constant
		ks2 = 8e-2;   // spring constant
		ks3 = 8e-2;   // spring constant
		print_interval = 100000.0;  // printouts
		surfdist = 0.1;         // for identifying surface particles
		Pressure = 0.0;         // apply presure
	} else {
		FILE *fpi; // read in a parameter file
		fpi = fopen(argv[1], "r");
		char line[300];
		char froot[100];
		fgets(line, 300, fpi);
		sscanf(line, "%s", froot);
		fileroot = froot;
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &dt);
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &tmax);
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &print_interval);
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &b_distance);
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &mush_fac1);
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &mush_fac2);
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &mush_fac3);
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &ks1);
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &ks2);
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &ks3);
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &def_gamma);
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &gamma_fac);
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &t_damp);
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &surfdist);
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &Pressure);

		printf("parm file read in\n");

	}

/// end of things to set /////////////////////////

	n_body_sim->dt = dt;            // set integration timestep
	const double boxsize = 1.1 * rcube;    // display window
	reb_configure_box(n_body_sim, boxsize, 1, 1, 1);
	n_body_sim->softening = b_distance / 100.0;	// Gravitational softening length
// viewer +x to right, +y to up, z back and forth along line of sight

	struct spring spring_mush; // spring parameters for mush
	// properties of springs
	spring_mush.gamma = def_gamma; // damping coefficient
	spring_mush.k = ks1; // spring constant
	spring_mush.k_heat = 1.0;    // heat transport coefficient
	double mush_distance1 = b_distance * mush_fac1;
	double mush_distance2 = b_distance * mush_fac2;
	double mush_distance3 = b_distance * mush_fac3;
	// distance for connecting and reconnecting springs

	char fname[200];
	sprintf(fname, "%s_run.txt", fileroot.c_str()); // for simulation info

	num_springs = 0; // start with no springs

	// create particle distribution

	rand_rectangle(n_body_sim, b_distance, rcube, rcube, rcube, mcube);

	// make springs, all pairs connected within interparticle distance mush_distance
	connect_springs_dist(n_body_sim, mush_distance1, 0, n_body_sim->N, spring_mush);
	spring_mush.k = ks2; // add longer springs!
	connect_springs_dist(n_body_sim, mush_distance2, 0, n_body_sim->N, spring_mush);
	spring_mush.k = ks3; // add longer springs!
	connect_springs_dist(n_body_sim, mush_distance3, 0, n_body_sim->N, spring_mush);

	reb_springs(n_body_sim); // pass spring index list to display
	set_gamma(def_gamma*gamma_fac);  // start with enhanced damping

	n_body_sim->heartbeat = heartbeat;

	if (tmax == 0.0) // start the integration!!!!
		reb_integrate(n_body_sim, INFINITY);
	else
		reb_integrate(n_body_sim, tmax);
}

void heartbeat(struct reb_simulation *const r) {
	static int index = 0;
	if (reb_output_check(r, 10.0 * r->dt)) {
		reb_output_timing(r, 0); // print time of simulation run on screen
	}
	if (abs(r->t - t_damp) < 0.9 * r->dt) {
		set_gamma(def_gamma);
	}
	// damp initial bounce only , end damping
	// reset gamma only at t near t_damp

	// stuff to do every timestep
	// nothing!

	if (reb_output_check(r, print_interval)) {
		write_particles(r, fileroot, index); //   output particle positions
		index++;
	}

}

void additional_forces(reb_simulation *n_body_sim) {
	zero_accel(n_body_sim);
	spring_forces(n_body_sim); // spring forces
	top_push(n_body_sim); // top force
	table_top(n_body_sim); // bottom force
}

// make a spring index list to pass to viewer
void reb_springs(struct reb_simulation *const r) {
	r->NS = num_springs;
	r->springs_i = (int*)malloc(num_springs * sizeof(int));
	r->springs_j = (int*)malloc(num_springs * sizeof(int));
	for (int i = 0; i < num_springs; i++) {
		r->springs_i[i] = springs[i].particle_1;
		r->springs_j[i] = springs[i].particle_2;
	}
}

// I need to mark top surface in the first call as it can move!
void top_push(struct reb_simulation *const n_body_sim) {
	static int *surflist; // list of top surface  particles!
	static int Nsurf = 0;  // number of particles on top surface
	static int first = 0;
	if (first == 0) { // first time the routine is called
		first = 1;
		surflist = (int*)malloc(sizeof(int) * n_body_sim->N);
		for (int i = 0; i < n_body_sim->N; i++) {
			if (n_body_sim->particles[i].y > 0.5 - surfdist) {
				surflist[Nsurf] = i;
				Nsurf++;
			}
		}
		printf("Nsurf top =%d\n", Nsurf);
	}
	for (int j = 0; j < Nsurf; j++) {
		int i = surflist[j];
		n_body_sim->particles[i].ay -= Pressure / Nsurf
				/ n_body_sim->particles[i].m;
	}
}

void table_top(struct reb_simulation *const r) {
	static int *surflist; // list of bottom surface  particles!
	static int Nsurf = 0;  // number of particles on bottom surface
	static int first = 0;
	if (first == 0) { // first time the routine is called
		first = 1;
		surflist = (int*)malloc(sizeof(int) * r->N);
		for (int i = 0; i < r->N; i++) {
			if (r->particles[i].y < -0.5 + surfdist) {
				surflist[Nsurf] = i;
				Nsurf++;
			}
		}
		printf("Nsurf bottom=%d\n", Nsurf);
	}
	for (int j = 0; j < Nsurf; j++) {
		int i = surflist[j];
		r->particles[i].ay = 0.0;
	}
}
