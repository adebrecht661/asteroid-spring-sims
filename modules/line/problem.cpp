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
#include "heat.h"

using std::string;
using std::vector;

// needed for spring library and hacked rebound viewer
int num_springs = 0;
vector<spring> springs;  // springs structure
extern vector<node> nodes;

void heartbeat(struct reb_simulation *const r);  // runtime subroutine
void gen_particles(struct reb_simulation *const r, int NN); // generate particles
void sphere_forces(struct reb_simulation *const r); // interaction forces
void pulse_base(struct reb_simulation *const r); // boundary condition

// icky globals   because rebound doesn't make it easy to pass anything
string fileroot;   // for output files
string pfilename;  // for a particular particle output
int p_particle; // printing out info on this particle
double print_interval;   // outputs on this time interval
double tp_print;   // outputs for p_particler on this time interval
double L_overlap; // for springs, overlap distance
double k_spring; // spring constant
double alpha;   // damping constant being used now
double alpha_run;  // damping during run
double alpha_damp; // damping at beginning of run
double t_damp;   // time for early higher damping
double G_grav;   //  downward gravity
double py0;     //  if first particle is below this then push it up
double A_pulse;  // pulse amplitude
double tau_pulse;  // pulse time width
double tau_wait;  // time to wait before doing another pulse

int num_perts = 0;

// passed to rebound so you can have forces other than gravity
void additional_forces(struct reb_simulation *r) {
	zero_accel(r); // because gravity routines in rebound set acceleration and if you don't
				   // call a gravity routine you need to zero the acceleration array
	sphere_forces(r); // sof sphere forces

}

// primarily set things up so we can call integrator and display
int main(int argc, char *argv[]) {
	struct reb_simulation *const r = reb_create_simulation();
	// Setup constants
	r->integrator = reb_simulation::REB_INTEGRATOR_LEAPFROG;
	r->gravity = reb_simulation::REB_GRAVITY_NONE;
	r->boundary = reb_simulation::REB_BOUNDARY_NONE;
	r->G = 0;
	r->additional_forces = additional_forces; // setup callback function for additional forces
	double rcube = 1.0;          // a length scale
	double tmax = 0.0;           // if 0 integrate forever

	double dt;             // timestep
	int NN;  // number of particles in the line

// things to set! ////////////////////// could be read in with parameter file

	if (argc == 1) { // if you run it without passing a filename on command line
		fileroot = "t1";   // to make output files
		dt = 1e-6;    // Timestep
		NN = 50; // number of particles
		print_interval = 1.0e6;
		tp_print = 1.0e6;
		k_spring = 1.0e5;
		G_grav = 1.000;
		alpha_run = 1.0e2;
		alpha_damp = 1.0e4;
		t_damp = 3.0;
		py0 = -1.7;
		A_pulse = 0.3;
		tau_pulse = 1.0;
		tau_wait = 6.0;
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
		sscanf(line, "%lf", &tp_print);
		fgets(line, 300, fpi);
		sscanf(line, "%d", &NN);
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &k_spring);
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &G_grav);
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &alpha_run);
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &alpha_damp);
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &t_damp);
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &py0);
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &A_pulse);
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &tau_pulse);
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &tau_wait);

		printf("parm file read in\n");

	}
	alpha = alpha_damp;

/// end of things to set /////////////////////////
	gen_particles(r, NN); // generate some particles

	r->dt = dt;            // set integration timestep
	const double boxsize = 1.1 * rcube;    // display window
	reb_configure_box(r, boxsize, 1, 2, 1);
// viewer +x to right, +y to up, z back and forth along line of sight
// I have put all motion in the y direction to be consistent with these
// orientations

//   FILE *fpr;
//   char fname[200];
//   sprintf(fname,"%s_run.txt",froot); // for simulation info
//   fpr = fopen(fname,"w");

	pfilename = fileroot + "_p" + std::to_string((int) NN / 2) + ".txt";

	r->heartbeat = heartbeat; // tell rebound about the hearbeat routine

	if (tmax == 0.0) // start the integration!!!! (and display)
		reb_integrate(r, INFINITY);
	else
		reb_integrate(r, tmax);
}

// this is called every timestep while the integration is running
void heartbeat(struct reb_simulation *const r) {
	struct reb_particle *particles = r->particles;
	static int index = 0;
	static FILE *fpp;
	if (index == 0) {
		fpp = fopen(pfilename.c_str(), "w"); // open it once!
	}

	if (reb_output_check(r, 10.0 * r->dt)) {
		reb_output_timing(r, 0); // print time of simulation run on screen
	}

	if ((r->t > 0.0) && (reb_output_check(r, t_damp))) {
		alpha = alpha_run;
	}

	if (reb_output_check(r, print_interval)) {
		write_particles(r, fileroot, index); //   output particle positions
		index++;
	}

	// output a lot of information on one particle in the middle
	if (r->t - t_damp > 0) {
		if (reb_output_check(r, tp_print)) {
			fprintf(fpp, "%.6e ", r->t - t_damp);
			fprintf(fpp, "%.6f ", particles[p_particle].y);
			fprintf(fpp, "%.6f ", particles[p_particle].vy);
			fprintf(fpp, "%.6f ", particles[p_particle].ay);
			fprintf(fpp, "\n");
		}
	}

	pulse_base(r); // boundary condition

}

// generate particles in a line
void gen_particles(struct reb_simulation *const r, int NN) {
	struct reb_particle pt;  // particle structure
	double dz = 1.0 / NN; // spacing
	pt.m = dz; // mass, sums to 1
	pt.x = 0.0;
	pt.y = 0.0;
	pt.z = 0.0;  // position
	pt.vx = 0.0;
	pt.vy = 0.0;
	pt.vz = 0.0;  // velocity
	pt.ax = 0.0;
	pt.ay = 0.0;
	pt.az = 0.0;  // acceleration
	pt.r = dz / 2; // display radius
	for (int i = 0; i < NN; i++) {
		double z0 = dz * i;
		pt.y = z0 + py0;   // up is y!
		reb_add(r, pt);
	}
	L_overlap = 1.0 / NN;   // set global variable for interaction forces

	// If node vector has already been initialized, add nodes for new particles
	if (!nodes.empty()) {
		node zero_node;
		zero_node.cv = -1;
		zero_node.is_surf = true;
		zero_node.temp = -1;
		nodes.resize(nodes.size() + NN, zero_node);
		std::cout
				<< "Caution: nodes vector size increased, but new nodes were initialized with nonsense.";
	}
}

// needs L_overlap and k_spring and alpha and G_grav
// aply particle interaction forces
void sphere_forces(struct reb_simulation *const r) {
	struct reb_particle *particles = r->particles;
	for (int i = 0; i < r->N - 1; i++) {
		int j = i + 1;
		double mii = particles[i].m;
		double mjj = particles[j].m;
		double dd = fabs(particles[i].y - particles[j].y);
		if (dd < L_overlap) { // only apply if particles are close together
			double facc = k_spring * pow(fabs(dd - L_overlap), 1.5); // here is our spring force!
			particles[i].ay -= facc / mii;
			r->particles[j].ay += facc / mjj;
			double ddv = particles[i].vy - particles[j].vy;
			particles[i].ay -= alpha * ddv / mii;
			r->particles[j].ay += alpha * ddv / mjj; // damping
		}

	}
	for (int i = 0; i < r->N; i++) {
		particles[i].ay -= G_grav;  // gravity acceleration
	}
	particles[0].ay = 0.0; // don't let this guy move?
}

// give base particle a some kind of pulse?
// and fix it if outside of pulse
void pulse_base(struct reb_simulation *const r) {
	struct reb_particle *particles = r->particles;

	if (r->t < t_damp) {
		particles[0].y = py0;   // prevent first particle from going too low
		particles[0].vy = 0.0;
		particles[0].ay = 0.0;
		return;
	}

	double dtau = r->t - t_damp;
	double freq = 2.0 * M_PI / tau_pulse;

// apply a pulse then wait then apply another pulse etc
	int nps = (int) (dtau / (tau_pulse + tau_wait));
	dtau -= nps * (tau_pulse + tau_wait);

	if (dtau < tau_pulse) {
		double bpos = A_pulse * (1.0 - cos(dtau * freq)) / 2.0; // goes from 0 to Apulse at dtau=tau/2
		double bposdot = A_pulse * freq * sin(dtau * freq) / 2.0;
		// double bposddot = A_pulse*freq*freq*cos(dtau*freq)/2.0;
		particles[0].y = py0 + bpos;
		particles[0].vy = bposdot;
		// particles[0].ay = bposddot;
		particles[0].ay = 0.0;
		return;
	}
	if (particles[0].y < py0) {   // prevent first particle from going too low
		particles[0].y = py0;
		particles[0].vy = 0.0;
		particles[0].ay = 0.0;
	}

}
