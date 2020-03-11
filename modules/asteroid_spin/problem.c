/**
 * resolved mass spring model
 * using the leap frog integrator. 
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "rebound.h"
#include "tools.h"
#include "output.h"
#include "spring.h"

int NS;  // global numbers of springs
struct spring *springs;
void reb_springs();  // to pass springs to display
void read_particles_i();

double gamma_all; // for gamma  of all springs
double t_damp;    // end faster damping, relaxation
double t_print;   // for table printout 
char froot[30];   // output files
int npert = 0;

// passed to rebound so you can have forces other than gravity
void additional_forces(struct reb_simulation *r) {
	zero_accel(r); // because gravity routines in rebound set acceleration and if you don't
				   // call a gravity routine you need to zero the acceleration array
	spring_forces(r); // springs

}

void heartbeat(struct reb_simulation *const r);

int main(int argc, char *argv[]) {
	struct reb_simulation *const r = reb_create_simulation();
	struct spring spring_mesh; // spring parameters for mush
	// Setup constants
	r->integrator = REB_INTEGRATOR_LEAPFROG;
	// r->gravity	= REB_GRAVITY_BASIC;
	r->gravity = REB_GRAVITY_NONE;
	r->boundary = REB_BOUNDARY_NONE;
	r->G = 1;
	r->additional_forces = additional_forces; // spring forces
	// setup callback function for additional forces
	// double mball = 1.0;          // total mass of ball
	double rball = 1.0;          // radius of a ball
	double tmax = 0.0;  // if 0 integrate forever

// things to set!  can be read in with parameter file
	double dt = 0.000;
	double omegaz = 0.0;
	double ks = 0.0;
	double gamma_fac = 0;
	double mesh_distance = 0.0; //  max distance for connecting springs

	if (argc == 1) {
		strcpy(froot, "a2");   // to make output files
		dt = 1e-2;    // Timestep
		tmax = 0.0;     // max integration time
		t_print = 100.0;     // printouts
		ks = 0.005;   // spring constant
		// spring damping
		gamma_all = 1.0;    // final damping coeff
		gamma_fac = 5.0;    // factor initial gamma is higher that gamma_all
		t_damp = 1.0;    // gamma from initial gamma
						 // to gamma_all for all springs at this time
		omegaz = 0.8;    // initial spin
		mesh_distance = 0.1; //  max distance for connecting springs
	}

/// end parameters of things to set /////////////////////////

	r->dt = dt; // set integration timestep
	const double boxsize = 3.2 * rball;    // display
	reb_configure_box(r, boxsize, 1, 1, 1);
	r->softening = 1e-6;	// Gravitational softening length

	// properties of springs
	spring_mesh.gamma = gamma_fac * gamma_all; // initial damping coefficient
	spring_mesh.ks = ks; // spring constant

	FILE *fpr;
	char fname[200];
	sprintf(fname, "%s_run.txt", froot);
	fpr = fopen(fname, "w");

	NS = 0; // start with no springs

	// read in a particle distribution  here!!!!
	read_particles_i(r, "input.txt");  // read in routine below
	// format is  &x, &y, &z, &vx, &vy, &vz, &rad, &m);

	int i_low = 0;  // index range for resolved body
	int i_high = r->N;

	subtractcom(r, i_low, i_high);  // move reference frame to resolved body

	// rotate_to_principal(r, il, ih); // rotate to principal axes

	// spin it
	subtractcov(r, i_low, i_high); // center of velocity subtracted
	spin(r, i_low, i_high, 0.0, 0.0, omegaz); // change one of these zeros to tilt it!
	// can spin about non principal axis
	subtractcov(r, i_low, i_high); // center of velocity subtracted
	double speriod = fabs(2.0 * M_PI / omegaz);
	printf("spin period %.6f\n", speriod);
	fprintf(fpr, "spin period %.6f\n", speriod);

	// make springs, all pairs connected within interparticle distance mesh_distance
	connect_springs_dist(r, mesh_distance, i_low, i_high, spring_mesh);

	double ddr = 0.4; // mini radius  for computing Young modulus
	double Emush = Young_mush(r, i_low, i_high, 0.0, ddr); // compute from springs out to ddr
	double Emush_big = Young_mush_big(r, i_low, i_high);
	printf("ddr = %.3f mush_distance =%.3f\n", ddr, mesh_distance);
	printf("Youngs_modulus %.6f\n", Emush);
	printf("Youngs_modulus_big %.6f\n", Emush_big);
	fprintf(fpr, "Youngs_modulus %.6f\n", Emush);
	fprintf(fpr, "Youngs_modulus_big %.6f\n", Emush_big);
	fprintf(fpr, "mush_distance %.4f\n", mesh_distance);
	double LL = mean_L(r);  // mean spring length
	printf("mean_L  %.4f\n", LL);
	fprintf(fpr, "mean_L %.4f\n", LL);

	// ratio of numbers of particles to numbers of springs for resolved body
	double Nratio = (double) NS / (i_high - i_low);
	printf("N=%d NS=%d NS/N=%.1f\n", r->N, NS, Nratio);
	fprintf(fpr, "N %d\n", r->N);
	fprintf(fpr, "NS %d\n", NS);
	fprintf(fpr, "NS/N %.1f\n", Nratio);

	fclose(fpr);

	reb_springs(r); // pass spring index list to display
	r->heartbeat = heartbeat;
	centerbody(r, i_low, i_high);  // move reference frame to resolved body

	// max integration time
	if (tmax == 0.0)
		reb_integrate(r, INFINITY);
	else
		reb_integrate(r, tmax);
}

void heartbeat(struct reb_simulation *const r) {
	static int first = 0;
	static char extendedfile[50];
	if (first == 0) {
		first = 1;
		sprintf(extendedfile, "%s_ext.txt", froot);
	}
	if (reb_output_check(r, 10.0 * r->dt)) {
		reb_output_timing(r, 0);
	}
	if (fabs(r->t - t_damp) < 0.9 * r->dt)
		set_gamma(gamma_all, 0, r->N - npert);
	// damp initial bounce only
	// reset gamma only at t near t_damp

	// stuff to do every timestep
	centerbody(r, 0, r->N - npert); // move reference frame to resolved body for display
	if (reb_output_check(r, t_print)) {
		int index = (int) (r->t / t_print);
		if (index > 0) {
			write_springs(r, froot, index);
			write_particles(r, froot, index);
		}
	}

}

// make a spring index list for display
void reb_springs(struct reb_simulation *const r) {
	r->NS = NS;
	r->springs_ii = malloc(NS * sizeof(int));
	r->springs_jj = malloc(NS * sizeof(int));
	for (int i = 0; i < NS; i++) {
		r->springs_ii[i] = springs[i].i;
		r->springs_jj[i] = springs[i].j;
	}
}

void read_particles_i(struct reb_simulation *const r, char *filename) {
	printf("\n reading in particles %s\n", filename);
	FILE *fpi;
	fpi = fopen(filename, "r");
	char string[300];
	struct reb_particle pt;
	pt.ax = 0.0;
	pt.ay = 0.0;
	pt.az = 0.0;
	double x, y, z, vx, vy, vz, m, rad;
	while (fgets(string, 300, fpi) != NULL) {
		sscanf(string, "%lf %lf %lf %lf %lf %lf %lf %lf\n", &x, &y, &z, &vx,
				&vy, &vz, &rad, &m);
		pt.x = x;
		pt.y = y;
		pt.z = z;
		pt.vx = vx;
		pt.vy = vy;
		pt.vz = vz;
		pt.r = rad;
		pt.m = m;
		reb_add(r, pt);
	}
	fclose(fpi);
	printf("read_particles: N=%d\n", r->N);
}

