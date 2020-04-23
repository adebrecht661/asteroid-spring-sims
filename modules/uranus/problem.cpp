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
#include "orb.h"

using std::string;
using std::vector;

int num_springs;
vector<spring> springs;
void reb_springs(struct reb_simulation *const r); // to pass springs to display

double def_gamma; // for gamma  of all springs
double t_damp;    // end faster damping, relaxation
double print_interval;   // for table printout
string fileroot;   // output files
int num_perts = 0; // number of point masses

int icentral = -1; // central mass location

#define NPMAX 10  // maximum number of point masses
double itaua[NPMAX], itaue[NPMAX]; // inverse of migration timescales
double itmig[NPMAX];  // inverse timescale to get rid of migration

void heartbeat(struct reb_simulation *const r);

void additional_forces(struct reb_simulation *r) {
	spring_forces(r); // spring forces
}

// check for encounters between planets
void ck_encounter(struct reb_simulation *const r);

int main(int argc, char *argv[]) {
	struct reb_simulation *const r = reb_create_simulation();
	struct spring spring_mush; // spring parameters for mush
	// Setup constants
	r->integrator = reb_simulation::REB_INTEGRATOR_LEAPFROG;
	r->gravity = reb_simulation::REB_GRAVITY_BASIC;
	r->boundary = reb_simulation::REB_BOUNDARY_NONE;
	r->G = 1;
	r->additional_forces = additional_forces; // setup callback function for additional forces
	double mball = 1.0;          // total mass of ball
	double rball = 1.0;          // radius of a ball
	double tmax = 0.0;  // if 0 integrate forever

// things to set!  can be read in with parameter file
	double dt;
	double b_distance, omegaz, ks, mush_fac, gamma_fac;
	double ratio1, ratio2, obliquity_deg;
	int lattice_type;
	double rad[NPMAX], mp[NPMAX];
	double aa[NPMAX], ee[NPMAX], ii[NPMAX];
	double longnode[NPMAX], argperi[NPMAX], meananom[NPMAX];
	int npointmass = 0;

	if (argc == 1) {
		fileroot = "t1";   // to make output files
		dt = 1e-3;    // Timestep
		tmax = 0.0;     // max integration time
		print_interval = 1.0;     // printouts for table
		lattice_type = 0;    // 0=rand 1=hcp
		b_distance = 0.15;    // min separation between particles
		mush_fac = 2.3; // ratio of smallest spring distance to minimum interparticle dist
		ks = 8e-2;   // spring constant
		// spring damping
		def_gamma = 1.0;    // final damping coeff
		gamma_fac = 5.0;    // factor initial gamma is higher that gamma_all
		t_damp = 1.0;    // gamma from initial gamma
						 // to gamma_all for all springs at this time
		ratio1 = 0.7;    // axis ratio resolved body  y/x  b/a
		ratio2 = 0.5;    // axis ratio c/a
		omegaz = 0.2;    // initial spin
		obliquity_deg = 0.0;  // obliquity

		npointmass = 1;  // add one point mass
		int ip = 0;   // index
		mp[ip] = 1.0;   // mass
		rad[ip] = 0.0;   // display radius
		itaua[ip] = 0.0;   // inverse drift rate in a
		itaue[ip] = 0.0;   // inverse drift rate in e
		itmig[ip] = 0.0;   // get rid of drift rate in inverse of this time
		// orbit
		aa[ip] = 7.0;   // distance of m1 from resolved body, semi-major orbital
		ee[ip] = 0.0;    // initial eccentricity
		ii[ip] = 0.0;    // initial inclination
		argperi[ip] = 0.0;    // initial orbtal elements
		longnode[ip] = 0.0;    // initial
		meananom[ip] = 0.0;    // initial

	} else {
		FILE *fpi;
		fpi = fopen(argv[1], "r");
		char line[300];
		char froot[100];
		fgets(line, 300, fpi);
		sscanf(line, "%s", froot);   // fileroot for outputs
		fileroot = froot;
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &dt);    // timestep
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &tmax);  // integrate to this time
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &print_interval); // output timestep
		fgets(line, 300, fpi);
		sscanf(line, "%d", &lattice_type); // particle lattice type
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &b_distance); // min interparticle distance
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &mush_fac);   // sets max spring length
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &ks);         // spring constant
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &gamma_fac);  // factor initial gamma is higher
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &def_gamma);  // damping final
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &t_damp);     // time to switch
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &ratio1);     // axis ratio for body b/a
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &ratio2);     // second axis ratio   c/a
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &omegaz);     // initial body spin
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &obliquity_deg); // obliquity degrees
		fgets(line, 300, fpi);
		sscanf(line, "%d", &npointmass); // number of point masses
		for (int ip = 0; ip < npointmass; ip++) {
			fgets(line, 300, fpi);
			sscanf(line, "%lf %lf %lf %lf %lf", mp + ip, rad + ip, itaua + ip,
					itaue + ip, itmig + ip);
			fgets(line, 300, fpi);
			sscanf(line, "%lf %lf %lf %lf %lf %lf", aa + ip, ee + ip, ii + ip,
					longnode + ip, argperi + ip, meananom + ip);
		}

	}
	double obliquity = obliquity_deg * M_PI / 180.0;   // in radians
	num_perts = 0;

/// end parameteres things to set /////////////////////////

	r->dt = dt; // set integration timestep
	const double boxsize = 3.2 * rball;    // display
	reb_configure_box(r, boxsize, 1, 1, 1);
	r->softening = b_distance / 100.0;	// Gravitational softening length

	// properties of springs
	spring_mush.gamma = gamma_fac * def_gamma; // initial damping coefficient
	spring_mush.k = ks; // spring constant
	// spring_mush.smax      = 1e6; // not used currently
	double mush_distance = b_distance * mush_fac;
	// distance for connecting and reconnecting springs

	FILE *fpr;
	char fname[200];
	sprintf(fname, "%s_run.txt", fileroot.c_str());
	fpr = fopen(fname, "w");

	num_springs = 0; // start with no springs

// for volume to be the same, adjusting here!!!!
	double volume_ratio = pow(rball, 3.0) * ratio1 * ratio2; // neglecting 4pi/3 factor
	double vol_radius = pow(volume_ratio, 1.0 / 3.0);

	rball /= vol_radius; // volume radius used to compute body semi-major axis
	// assuming that semi-major axis is rball
	// fprintf(fpr,"vol_ratio %.6f\n",volume_ratio); // with respect to 4pi/3
	fprintf(fpr, "a %.3f\n", rball);
	fprintf(fpr, "b %.3f\n", rball * ratio1);
	fprintf(fpr, "c %.3f\n", rball * ratio2);
	volume_ratio = pow(rball, 3.0) * ratio1 * ratio2; // neglecting 4pi/3 factor
	fprintf(fpr, "vol_ratio %.6f\n", volume_ratio); // with respect to 4pi/3
	// so I can check that it is set to 1

	// create resolved body particle distribution
	if (lattice_type == 0) {
		// rand_football_from_sphere(r,b_distance,rball,rball*ratio1, rball*ratio2,mball );
		rand_ellipsoid(r, b_distance, rball, rball * ratio1, rball * ratio2,
				mball);
	}
	if (lattice_type == 1) {
		hcp_ellipsoid(r, b_distance, rball, rball * ratio1, rball * ratio2,
				mball);
	}
	if (lattice_type == 2) {
		cubic_ellipsoid(r, b_distance, rball, rball * ratio1, rball * ratio2,
				mball);
	}

	int il = 0;  // index range for resolved body
	int ih = r->N;

	subtract_com(r, il, ih);  // move reference frame to resolved body
	// spin it
	subtract_cov(r, il, ih); // center of velocity subtracted
	spin_body(r, il, ih, Vector( { 0.0, 0.0, omegaz })); // change one of these zeros to tilt it!
	// can spin about non principal axis
	subtract_cov(r, il, ih); // center of velocity subtracted
	double speriod = fabs(2.0 * M_PI / omegaz);
	printf("spin period %.6f\n", speriod);
	fprintf(fpr, "spin period %.6f\n", speriod);
	double E_gamma = 0.0; // Euler angles
	double E_beta = obliquity;
	if (npointmass > 0) {
		E_gamma = longnode[0]; // relative to orbit
		E_beta += ii[0];
		// printf("hello %.2f %.2f %.2f\n",obliquity, E_gamma,E_beta);
	}
	rotate_body(r, il, ih, 0.0, E_beta, E_gamma); // tilt by obliquity

	// make springs, all pairs connected within interparticle distance mush_distance
	connect_springs_dist(r, mush_distance, 0, r->N, spring_mush);

	// assume minor semi is rball*ratio2
	double ddr = rball * ratio2 - 0.5 * mush_distance;
	ddr = 0.4; // mini radius  for computing Young modulus
	double Emush = Young_mesh(r, il, ih, 0.0, ddr); // compute from springs out to ddr
	double Emush_big = Young_full_mesh();
	printf("ddr = %.3f mush_distance =%.3f \n", ddr, mush_distance);
	printf("Young's modulus %.6f\n", Emush);
	printf("Young's modulus big %.6f\n", Emush_big);
	fprintf(fpr, "Young's_modulus %.6f\n", Emush);
	fprintf(fpr, "Young's_modulus big %.6f\n", Emush_big);
	fprintf(fpr, "mush_distance %.4f\n", mush_distance);
	double LL = mean_spring_length();  // mean spring length
	printf("mean L = %.4f\n", LL);
	fprintf(fpr, "mean_L  %.4f\n", LL);
	// relaxation timescale
	// note no 2.5 here!
	double tau_relax = 1.0 * def_gamma * 0.5 * (mball / (r->N - num_perts))
			/ spring_mush.k; // Kelvin Voigt relaxation time
	printf("relaxation time %.3e\n", tau_relax);
	fprintf(fpr, "relaxation_time  %.3e\n", tau_relax);

	double om = 0.0;
	if (npointmass > 0) {
		// set up rest of point masses
		for (int ip = 0; ip < npointmass; ip++) {
			OrbitalElements orb_el;
			orb_el.a = aa[ip];
			orb_el.e = ee[ip];
			orb_el.i = ii[ip];
			orb_el.long_asc_node = longnode[ip];
			orb_el.arg_peri = argperi[ip];
			orb_el.mean_anom = meananom[ip];
			om = add_pt_mass_kep(r, il, ih, icentral, mp[ip], rad[ip], orb_el);
			fprintf(fpr, "resbody mm=%.3f period=%.2f\n", om, 2.0 * M_PI / om);
			printf("resbody mm=%.3f period=%.2f\n", om, 2.0 * M_PI / om);
			// note central body mass (only changes on first run - can probably do this better
			icentral = ih;
		}
		num_perts = npointmass;
	}

	// this is probably not correct if obliquity is greater than pi/2
	double barchi = 2.0 * fabs(om - omegaz) * tau_relax; // initial value of barchi
	double posc = 0.5 * 2.0 * M_PI / fabs(om - omegaz); // for oscillations!
	fprintf(fpr, "barchi  %.4f\n", barchi);
	printf("barchi %.4f\n", barchi);
	fprintf(fpr, "posc %.6f\n", posc);

// double na = om*aa;
// double adot = 3.0*m1*na/pow(aa,5.0); // should approximately be adot
// fprintf(fpr,"adot %.3e\n",adot);

	// ratio of numbers of particles to numbers of springs for resolved body
	double Nratio = (double) num_springs / (ih - il);
	printf("N=%d  NS=%d NS/N=%.1f\n", r->N, num_springs, Nratio);
	fprintf(fpr, "N=%d  NS=%d NS/N=%.1f\n", r->N, num_springs, Nratio);
	fclose(fpr);

	reb_springs(r); // pass spring index list to display
	r->heartbeat = heartbeat;
	center_sim(r, il, ih);  // move reference frame to resolved body

	// max integration time
	if (tmax == 0.0)
		reb_integrate(r, INFINITY);
	else
		reb_integrate(r, tmax);
}

#define NSPACE 50
void heartbeat(struct reb_simulation *const r) {
	static int first = 0;
	static char extendedfile[50];
	static char pointmassfile[NPMAX * NSPACE];
	if (first == 0) {
		first = 1;
		sprintf(extendedfile, "%s_ext.txt", fileroot.c_str());
		for (int i = 0; i < num_perts; i++) {
			sprintf(pointmassfile + i * NSPACE, "%s_pm%d.txt", fileroot.c_str(), i);
		}
	}
	if (reb_output_check(r, 10.0 * r->dt)) {
		reb_output_timing(r, 0);
	}
	if (fabs(r->t - t_damp) < 0.9 * r->dt)
		set_gamma(def_gamma);
	// damp initial bounce only
	// reset gamma only at t near t_damp

	// stuff to do every timestep
	center_sim(r, 0, r->N - num_perts); // move reference frame to resolved body for display

	ck_encounter(r); // check for collisions

	// dodrifts!!!!
	int il = 0; // resolved body index range
	int ih = r->N - num_perts;
	if (num_perts > 0) {
		for (int i = 0; i < num_perts; i++) {
			int ip = icentral + i;  // which body drifting
			double migfac = exp(-1.0 * r->t * itmig[i]);
			if (i == 0)  // it is central mass, so drift resolved body
				drift_resolved(r, r->dt, itaua[i] * migfac, itaue[i] * migfac,
						icentral, il, ih);
			else
				// it is another point mass, drifts w.r.t to icentral
				drift_bin(r, r->dt, itaua[i] * migfac, itaue[i] * migfac,
						icentral, ip);
		}

	}

	if (reb_output_check(r, print_interval)) {
		write_resolved_no_E(r, il, ih, extendedfile); // orbital info and stuff
		if (num_perts > 0)
			for (int i = 0; i < num_perts; i++) {
				int ip = icentral + i;
				write_pt_mass(r, ip, i, pointmassfile + i * NSPACE);
			}
	}

}

// make a spring index list for display
void reb_springs(struct reb_simulation *const r) {
	r->NS = num_springs;
	r->springs_i = (int*) malloc(num_springs * sizeof(int));
	r->springs_j = (int*) malloc(num_springs * sizeof(int));
	for (int i = 0; i < num_springs; i++) {
		r->springs_i[i] = springs[i].particle_1;
		r->springs_j[i] = springs[i].particle_2;
	}
}

// check for encounters between planets and between planets and extended body
void ck_encounter(struct reb_simulation *const r) {
	if (num_perts > 1) {
		// check for encounter between a planet or star and extended body
		double xc = 0.0;
		double yc = 0.0;
		double zc = 0.0;
		// compute_com(r,0, icentral, &xc, &yc, &zc); // not necessary as centerbody called beforhand
		for (int i = 0; i < num_perts; i++) { // includes central star
			int ip_i = icentral + i;
			double x2 = pow(r->particles[ip_i].x - xc, 2.0);
			double y2 = pow(r->particles[ip_i].y - yc, 2.0);
			double z2 = pow(r->particles[ip_i].z - zc, 2.0);
			double dist = sqrt(x2 + y2 + z2);
			double sumrad = r->particles[ip_i].r + 1.0;
			if (dist < sumrad) { // smash!!!
				printf("\n collision %d + extended\n", i);
				reb_exit("collision");
			}
		}
		// check for encounter between planets (including star)
		for (int i = 0; i < num_perts; i++) { // includes central star
			int ip_i = icentral + i;
			for (int j = i + 1; j < num_perts; j++) {
				int ip_j = icentral + j;
				double x2 = pow(r->particles[ip_i].x - r->particles[ip_j].x,
						2.0);
				double y2 = pow(r->particles[ip_i].y - r->particles[ip_j].y,
						2.0);
				double z2 = pow(r->particles[ip_i].z - r->particles[ip_j].z,
						2.0);
				double dist = sqrt(x2 + y2 + z2);
				double sumrad = r->particles[ip_i].r + r->particles[ip_j].r;
				if (dist < sumrad) { // smash!!!
					printf("\n collision %d %d\n", ip_i, ip_j);
					reb_exit("collision");
				}
			}
		}
	}
}
