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
 * This lets you read in the Bennu shape model and give it an impact
 * no external point mass particles in this routine
 */

#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "libconfig.h++"
extern "C" {
#include "rebound.h"
#include "gravity.h"
}
#include "matrix_math.h"
#include "input_spring.h"
#include "output_spring.h"
#include "springs.h"
#include "physics.h"
#include "stress.h"
#include "shapes.h"
#include "heat.h"
#include "orb.h"

using namespace libconfig;

using std::string;
using std::vector;

//void reb_calculate_acceleration(reb_simulation *r); // rebound routine

int num_springs;  // number of springs
spring *springs; // global array containing springs
void reb_springs(reb_simulation *const r); // to pass springs to display so springs can be seen in viewer
int *surfp; // is allocated with marksurface routine! surface particles!

double gamma_all;
double gamma_fac; // for adjustment of gamma of all springs
double t_damp;    // end faster damping, relaxation
int ndt_print;   // for printouts
char froot[30];   // output files
string fileroot;
int num_perts; // number of point mass perturbers, not used here

// globals for a pulse
void addforcepulse(reb_simulation *const r, double lat, double longi,
		double dangle, double famp, double taus, int *surfarr); // do a pulse on the surface particles
double lat_psh;    // latitute in radius of pulse
double long_psh;   // longitude in radius of pulse
double dangle_psh; // angular distance in radius of pulse applied
double famp_psh;   // force amplitude of pulse
double dt_psh;     // length of pulse in time and starting at t_damp

// Global scales
double mass_scale, time_scale, length_scale, temp_scale, omega_scale, vel_scale,
		p_scale, L_scale, a_scale, F_scale, E_scale, dEdt_scale, P_scale;

void reset_surface_radii(reb_simulation *const r, double b_distance); // post impact reset of surface

void heartbeat(reb_simulation *const r);

void additional_forces(reb_simulation *r) {
	spring_forces(r); // spring forces
	if ((r->t > t_damp) && (r->t <= t_damp + dt_psh)) {
		addforcepulse(r, lat_psh, long_psh, dangle_psh, famp_psh, dt_psh,
				surfp); // pulse on surface
	}
}

int main(int argc, char *argv[]) {
	reb_simulation *const r = reb_create_simulation();
	spring spring_mush; // spring parameters for mush
	// Setup constants
	r->integrator = reb_simulation::REB_INTEGRATOR_LEAPFROG;
	r->gravity = reb_simulation::REB_GRAVITY_BASIC;
	r->boundary = reb_simulation::REB_BOUNDARY_NONE;
	r->G = 1;
	r->additional_forces = additional_forces; // setup callback function for additional forces
	double mball = 1.0;          // total mass of ball
	double tmax = 0.0;  // if 0 integrate forever

// things to set!  can be read in with parameter file
	double dt, b_distance, omegaz, ks, mush_fac;
	double lat_psh_deg = 0, long_psh_deg = 0, dangle_psh_deg = 0;
	double obliq_deg = 0, surfdist = 0;
	int surface_type;
	num_perts = 0;
	double ratio1 = 1.0;
	double ratio2 = 1.0;

	if (argc == 1) {
		fileroot = "t1";    // to make output files
		dt = 2e-3;    // Timestep
		tmax = 0;              // max integration time, if 0 then infinity
		ndt_print = 20;        // printouts each one this times dt
		surface_type = 0;      // 0=bennu shape mode, 1=football,2=cone
		b_distance = 0.090;   // min separation between particles
		mush_fac = 2.3; // ratio of smallest spring distance to minimum interparticle dist
		omegaz = 0.7;     // initial spin
		obliq_deg = 0.0;
		// spring damping, gamma in units of 1/time
		gamma_fac = 2.0; // initial factor for initial damping value for springs
		gamma_all = 1.0e-2;  // final damping coeff
		t_damp = 1e-3;    // gamma to final values for all springs at this time
		ks = 10.0;   // spring constant
		surfdist = 0.15;       // for identifying surface particles
		lat_psh_deg = 0.0;    // latitude of pulse in degrees
		long_psh_deg = 0.0;    // longitude of pulse in degrees
		dangle_psh_deg = 20.0;  // angular width distance of pulse in degrees
		famp_psh = 1.0e1;     // amplitude total of pulse force
		dt_psh = 1e-1;        // length of pulse

	} else {
		FILE *fpi;
		fpi = fopen(argv[1], "r");
		char line[300];
		fgets(line, 300, fpi);
		sscanf(line, "%s", froot); // file root for outputs
		fileroot = froot;
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &dt);  // timestep
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &tmax); // max time for integration
		fgets(line, 300, fpi);
		sscanf(line, "%d", &ndt_print);  // between printouts in timesteps
		fgets(line, 300, fpi);
		sscanf(line, "%d", &surface_type);  // shape model 0=bennu
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &b_distance);  // min interparticle spacing
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &mush_fac);   // sets max spring length
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &omegaz);    // spin
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &obliq_deg);  // obliquity or tilt
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &gamma_fac); // multiply to get initial damping relaxation
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &gamma_all);  // sets final damping parameter
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &t_damp);    // time for initial relaxation
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &ks);        // spring constant
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &surfdist);  // for finding surface
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &lat_psh_deg);   // latitude of pulse (degrees)
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &long_psh_deg);  // longitude of pulse (degrees)
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &dangle_psh_deg); // angular width of pulse (degrees)
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &famp_psh);    // amplitude of pulse
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &dt_psh);      // length of time of pulse
		fgets(line, 300, fpi);
		sscanf(line, "%lf %lf", &ratio1, &ratio2); // used if not using bennu shape model
	}

/// end of things to set /////////////////////////

	r->dt = dt; // set integration timestep
	r->softening = b_distance / 100.0;	// Gravitational softening length

// degrees to radians
	double obliquity = obliq_deg * M_PI / 180.0;
	lat_psh = lat_psh_deg * M_PI / 180.0;
	long_psh = long_psh_deg * M_PI / 180.0;
	dangle_psh = dangle_psh_deg * M_PI / 180.0;

	// properties of springs
	spring_mush.gamma = gamma_fac * gamma_all; // initial damping coefficient
	spring_mush.k = ks; // spring constant
	double mush_distance = b_distance * mush_fac;
	// distance for connecting and reconnecting springs
	printf("max spring length =%.2f min interparticle distance=%.2f\n",
			mush_distance, b_distance);

	FILE *fpr;
	char fname[200];
	sprintf(fname, "%s_run.txt", froot);
	fpr = fopen(fname, "w");

	num_springs = 0; // start with no springs

	if (surface_type == 0) { // using bennu shape model here
		// read in shape file
		string filename = "101955bennu.tab"; // the shape model!
		// read in vertex file
		read_vertices(r, filename);
		int N_bennu = r->N;
		std::cout << filename << " read in \n";
		// correct from km to mean radius.   mean radius 246m is  0.246km
		stretch(r, 0, r->N, 1.0 / 0.246);
		// fill shape with particles
		rand_shape(r, b_distance, 1.0); // fill in particles into the shape model
		printf("shape filled \n");

		// find surface particles
		mark_surf_shrink_int_shape(r, 0, N_bennu, surfdist); // number here is distance from nearest surface particle
		// this is done with shape model in place
		rm_particles(r, 0, N_bennu); // remove shape vertices
		printf("shape vertices deleted \n");
	}
	if (surface_type == 1) {  //random football shape model
		double rball = 1.0;
		double volume_ratio = pow(rball, 3.0) * ratio1 * ratio2; // neglecting 4pi/3 factor
		//  vol of a triax ellipsoid is abc*4pi/3
		double vol_radius = pow(volume_ratio, 1.0 / 3.0);
		rball /= vol_radius; // volume radius used to compute body semi-major axis
		// assuming that semi-major axis is rball
		rand_ellipsoid(r, b_distance, rball, rball * ratio1, rball * ratio2,
				mball);
		mark_surf_shrink_int_ellipsoid(r, surfdist, rball, rball * ratio1,
				rball * ratio2);
	}
	if (surface_type == 2) {  //random cone shape model
		double rball = 1.0; // h = r*ratio1
		double volume_ratio = pow(rball, 3.0) * ratio1 * 0.5; // vol of a cone is 2 pi/3 r^2 h
		double vol_radius = pow(volume_ratio, 1.0 / 3.0);
		rball /= vol_radius; // volume radius used to compute r,h of cone
		rand_cone(r, b_distance, rball, rball * ratio1, mball); // h is now ratio1*rball
		mark_surf_shrink_int_cone(r, surfdist, rball, rball * ratio1);
		// exit(0);
	}

	int il = 0;
	int ih = r->N;
	subtract_com(r, il, ih);  // move reference frame to resolved body
	subtract_cov(r, il, ih);  // subtract center of velocity
	printf("centerbody, subtractcov\n");

	// spin it
	spin_body(r, il, ih, Vector({0.0, 0.0, omegaz}));  // you can change one of these to tilt!
	subtract_cov(r, il, ih);  // subtract center of velocity
	double speriod = fabs(2.0 * M_PI / omegaz);
	printf("spin period %.6f\n", speriod);
	fprintf(fpr, "spin period %.6f\n", speriod);
	fprintf(fpr, "omegaz %.6f\n", omegaz);

	if (obliquity != 0.0)
		rotate_body(r, il, ih, 0.0, obliquity, 0.0); // tilt by obliquity in radians

	// make springs, all pairs connected within interparticle distance mush_distance
	connect_springs_dist(r, mush_distance, 0, r->N, spring_mush);
	printf("springs connected \n");

	// double minr = min_radius(r,il,ih);
	double maxr = max_radius(r, il, ih);
	// printf("minr=%lf maxr=%lf\n",minr,maxr); // min and max radii of particles
	const double boxsize = 2.5 * maxr;    // display
	reb_configure_box(r, boxsize, 1, 1, 1); // sets size of display view

	double ddr = 0.6 * maxr;
	double Emush = Young_mesh(r, il, ih, 0.0, ddr); // compute Youngs modulus
	printf("Young's modulus %.6f\n", Emush);
	fprintf(fpr, "Young's_modulus %.6f\n", Emush);
	printf("ddr = %.3f mush_distance =%.3f \n", ddr, mush_distance);
	fprintf(fpr, "max spring length %.4f\n", mush_distance);
	double LL = mean_spring_length();
	printf("mean L = %.4f\n", LL); // mean spring length
	fprintf(fpr, "mean_L  %.4f\n", LL);

	// factor of 0.5 is due to reduced mass being used in calculation
	double tau_relax = 1.0 * gamma_all * 0.5 * (mball / (r->N - 1))
			/ spring_mush.k; // Kelvin Voigt relaxation time
	printf("relaxation time %.3e\n", tau_relax);
	fprintf(fpr, "relaxation_time  %.3e\n", tau_relax);

	double Nratio = (double) num_springs / (double) r->N;
	printf("N=%d  NS=%d NS/N=%.1f\n", r->N, num_springs, Nratio);
	fprintf(fpr, "N=%d  NS=%d NS/N=%.1f\n", r->N, num_springs, Nratio);
	int Nsurf = 0;
	for (int i = 0; i < r->N; i++)
		Nsurf += surfp[i];
	printf("Nsurf=%d  \n", Nsurf);
	fprintf(fpr, "Nsurf=%d  \n", Nsurf); // number of surface particles

	fclose(fpr);

	reb_springs(r); // pass spring index list to display
	r->heartbeat = heartbeat; // set up integration

	// calculate surface gravity alone at beginning of simulation
	reb_calculate_acceleration(r);
	string filename = fileroot + "_surf_nosprings.txt";
	write_surf_part(r, 0, r->N - num_perts, filename); // in case you want to know about surface

	if (tmax == 0.0) // integrate!!!!!
		reb_integrate(r, INFINITY);
	else
		reb_integrate(r, tmax);
}

// this is stuff done during the simulation integration
void heartbeat(reb_simulation *const r) {
	static int ifile = 0;
	if (reb_output_check(r, 10.0 * r->dt)) {
		reb_output_timing(r, 0);
	}
	if (fabs(r->t - t_damp) < 0.9 * r->dt)
		set_gamma(gamma_all);
	// damp initial bounce only
	// reset gamma only at t near t_damp
	// divide all gamma's by gamma_fac

	// printouts!
	if (r->t > t_damp) {
		if (reb_output_check(r, ndt_print * r->dt)) {
			string fname = froot + zero_pad_int(5, ifile);
			write_surf_part(r, 0, r->N - num_perts, fname);
			ifile++;
		}
	}

	subtract_com(r, 0, r->N - num_perts);  // move reference frame, position only

}

// make a spring index list, to pass to display
void reb_springs(reb_simulation *const r) {
	r->NS = num_springs;
	r->springs_i = (int*) malloc(num_springs * sizeof(int));
	r->springs_j = (int*) malloc(num_springs * sizeof(int));
	for (int i = 0; i < num_springs; i++) {
		r->springs_i[i] = springs[i].particle_1;
		r->springs_j[i] = springs[i].particle_2;
	}
}

// do an impact by pushing surface particles
//   push on surface particles radially,
//   push is centered at lat,long,  in radians
//   all particles on surface with angular size dangle in radians
// with total force amplitude famp? distributed among all particles int the angular distance
// push direction is radially inward, though if famp<0 you could make it outward
// this routine is meant to be called for a specific length of time
// this routine changes surface particle accelerations so is an additional force
// if famp is constant then you get a tophat force pulse
// here we have a cosine function for the pulse
void addforcepulse(reb_simulation *const r, double lat, double longi,
		double dangle, double famp, double taus, int *surfarr) {
	static int first = 0;
	static double t0 = 0.0;
	static double r0 = 0.0; // store the initial particle radius for display
	if (first == 0) {
		first = 1;
		t0 = r->t; // keep track of the start of the pulse
		for (int i = 0; i < r->N - num_perts; i++) {
			if (surfarr[i] == 1) {
				r0 = r->particles[i].r; // keep track of radius of surface particles for display
				i = r->N;
			}
		}
	}
	double x0 = cos(longi) * cos(lat);
	double y0 = sin(longi) * cos(lat);
	double z0 = sin(lat);
	double tpulse = r->t - t0;  // time since start of pulse

	int npush = 0;
	for (int i = 0; i < r->N - num_perts; i++) { // first figure out how many particles we will push
		if (surfarr[i] == 1) { // particle is on the surface
			double x = r->particles[i].x;
			double y = r->particles[i].y;
			double z = r->particles[i].z;
			double rad = sqrt(x * x + y * y + z * z);
			double rhatdotr0 = (x * x0 + y * y0 + z * z0) / rad;
			double angdist = acos(rhatdotr0); // in [0,pi] is positive
			if (angdist < dangle) { // yes we will push the particle!
				npush++; // number of particles that are pushed
			}
		}
	}
	printf("addforcepulse: npush=%d\n", npush);
	for (int i = 0; i < r->N - num_perts; i++) {
		if (surfarr[i] == 1) { // particle is on the surface
			double x = r->particles[i].x;
			double y = r->particles[i].y;
			double z = r->particles[i].z;
			double m = r->particles[i].m;
			double rad = sqrt(x * x + y * y + z * z);
			double rhatdotr0 = (x * x0 + y * y0 + z * z0) / rad;
			double angdist = acos(rhatdotr0); // in [0,pi] is positive
			double fac = famp / (m * npush * rad); // dv = F*dt/m  with total force applied F = famp
			// time dependence of pulse here
			fac *= 1.0 - cos(2.0 * M_PI * tpulse / taus); // smooth the pulse shape
			if (angdist < dangle) { // we now push the surface particle!
				r->particles[i].ax -= x * fac; // radial push, accellerations are changed
				r->particles[i].ay -= y * fac;
				r->particles[i].az -= z * fac;
				if (tpulse < taus - 2.1 * r->dt)
					r->particles[i].r = 0.09; // something weird so we can see where we pushed
				else
					r->particles[i].r = r0; // return to normal at the end of the pulse
			}
		}
	}
	// subtractcov(r,0,r->N);  // subtract center of velocity
	if (tpulse >= taus - 2.1 * r->dt) {
		for (int i = 0; i < r->N - num_perts; i++) {
			if (surfarr[i] == 1) // particle is on the surface
				r->particles[i].r = r0; // return to normal at the end of the pulse
		}
	}
}

// reset surface particle radii post impact
void reset_surface_radii(reb_simulation *const r, double b_distance) {
	for (int i = 0; i < r->N; i++) {
		if (surfp[i] == 1)
			r->particles[i].r = b_distance / 4;
	}
}
