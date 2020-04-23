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
double gammafac2;
double gamma_fac;
double t_damp;    // end faster damping, relaxation
double print_interval;   // for table printout
string fileroot;   // output files
int num_perts = 0;

int icentral = -1; // central mass location

#define NPMAX 10  // maximum number of point masses
double itaua[NPMAX], itaue[NPMAX]; // inverse of migration timescales
double itmig[NPMAX];  // inverse timescale to get rid of migration

void heartbeat(struct reb_simulation *const r);

void additional_forces(struct reb_simulation *r) {
	spring_forces(r); // spring forces
}

int main(int argc, char *argv[]) {
	struct reb_simulation *const r = reb_create_simulation();
	struct spring spring_mush; // spring parameters for mush
	// Setup constants
	r->integrator = reb_simulation::REB_INTEGRATOR_LEAPFROG;
	r->gravity = reb_simulation::REB_GRAVITY_BASIC;
	r->boundary = reb_simulation::REB_BOUNDARY_NONE;
	r->G = 1.00;	//gravitational constant, should be 1
	r->additional_forces = additional_forces; // setup callback function for additional forces
	double mball = 1.0;          // total mass of ball
	double rball = 1.0;          // radius of a ball
	double tmax = 0.0;  // if 0 integrate forever

// things to set!  can be read in with parameter file
	double dt;
	double b_distance, omegax, omegay, omegaz, ks, mush_fac;
	double ratio1, ratio2, obliquity_deg;
	int lattice_type;
	double rad[NPMAX], mp[NPMAX];
	double aa[NPMAX], ee[NPMAX], ii[NPMAX];
	double longnode[NPMAX], argperi[NPMAX], meananom[NPMAX];
	int npointmass = 0;
	double kfac = 1.0;
	double Rshell = 0.0;

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
		gamma_fac = 5.0;    // factor initial gamma is higher that def_gamma
		t_damp = 1.0;    // gamma from initial gamma
						 // to def_gamma for all springs at this time
		ratio1 = 0.7;    // axis ratio resolved body  y/x  b/a
		ratio2 = 0.5;    // axis ratio c/a
		omegaz = 0.2;    // initial spin
		omegax = 0.0;    // initial spin
		omegay = 0.0;    // initial spin
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
		gammafac2 = 1.0;

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
		sscanf(line, "%lf", &omegax);     // initial body spin
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &omegay);     // initial body spin
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &omegaz);     // initial body spin
		fgets(line, 300, fpi);
		sscanf(line, "%lf %lf", &kfac, &gammafac2); // change inside Rshell  by this factor
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &Rshell);     // Rshell

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
	// double obliquity = obliquity_deg*M_PI/180.0;   // in radians
	num_perts = 0;

/// end parameteres things to set /////////////////////////

	r->dt = dt; // set integration timestep
	const double boxsize = 3.8 * rball;    // display
	reb_configure_box(r, boxsize, 1, 1, 1);
	// r->softening      = b_distance/100.0;	// Gravitational softening length
	r->softening = 0.0;

	// properties of springs
	spring_mush.gamma = gamma_fac * def_gamma; // initial damping coefficient
	spring_mush.k = ks; // spring constant
	// spring_mush.smax      = 1e6; // not used currently
	double mush_distance = b_distance * mush_fac;
	// distance for connecting and reconnecting springs

	FILE *fpr;
	char fname[200];
	sprintf(fname, "%s_run.txt", fileroot.c_str());  // store parameters for simulation
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

	// if ratio1 ==  1.0 then oblate
	// both ratios should be <1
	// if ratio1 == ratio2 then prolate

	// create resolved body particle distribution
// the body is orientated with principal axes along xyz axes
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
	if (lattice_type == 3) { // use Bennu  shape model
		string filename;
		filename = "101955bennu.tab"; // the shape model!
		// read in vertex file
		read_vertices(r, filename);
		int N_bennu = r->N; // needed to remove shape vertices
		std::cout << filename << " read in\n";
		// correct from km to mean radius.   mean radius 246m is  0.246km
		stretch(r, 0, r->N, 1.0 / 0.246);
		// fill shape with particles
		rand_shape(r, b_distance, 1.0); // fill in particles into the shape model
		printf("shape filled \n");
		rm_particles(r, 0, N_bennu); // remove shape vertices
		printf("shape vertices deleted \n");

	}

// for oblate z is up is narrow symmetry axis, xy  in plane is circular
// for prolate x is long symmetry axis in the plane, z is still up

	int il = 0;  // index range for resolved body
	int ih = r->N;

	subtract_com(r, il, ih);  // move reference frame to resolved body
	subtract_cov(r, il, ih); // center of velocity subtracted
	// spin it
	// this routine spins the body with Omega three components
	spin_body(r, il, ih, Vector( { omegax, omegay, omegaz }));
	// can spin about non principal axis
	subtract_cov(r, il, ih); // center of velocity subtracted

	// rotate_body(r, il, ih, 0.0, -atan2(llz,lly),0);
	// rotate by Euler angles, including velocities

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
	printf("mean spring length = %.4f\n", LL);
	fprintf(fpr, "mean_L (mean spring length)  %.4f\n", LL);
	// relaxation timescale
	// note no 2.5 here!
	double tau_relax = 1.0 * def_gamma * 0.5 * (mball / (r->N - num_perts))
			/ spring_mush.k; // Kelvin Voigt relaxation time
// factor of 0.5 is consistent with damping force law use of reduced mass
// this is equation 31 by Frouard+16
	printf("relaxation time %.3e\n", tau_relax);
	fprintf(fpr, "relaxation_time  %.3e\n", tau_relax);
	double viscosity = tau_relax * Emush_big / 2.5;
	printf("viscosity   %.3e\n", viscosity);
	fprintf(fpr, "viscosity   %.3e\n", viscosity);

// hard/soft interior!
	printf("Rshell = %.2f\n", Rshell);

	adjust_spring_props(r, kfac * ks, gammafac2 * spring_mush.gamma, 0.0, 0.0,
			Rshell);

	//compute moment of inertia matrix
	Matrix inertia_mat = mom_inertia(r, 0, r->N);
	// get its eigenvalues eig1>eig2>eig3
	double eigs[3];
	eigenvalues(inertia_mat, eigs);
	printf("sqrt I1/I3  %.3f should be an axis ratio \n",
			sqrt(eigs[2] / eigs[0]));
	printf("sqrt I1/I2  %.3f should be an axis ratio \n",
			sqrt(eigs[2] / eigs[1]));
	fprintf(fpr, "sqrt I1/I3  %.3f should be an axis ratio \n",
			sqrt(eigs[2] / eigs[0]));
	fprintf(fpr, "sqrt I1/I2  %.3f should be an axis ratio \n",
			sqrt(eigs[2] / eigs[1]));

	// compute angular momentum vector body with respect to its center
	// of mass position and velocity
	Vector L = measure_L(r, il, ih);

	// compute angle between angular momentum and principal axis and precession rate
	double theta = 0.0;
	double omega_3 = 0.0;
	double omega_prec = 0.0; // precession rate
	if (ratio1 == 1) { // oblate
		double hsquare = ratio2 * ratio2; // is (c/a)^2
		if (lattice_type == 3)
			hsquare = eigs[2] / eigs[0];
		theta = acos(L.getZ() / L.len()); // axis of symmetry is z axis
		omega_3 = L.len() / eigs[0]; // eig1 = I3 is largest moment of inertia, corresponding to smallest body axis
		omega_prec = (1.0 - hsquare) / (1.0 + hsquare);
		omega_prec *= omega_3 * cos(theta);
		printf("oblate\n");
	}
	if (ratio1 == ratio2) { // prolate
		double hsquare = 1.0 / (ratio2 * ratio2); // is a/c
		theta = acos(L.getX() / L.len()); // axis of symmetry is x axis
		omega_3 = L.len() / eigs[2]; // eig3 is smallest moment,  is I parallel, largest body axis
		omega_prec = (1.0 - hsquare) / (1.0 + hsquare);
		omega_prec *= omega_3 * cos(theta);
		printf("prolate\n");
	}
	printf("theta (deg)= %.6f (radians)= %.6f\n", theta * 180.0 / M_PI, theta);
	fprintf(fpr, "theta = %.6f radians, %.6f (degrees)\n", theta,
			theta * 180.0 / M_PI);
	printf("omega_prec= %.3e\n", omega_prec);
	fprintf(fpr, "omega_prec= %.3e\n", omega_prec);
	printf("omega_3= %.6f\n", omega_3);
	fprintf(fpr, "omega_3= %.6f\n", omega_3);

// do rotations afterwards!
// -------------------------------
	// I would like to orient the spinning body so that angular momentum is up
	// Euler angles: alpha,beta,gamma
	//   first angle alpha: about z in xy plane
	//   second angle beta: about x' axis in yz plane
	//   lastly gamma: rotate about z'' axis in xy plane
	// in our viewer up is y direction

	// compute angular momentum vector body with respect to its center
	// of mass position and velocity
	// double llx,lly,llz;

	L = measure_L(r, il, ih);
	// printf("llx lly llz = %.2f %.2f %.2f\n",llx,lly,llz);
	rotate_body(r, il, ih, 0.0, -atan2(L.getZ(), L.getY()), 0); // rotate by Euler angles, including velocities
	L = measure_L(r, il, ih);
	// printf("llx lly llz = %.3f %.3f %.3f\n",llx,lly,llz);
	rotate_body(r, il, ih, -atan2(L.getY(), L.getX()), 0, 0);
	L = measure_L(r, il, ih);
	// printf("llx lly llz = %.3f %.3f %.3f\n",llx,lly,llz);
	rotate_body(r, il, ih, 0, 0, M_PI / 2.0);
	L = measure_L(r, il, ih);
	printf("llx lly llz = %.3f %.3f %.3f\n", L.getX(), L.getY(), L.getZ());
	// ---------------- body should now be oriented with L up (in y direction)

	double barchi = 1.0 * fabs(omega_prec) * tau_relax; // initial value of barchi
	fprintf(fpr, "barchi  %.6e\n", barchi);
	printf("barchi %.6e\n", barchi);

	// ratio of numbers of particles to numbers of springs for resolved body
	double Nratio = (double) num_springs / (ih - il);
	printf("N=%d  NS=%d NS/N=%.1f\n", r->N, num_springs, Nratio);
	fprintf(fpr, "N=%d  NS=%d NS/N=%.1f\n", r->N, num_springs, Nratio);
	fclose(fpr);

	// r->particles[0].x += 1.1; // kick one particle
	reb_springs(r); // pass spring index list to display
	r->heartbeat = heartbeat;
	subtract_cov(r, il, ih); // center of velocity subtracted
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
	// static char pointmassfile[NPMAX*NSPACE];
	if (first == 0) {
		first = 1;
		string extendedfile = fileroot + "_ext.txt";
	}
	if (reb_output_check(r, 10.0 * r->dt)) {
		reb_output_timing(r, 0);
	}
	if (fabs(r->t - t_damp) < 0.9 * r->dt) { // set_gamma(r,def_gamma);
		set_gamma(def_gamma);
	}
	// damp initial bounce only
	// reset gamma only at t near t_damp

	// stuff to do every timestep
	center_sim(r, 0, r->N); // move reference frame to resolved body for display

	if (reb_output_check(r, print_interval)) {
		write_resolved_with_E(r, 0, r->N, extendedfile, 0.0); // orbital info and stuff
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
