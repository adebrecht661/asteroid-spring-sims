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
#include "libconfig.h++"
extern "C" {
#include "rebound.h"
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

extern vector<double> tot_power;

int num_springs = 0;  // number of springs
spring *springs;  // springs structure
#define NPMAX 10  // maximum number of point masses
double itaua[NPMAX], itaue[NPMAX];
int *surfp; // is allocated with marksurface routine!

double gamma_hot;
double gamma_fac; // for adjustment of gamma of all springs
double t_damp;    // end faster damping, relaxation
double t_print;   // for printouts
double t_heat;    // for heat  printout
string fileroot;
char froot[30];   // output files
int num_perts; // number of point mass perturbers
double powerfac; // fraction of heat to center of spring rather than nodes
double dote_radg;  // radiogenic heating rate
// energy per unit mass per unit time
void print_run_double(double quantity, string label, std::ofstream *file); // generic routine for printing information about run

int icentral = -1; // central mass location, is -1 if no central mass

// Global scales
double mass_scale, time_scale, length_scale, temp_scale, omega_scale, vel_scale,
		p_scale, L_scale, a_scale, F_scale, E_scale, dEdt_scale, P_scale;

void heartbeat(reb_simulation *const r);

void reb_springs(reb_simulation *const r);
void additional_forces(reb_simulation *r) {
	spring_forces(r); // spring forces
}

spring spring_mush_hot; // spring parameters for mush
spring spring_mush_cold; // spring parameters for mush
double Ttrans;  // Transition temperature , two temp model

int main(int argc, char *argv[]) {
	reb_simulation *const r = reb_create_simulation();
	// Setup constants
	r->integrator = reb_simulation::REB_INTEGRATOR_LEAPFROG;
	r->gravity = reb_simulation::REB_GRAVITY_BASIC;
	r->boundary = reb_simulation::REB_BOUNDARY_NONE;
	r->G = 1;
	r->additional_forces = additional_forces; // setup callback function for additional forces
	double mball = 1.0;          // total mass of ball
	double rball = 1.0;          // radius of a ball
	double tmax = 0.0;  // if 0 integrate forever

// things to set! ////////////////////// could be read in with parameter file
	double dt;             // timestep
	double b_distance;     // mininum interparticle spacing
	double omegaz;         // spin
	double mush_fac;       // mush_fac*b_distance is maximum spring length
	double surfdist;       // for identifying surface particles
	double ratio1, ratio2;  // body axis ratios
	double obliq_deg;
	double Tsurf = 0.0;   // surface temperature
	double Tinit = 0.0;   // initial interior node temperature
	dote_radg = 0.0;   // radiogenic heating rate
	int lattice_type = 0;  // for chosing body lattices
	double rad[NPMAX], mp[NPMAX]; // point mass radii and masses
	double aa[NPMAX], ee[NPMAX], ii[NPMAX];  // orbital elements
	double longnode[NPMAX], argperi[NPMAX], meananom[NPMAX]; // orbital elements
	int npointmass = 0;      // number of point masses
	num_perts = 0;               // number of point masses
	double k_heat_hot = 1.0;   // thermal transport coeffs
	double k_heat_cold = 0.0;
	gamma_hot = 0.01;    // damping parm
	double gamma_cold = 0.01;
	double ks_hot = 0.0;      // spring constant
	double ks_cold = 0.0;
	double cp_hot = 1.0;      // heat capacity
	double cp_cold = 1.0;
	Ttrans = 0.0;                // transition temperature two state model
	double R_shell = 0.0;
	double ba_shell = 1.0; // shell axis ratios
	double ca_shell = 1.0;
	double x_shell = 0.0; // shell offsets
	double y_shell = 0.0;
	double z_shell = 0.0;
	double ks_fac = 0.0; // lopsided parm for ks
	double fmfac = 0.1;   // density factor

	if (argc == 1) {
		fileroot = "t1";   // to make output files
		dt = 1e-3;    // Timestep
		lattice_type = 0;      // 0=rand 1=hcp
		b_distance = 0.15; // for creating random sphere, min separation between particles
		mush_fac = 2.3; // ratio of smallest spring distance to minimum interparticle dist
		omegaz = 0.2;     // initial spin
		// spring damping
		gamma_fac = 1.0; // initial factor for initial damping value for springs
		gamma_hot = 1.0;    // final damping coeff
		t_damp = 1.0;    // gamma to final values for all springs at this time
		ks_hot = 8e-2;   // spring constant

		ratio1 = 1.0; // shape of resolved body  y/x b/a
		ratio2 = 0.95; // z/x c/a
		t_print = 100000.0;  // printouts
		t_heat = 10000.0;    // heat printouts
		powerfac = 1.0;       // fraction of heat to center of springs
		obliq_deg = 0.0;        // obliquity in degrees
		surfdist = 0.1;         // for identifying surface particles
	} else {
		FILE *fpi; // read in a parameter file
		fpi = fopen(argv[1], "r");
		char line[300];
		fgets(line, 300, fpi);
		sscanf(line, "%s", froot);
		fileroot = froot;
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &dt);
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &tmax);
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &t_print);
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &t_heat);
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &powerfac);
		fgets(line, 300, fpi);
		sscanf(line, "%d", &lattice_type);
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &b_distance);
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &mush_fac);
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &ratio1);
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &ratio2);
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &omegaz);
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &obliq_deg);
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &gamma_fac);
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &t_damp);
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &surfdist);

		fgets(line, 300, fpi);
		sscanf(line, "%lf", &ks_hot);
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &ks_cold);
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &gamma_hot);
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &gamma_cold);
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &cp_hot);
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &cp_cold);
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &k_heat_hot);
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &k_heat_cold);
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &Ttrans);
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &Tinit);
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &dote_radg);
		fgets(line, 300, fpi);
		sscanf(line, "%lf", &R_shell);
		fgets(line, 300, fpi);
		sscanf(line, "%lf %lf", &ba_shell, &ca_shell);
		fgets(line, 300, fpi);
		sscanf(line, "%lf %lf %lf %lf %lf", &x_shell, &y_shell, &z_shell,
				&fmfac, &ks_fac);

		fgets(line, 300, fpi);
		sscanf(line, "%d", &npointmass);
		for (int ip = 0; ip < npointmass; ip++) {
			fgets(line, 300, fpi);
			sscanf(line, "%lf %lf %lf %lf", mp + ip, rad + ip, itaua + ip,
					itaue + ip);
			fgets(line, 300, fpi);
			sscanf(line, "%lf %lf %lf %lf %lf %lf", aa + ip, ee + ip, ii + ip,
					longnode + ip, argperi + ip, meananom + ip);
		}
		printf("parm file read in\n");

	}
	double obliquity = obliq_deg * M_PI / 180.0;
	if (powerfac > 1.0)
		powerfac = 1.0;
	if (powerfac < 0.0)
		powerfac = 1.0;

/// end of things to set /////////////////////////

	r->dt = dt;            // set integration timestep
	const double boxsize = 3.2 * rball;    // display window
	reb_configure_box(r, boxsize, 1, 1, 1);
	r->softening = b_distance / 100.0;	// Gravitational softening length
// viewer +x to right, +y to up, z back and forth along line of sight

	// properties of springs
	spring_mush_hot.gamma = gamma_hot; // initial damping coefficient
	spring_mush_cold.gamma = gamma_cold; // damping coefficient
	spring_mush_hot.k = ks_hot; // spring constant
	spring_mush_cold.k = ks_cold;
	spring_mush_hot.k_heat = k_heat_hot; // heat transport coefficient
	spring_mush_cold.k_heat = k_heat_cold;
	double mush_distance = b_distance * mush_fac;
	// distance for connecting and reconnecting springs

	string filename = fileroot + "_run.txt";
	std::ofstream outfile(filename, std::ios::out | std::ios::app);

// do you want volume to be the same? yes, adjusting here!!!!
	double volume_ratio = pow(rball, 3.0) * ratio1 * ratio2; // neglecting 4pi/3 factor
	double vol_radius = pow(volume_ratio, 1.0 / 3.0);

	rball /= vol_radius; // volume radius used to compute semi-major axis
// assuming that body semi-major axis is rball
	outfile << std::setprecision(3) << "a " << rball << "\n";
	outfile << std::setprecision(3) << "b " << rball * ratio1 << "\n";
	outfile << std::setprecision(3) << "c " << rball * ratio2 << "\n";
	volume_ratio = pow(rball, 3.0) * ratio1 * ratio2; // neglecting 4pi/3 factor
	outfile << std::setprecision(6) << "vol_ratio " << volume_ratio << "\n"; // with respect to 4pi/3
	// printf("vol_ratio %.6f\n",volume_ratio); // with respect to 4pi/3
	// so I can check that it is set to 1
	// working in units of volume equivalent sphere
	outfile << std::setprecision(6) << "fmfac " << fmfac << "\n"; //

	// create particle distribution
	if (lattice_type == 0) {
		rand_ellipsoid(r, b_distance, rball, rball * ratio1, rball * ratio2,
				mball);
		init_nodes(r, cp_hot, Tinit); // temperatures on nodes
		mark_surf_shrink_int_ellipsoid(r, surfdist, rball, rball * ratio1,
				rball * ratio2);
		init_surf_temp(Tsurf);
	}
	if (lattice_type == 1) {
		hcp_ellipsoid(r, b_distance, rball, rball * ratio1, rball * ratio2,
				mball);
		init_nodes(r, cp_hot, Tinit); // temperatures on nodes
		mark_surf_shrink_int_ellipsoid(r, surfdist, rball, rball * ratio1,
				rball * ratio2);
		init_surf_temp(Tsurf);
	}
	if (lattice_type > 1) {
		exit(0);
	}

	int il = 0;
	int ih = r->N;
	subtract_com(r, il, ih);   // move reference frame to resolved body
	subtract_cov(r, il, ih);  // subtract center of velocity

	// make springs, all pairs connected within interparticle distance mush_distance
	connect_springs_dist(r, mush_distance, 0, r->N, spring_mush_hot);

	tot_power.resize(num_springs, 0.0);

	// assume minor semi is rball*ratio2
	double ddr = rball * ratio2 - 0.5 * mush_distance;
	ddr = 0.5;
	double Emush = Young_mesh(r, il, ih, 0.0, ddr);
	double Emush_big = Young_full_mesh();
	print_run_double(mush_distance, "max spring length", &outfile);
	print_run_double(ddr, "ddr", &outfile);
	print_run_double(Emush, "Young's modulus hot", &outfile);
	print_run_double(Emush_big, "Young's modulus big hot", &outfile);
	print_run_double(Emush * ks_cold / ks_hot, "Young's modulus cold",
			&outfile);
	print_run_double(Emush_big * ks_cold / ks_hot, "Young's modulus big cold",
			&outfile);

	double K_T = therm_cond_mesh(r, il, ih, 0.0, ddr); // depends on k_heat
	print_run_double(K_T, "Thermal conductivity hot", &outfile);
	print_run_double(K_T * k_heat_cold / k_heat_hot,
			"Thermal conductivity cold", &outfile);
	print_run_double(mean_spring_length(), "Mean spring length", &outfile);

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
			if (icentral == -1) {
				double na = om * aa[ip];
				double adot = 3.0 * mp[ip] * na / pow(aa[ip], 5.0); // should approximately be adot
				outfile << std::setprecision(3) << "adot " << adot << "\n";
			}
			outfile << std::setprecision(3) << "resbody mm= " << om
					<< std::setprecision(2) << " period= " << 2.0 * M_PI / om
					<< "\n";
			printf("resbody mm=%.3f period=%.2f\n", om, 2.0 * M_PI / om);
			// note central body mass (only changes on first run - can probably do this better
			icentral = ih;
		}
		num_perts = npointmass;
	}

	// factor of 0.5 is due to reduced mass being used in calculation
	double tau_relax_hot = 1.0 * gamma_hot * 0.5 * (mball / (r->N - num_perts))
			/ spring_mush_hot.k; // Kelvin Voigt relaxation time
	double tau_relax_cold = 1.0 * gamma_cold * 0.5
			* (mball / (r->N - num_perts)) / spring_mush_cold.k; // Kelvin Voigt relaxation time
	print_run_double(tau_relax_hot, "relaxation time hot", &outfile);
	print_run_double(tau_relax_cold, "relaxation time cold", &outfile);

	double barchi_hot = 2.0 * fabs(om) * tau_relax_hot; // initial value of barchi
	double barchi_cold = 2.0 * fabs(om) * tau_relax_cold; // initial value of barchi
	print_run_double(barchi_hot, "barchi hot", &outfile);
	print_run_double(barchi_cold, "barchi cold", &outfile);

	double speriod = fabs(2.0 * M_PI / omegaz);
	print_run_double(speriod, "spin period", &outfile);
	print_run_double(omegaz, "omegaz", &outfile);

	double Nratio = (double) num_springs / (double) r->N;
	printf("N=%d  NS=%d NS/N=%.1f\n", r->N, num_springs, Nratio);
	outfile << "N= " << r->N << " NS=" << num_springs << " NS/N=" << Nratio
			<< "\n";
	outfile.close();

	// now spin it
	spin_body(r, il, ih, Vector( { 0.0, 0.0, omegaz })); // you can change one of these to tilt!
	// if (obliquity != 0.0)
	//   rotate_body(r, il, ih, 0.0, obliquity, 0.0); // tilt by obliquity in radians
	// note: if orbit is inclined this is not actually the obliquity

	Vector L = measure_L_origin(r, 0, r->N); // total angular momentum
	printf("llx lly llz = %.6e %.6e %.6e\n", L.getX(), L.getY(), L.getZ());

	rotate_origin(r, 0, r->N, 0, -atan2(L.getY(), L.getZ()), 0); // rotate by Euler angles, including velocities
	L = measure_L_origin(r, 0, r->N); // total angular momentum
	printf("llx lly llz = %.6e %.6e %.6e\n", L.getX(), L.getY(), L.getZ());

	rotate_origin(r, 0, r->N, 0, 0, M_PI / 2); // rotate by Euler angles, including velocities
	L = measure_L_origin(r, 0, r->N); // total angular momentum
	printf("llx lly llz = %.6e %.6e %.6e\n", L.getX(), L.getY(), L.getZ());

	rotate_origin(r, 0, r->N, 0, -atan2(L.getY(), L.getZ()), 0); // rotate by Euler angles, including velocities
	L = measure_L_origin(r, 0, r->N); // total angular momentum
	printf("llx lly llz = %.6e %.6e %.6e\n", L.getX(), L.getY(), L.getZ());
// right now seems like orbit starts with big guy in negative y

	if (obliquity != 0.0)
		rotate_body(r, il, ih, 0.0, obliquity, 0.0); // tilt by obliquity in radians

	if (R_shell > 0) {
		adjust_spring_props_ellipsoid(r, spring_mush_cold.k, spring_mush_cold.gamma,
				spring_mush_cold.k_heat, R_shell, ba_shell * R_shell,
				ca_shell * R_shell, Vector({x_shell, y_shell, z_shell}), false); // x0,y0,z0, not-inside
		// adjust springs outside an ellipsoid with semi-axes: a=R, b=ba*R, c=ca*R
		// and with center x_shell,y_shell,z_shell

		if (fabs(fmfac) > 1e-3)
			adjust_mass_ellipsoid(r, fmfac + 1.0, R_shell,
					ba_shell * R_shell, ca_shell * R_shell, Vector({x_shell,
					y_shell + 0.0, z_shell}), false); // x0,y0,z0, not-inside
		// adjust density by mfac!
		// adjust_nodes_cp(r, num_perts, nodevec,Tinit/2,cp_cold,1); // adjust below transition temp

// lopsided softening of springs in shell
		if (fabs(ks_fac) > 1e-3) {
			// Approximately the period of the function
			double freq = 1.0;
			double phi0 = 0.0; // which side is lopsided in spring strengths
			adjust_spring_props_ellipsoid_phase(r, ks_fac,
					0.0, 0.0, freq, phi0, R_shell, ba_shell * R_shell,
					ca_shell * R_shell, Vector({x_shell, y_shell, z_shell}), 0); // x0,y0,z0, not-inside
		}
	}

	reb_springs(r); // pass spring index list to display
	set_gamma(gamma_hot * gamma_fac);

	r->heartbeat = heartbeat;

	subtract_com(r, 0, r->N - num_perts); // move reference frame, position only
	if (tmax == 0.0)
		reb_integrate(r, INFINITY);
	else
		reb_integrate(r, tmax);
}

#define NSPACE 50
void heartbeat(reb_simulation *const r) {
	static int first = 0;
	static char extendedfile[50];
	static char twopfile[50];
	static char pointmassfile[NPMAX * NSPACE];
	// static int j0 = 0;
	if (first == 0) {
		first = 1;
		sprintf(extendedfile, "%s_ext.txt", froot);
		sprintf(twopfile, "%s_2p.txt", froot);
		for (int i = 0; i < num_perts; i++) {
			sprintf(pointmassfile + i * NSPACE, "%s_pm%d.txt", froot, i);
		}
		subtract_com(r, 0, r->N - num_perts); // move reference frame, position only
		// j0 = nearest_to_shape(r,0,r->N-num_perts,0.0,0.0,0.0); // find particle nearest origin
		// printf("j0=%d\n",j0);
	}

	if (reb_output_check(r, 10.0 * r->dt)) {
		reb_output_timing(r, 0);
	}
	if (abs(r->t - t_damp) < 0.9 * r->dt) {
		set_gamma(gamma_hot);
	}
	// damp initial bounce only
	// reset gamma only at t near t_damp

	// stuff to do every timestep
	subtract_com(r, 0, r->N - num_perts); // move reference frame, position only
	// don't store heat until after dampdown!
	if (r->t > t_damp) {

		// Don't think this is necessary?????? (except for debug)
		rec_power(r); // store heat dE/dt accumulated each timestep
		// in heatvec
		// the heat vector is zeroed and normalized only when printed out

		heat_int_nodes_tidal(r, r->dt);
		// stores tidal heat every timestep
		// applies 0.5*dEdt from each spring to the node (does not use heatvec)
		// tidal heating only right now
		// nodevec[j0].temp = 1.0; // fix temperature of node near origin

		transport_heat(r, r->dt); // adjust temps over network

	}
	if (reb_output_check(r, t_heat)) { // heat files
		int ndt = (int) (t_heat / r->dt);  // how long tidal heat is stored up
		string hfile = heat_filename(r, froot, t_heat);
		write_heat(r, hfile, ndt, powerfac); // heat info printed out!

		string nfile = node_filename(r, froot, t_heat);
		write_nodes(r, nfile); // temperature info printed out!
		if (r->t > t_damp) {
			// adjust_spring_temp_ts(r, nodevec,
			//   spring_mush_hot, spring_mush_cold, Ttrans);
		}
	}

	if (reb_output_check(r, t_print)) {
		write_resolved_no_E(r, 0, r->N - num_perts, extendedfile); // orbital info and stuff
		write_resolved_2nodes(r, 0, r->N - num_perts, twopfile);
		// exit(0);
		if (num_perts > 0)
			for (int i = 0; i < num_perts; i++) {
				int ip = icentral + i;
				write_pt_mass(r, ip, i, pointmassfile + i * NSPACE);
			}
	}

}

// make a spring index list
void reb_springs(reb_simulation *const r) {
	r->NS = num_springs;
	r->springs_i = (int*) malloc(num_springs * sizeof(int));
	r->springs_j = (int*) malloc(num_springs * sizeof(int));
	for (int i = 0; i < num_springs; i++) {
		r->springs_i[i] = springs[i].particle_1;
		r->springs_j[i] = springs[i].particle_2;
	}
}

// Print doubles to file and standard out
void print_run_double(double quantity, string label, std::ofstream *outfile) {
	// Set precision base on size of quantity
	if ((abs(log10(quantity)) > 4)) {
		std::cout << label << std::setprecision(4) << " " << quantity << "\n";
		*outfile << label << std::setprecision(4) << " " << quantity << "\n";
	} else {
		std::cout << label << std::setprecision(3) << " " << quantity << "\n";
		*outfile << label << std::setprecision(3) << " " << quantity << "\n";
	}
}
