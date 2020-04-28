#ifdef __cplusplus
# 	ifdef __GNUC__
#		define restrict __restrict__
#	else
#		define restrict
#	endif
#endif

/*
 * output.cpp
 *
 * Output routines
 *
 *  Created on: Apr 15, 2020
 *      Author: alex
 */

#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
extern "C" {
#include "rebound.h"
}
#include "matrix_math.h"
#include "physics.h"
#include "springs.h"
#include "stress.h"
#include "heat.h"
#include "shapes.h"
#include "orb.h"
#include "input_spring.h" // For padding function - nowhere good for misc functions
#include "output_spring.h"

using std::string;
using std::vector;

extern vector<spring> springs;
extern int num_springs;
extern int num_perts;
extern vector<stress_tensor> stresses;
extern vector<node> nodes;
extern vector<double> tot_power;

/********************/
/* Output functions */
/********************/

// Write out positions and velocities of surface particles
// Can write at multiple times
void write_surf_part(reb_simulation *const n_body_sim, int i_low, int i_high,
		string filename) {
	// Get particle info
	reb_particle *particles = n_body_sim->particles;

	// Open output file to append (file is created if it doesn't exist)
	std::ofstream outfile(filename, std::ios::out | std::ios::app);

	// Get time of first print for constant offset for future prints
	static bool first = true;
	static double t_off = 0.0;
	if (first) {
		first = false;
		t_off = n_body_sim->t;
		std::cout << std::setprecision(6) << "print_surf: t_off = " << t_off
				<< "\n";

		// Print header only once
		outfile << "# t\tx\ty\tz\tvx\tvy\tvz\tax\tay\taz\tr\tm\n";
	}

	// Write actual current time
	outfile << std::setprecision(8) << "# t = " << n_body_sim->t << " #\n";

	// Print tab-separated particle info at time t - t_off
	for (int i = 0; i < n_body_sim->N; i++) {
		if (nodes[i].is_surf)
			outfile << std::setprecision(6) << n_body_sim->t - t_off << "\t"
					<< particles[i].x << "\t" << particles[i].y << "\t"
					<< particles[i].z << "\t" << particles[i].vx << "\t"
					<< particles[i].vy << "\t" << particles[i].vz << "\t"
					<< particles[i].ax << "\t" << particles[i].ay << "\t"
					<< particles[i].az << "\t" << particles[i].r << "\t"
					<< particles[i].m << "\n";
	}

	// Close file and report
	outfile.close();
	std::cout << "print_surf: Printed " << n_body_sim->N
			<< " particles to file " << filename << std::endl;
}

// Write out information for an extended body with indices [i_low, i_high)
void write_resolved_with_E(reb_simulation *const n_body_sim, int i_low,
		int i_high, string filename, double dEdt_ave) {

	// Open output file to append (created if it doesn't exist)
	std::ofstream outfile(filename, std::ios::out | std::ios::app);

	// Print header once
	static bool first = true;
	if (first) {
		first = false;
		outfile
				<< "# t\tx\ty\tz\tvx\tvy\tvz\tomx\tomy\tomz\tllx\tlly\tllz\tIxx\tIyy\tIzz\tIxy\tIyz\tIxz\tKErot\tPEspr\tPEgrav\tEtot\tdEdtnow\tdEdtave\n";
	}

	// Time
	outfile << std::setprecision(3) << n_body_sim->t << "\t";

	// Center of mass location and velocity
	Vector CoM = compute_com(n_body_sim, i_low, i_high);
	Vector CoV = compute_cov(n_body_sim, i_low, i_high);
	outfile << std::setprecision(5) << CoM.getX() << "\t" << CoM.getY() << "\t"
			<< CoM.getZ() << "\t";
	outfile << std::setprecision(5) << CoV.getX() << "\t" << CoV.getY() << "\t"
			<< CoV.getZ() << "\t";

	// Spin angular momentum
	// What's the difference between omega and L???????
	double eigs[3];
	Vector omega = body_spin(n_body_sim, i_low, i_high, eigs);
	Vector L = measure_L(n_body_sim, i_low, i_high);
	outfile << std::setprecision(4) << omega.getX() << "\t" << omega.getY()
			<< "\t" << omega.getZ() << "\t";
	outfile << std::setprecision(4) << L.getX() << "\t" << L.getY() << "\t"
			<< L.getZ() << "\t";

	// Moment of inertia
	Matrix inertia_mat = mom_inertia(n_body_sim, i_low, i_high);
	outfile << std::setprecision(5) << inertia_mat.getXX() << "\t"
			<< inertia_mat.getYY() << "\t" << inertia_mat.getZZ() << "\t"
			<< inertia_mat.getXY() << "\t" << inertia_mat.getYZ() << "\t"
			<< inertia_mat.getXZ() << "\t";

	// Non-translational kinetic energy
	double KErot = compute_rot_kin(n_body_sim, i_low, i_high);
	outfile << std::setprecision(5) << KErot << "\t";

	// Spring potential energy
	double SPE = spring_potential_energy(n_body_sim);
	outfile << std::setprecision(5) << SPE << "\t";

	// Gravitational potential energy
	double GPE = grav_potential_energy(n_body_sim, i_low, i_high);
	outfile << std::setprecision(10) << GPE << "\t";

	// Total internal energy
	outfile << std::setprecision(10) << (KErot + SPE + GPE) << "\t";

	// Total power dissipated at this time by spring damping
	outfile << std::setprecision(10) << dEdt_total(n_body_sim);

	// Average power dissipation by spring damping
	outfile << std::setprecision(10) << dEdt_ave << "\n";
	outfile.close();
}

// Write out information for an extended body with indices [i_low, i_high)
// No energy terms
void write_resolved_no_E(reb_simulation *const n_body_sim, int i_low,
		int i_high, string filename) {

	// Open output file to append (created if it doesn't exist)
	std::ofstream outfile(filename, std::ios::out | std::ios::app);

	// Print header once
	static bool first = true;
	if (first) {
		first = false;
		outfile
				<< "# t\tx\ty\tz\tvx\tvy\tvz\tomx\tomy\tomz\tllx\tlly\tllz\tIxx\tIyy\tIzz\tIxy\tIyz\tIxz\n";
	}

	// Time
	outfile << std::setprecision(3) << n_body_sim->t << "\t";

	// Center of mass location and velocity
	Vector CoM = compute_com(n_body_sim, i_low, i_high);
	Vector CoV = compute_cov(n_body_sim, i_low, i_high);
	outfile << std::setprecision(5) << CoM.getX() << "\t" << CoM.getY() << "\t"
			<< CoM.getZ() << "\t";
	outfile << std::setprecision(5) << CoV.getX() << "\t" << CoV.getY() << "\t"
			<< CoV.getZ() << "\t";

	// Spin angular momentum
	// What's the difference between omega and L???????
	double eigs[3];
	Vector omega = body_spin(n_body_sim, i_low, i_high, eigs);
	Vector L = measure_L(n_body_sim, i_low, i_high);
	outfile << std::setprecision(4) << omega.getX() << "\t" << omega.getY()
			<< "\t" << omega.getZ() << "\t";
	outfile << std::setprecision(4) << L.getX() << "\t" << L.getY() << "\t"
			<< L.getZ() << "\t";

	// Moment of inertia
	Matrix inertia_mat = mom_inertia(n_body_sim, i_low, i_high);
	outfile << std::setprecision(5) << inertia_mat.getXX() << "\t"
			<< inertia_mat.getYY() << "\t" << inertia_mat.getZZ() << "\t"
			<< inertia_mat.getXY() << "\t" << inertia_mat.getYZ() << "\t"
			<< inertia_mat.getXZ() << "\n";
	outfile.close();
}

// Write out information for two specific nodes in an extended body with indices [i_low, i_high)
// The first node has largest x value at initial time
// The second node has largest z value at initial time
// Same nodes are printed each time after they are chosen
void write_resolved_2nodes(reb_simulation *const n_body_sim, int i_low,
		int i_high, string filename) {
	// Get particle info
	reb_particle *particles = n_body_sim->particles;

	// Indices of particles with largest z and x values
	static int i_x = 0;
	static int i_z = 0;

	// Open output file to append (created if it doesn't exist)
	std::ofstream outfile(filename, std::ios::out | std::ios::app);

	// Write header once
	static bool first = true;
	if (first) {
		first = false;
		outfile << "# xnode \t\t\t znode #\n";
		outfile << "# t\tx\ty\tz\tvx\tvy\tvz\t\tx\ty\tz\tvx\tvy\tvz\n";

		// Get particles with largest x and z values (points should be well outside the system)
		i_x = nearest_to_point(n_body_sim, i_low, i_high, Vector( { 10.0, 0.0,
				0.0 }));
		i_z = nearest_to_point(n_body_sim, i_low, i_high, Vector( { 0.0, 0.0,
				10.0 }));
	}

	// Time
	outfile << std::setprecision(3) << n_body_sim->t << "\t";

	// x particle location and velocity
	outfile << std::setprecision(5) << particles[i_x].x << "\t"
			<< particles[i_x].y << "\t" << particles[i_x].z << "\t"
			<< particles[i_x].vx << "\t" << particles[i_x].vy << "\t"
			<< particles[i_x].vz << "\t\t";

	// z particle location and velocity
	outfile << std::setprecision(5) << particles[i_z].x << "\t"
			<< particles[i_z].y << "\t" << particles[i_z].z << "\t"
			<< particles[i_z].vx << "\t" << particles[i_z].vy << "\t"
			<< particles[i_z].vz << "\n";
	outfile.close();
}

// Write out information for a point mass with particle index i_pert
// pert_num is count of perturbing particle
void write_pt_mass(reb_simulation *const n_body_sim, int i_pert, int pert_num,
		string filename) {
	// Get particle info
	reb_particle *particles = n_body_sim->particles;

	// Record whether it's the first run for each particle
	const int max_perts = 100;
	static bool first[max_perts];

	// Initialize boolean array
	static bool first_run = true;
	if (first_run) {
		for (int i = 0; i < num_perts; i++) {
			first[i] = true;
		}
	}

	// Open output file to append (created if it doesn't exist)
	std::ofstream outfile(filename, std::ios::out | std::ios::app);

	// Write header once per particle
	if (first[pert_num]) {
		first[pert_num] = false;
		outfile << "# t\tx\ty\tz\tvx\tvy\tvz\tm\n";
	}

	// Time
	outfile << std::setprecision(3) << n_body_sim->t << "\t";

	// Particle info
	outfile << std::setprecision(5) << particles[i_pert].x << "\t"
			<< particles[i_pert].y << "\t" << particles[i_pert].z << "\t"
			<< particles[i_pert].vx << "\t" << particles[i_pert].vy << "\t"
			<< particles[i_pert].vz << "\t";
	outfile << std::setprecision(3) << particles[i_pert].m << "\n";
	outfile.close();
}

// Write information about resolved body orbiting one or more point particles
void write_resolved_orb(reb_simulation *const n_body_sim, string filename) {

	// Open output file to append (created if it doesn't exist)
	std::ofstream outfile(filename, std::ios::out | std::ios::app);

	// Write header once
	static bool first = true;
	if (first) {
		first = false;
		outfile
				<< "# t\ta\tn\te\ti\tomx\tomy\tomz\tA\tB\tC\tE\tlx\tly\tlz\tang\tlox\tloy\tloz\tIxx\tIyy\tIzz\tIxy\tIyz\tIxz\n";
	}

	// Get indices of resolved body
	int i_low = 0;
	int i_high = n_body_sim->N - num_perts;

	// Time
	outfile << std::setprecision(3) << n_body_sim->t << "\t";

	// Orbital properties of resolved body around perturbers
	double a, mean_motion, ecc, incl, Lorb;
	compute_orb(n_body_sim, i_low, i_high, &a, &mean_motion, &ecc, &incl,
			&Lorb);
	outfile << std::setprecision(6) << a << "\t";
	outfile << std::setprecision(4) << mean_motion << "\t" << ecc << "\t"
			<< incl << "\t";

	// Spin
	// Why both omega and L??????
	double eigs[3];
	Vector omega = body_spin(n_body_sim, i_low, i_high, eigs);
	outfile << std::setprecision(4) << omega.getX() << "\t" << omega.getY()
			<< "\t" << omega.getZ() << "\t";
	outfile << std::setprecision(3) << eigs[0] << "\t" << eigs[1] << "\t"
			<< eigs[2] << "\t";

	// Young's modulus
	double r_min = 0.0, r_max = 0.5;
	outfile << std::setprecision(3)
			<< Young_mesh(n_body_sim, i_low, i_high, r_min, r_max) << "\t";

	// Spin
	// Why both omega and L??????
	Vector L = measure_L(n_body_sim, i_low, i_high);
	outfile << std::setprecision(5) << L.getX() << "\t" << L.getY() << "\t"
			<< L.getZ() << "\t";

	// Get orbital angular momentum
	Vector L_o = compute_Lorb(n_body_sim, i_low, i_high);

	// Angle between spin and orbit
	outfile << std::setprecision(5) << acos(dot(L, L_o) / (L.len() * L_o.len()))
			<< "\t";

	// Output orbital angular momentum
	outfile << std::setprecision(5) << L_o.getX() << "\t" << L_o.getY() << "\t"
			<< L_o.getZ() << "\t";

	// Moment of inertia
	Matrix inertia_mat = mom_inertia(n_body_sim, i_low, i_high);
	outfile << std::setprecision(5) << inertia_mat.getXX() << "\t"
			<< inertia_mat.getYY() << "\t" << inertia_mat.getZZ() << "\t"
			<< inertia_mat.getXY() << "\t" << inertia_mat.getYZ() << "\t"
			<< inertia_mat.getXZ() << "\n";
	outfile.close();
}

// Write out resolved body and point particle info
void write_resolved_bin(reb_simulation *const n_body_sim, string filename) {
	// Get particle info
	reb_particle *particles = n_body_sim->particles;

	// Get particle with largest radius
	static int part_largest_r = -1;

	// Open output file to append (created if it doesn't exist)
	std::ofstream outfile(filename, std::ios::out | std::ios::app);

	// Write header once
	static bool first = true;
	if (first) {
		first = false;
		outfile
				<< "# resolved body to perts\t\t\tlargest radius to resolved\t\t\t1st pert to 2nd pert\t\t\tresolved body\n";
		outfile
				<< "# t\tx\ty\tz\tvx\tvy\tvz\t\tx\ty\tz\tvx\tvy\tvz\t\tx\ty\tz\tvx\tvy\tvz\t\tlx\tly\tlz\tIxx\tIyy\tIzz\tIxy\tIyz\tIxz\n";

		// Find particle with largest radius
		double r_max = 0.0;
		for (int i = 0; i < n_body_sim->N - num_perts; i++) {
			Vector x = { particles[i].x, particles[i].y, particles[i].z };
			double r = x.len();
			if (r > r_max) {
				part_largest_r = i;
				r_max = r;
			}
		}
	}

	// Indices of resolved body
	int i_low = 0;
	int i_high = n_body_sim->N - num_perts;

	// Indices of perturbing masses
	int i_pert_low = n_body_sim->N - num_perts;
	int i_pert_high = n_body_sim->N;

	// Time
	outfile << std::setprecision(3) << n_body_sim->t << "\t";

	// Center of mass of resolved body location and velocity
	Vector CoM = compute_com(n_body_sim, i_low, i_high);
	Vector CoV = compute_cov(n_body_sim, i_low, i_high);

	// Center of mass of perturbing particles location and velocity
	Vector CoM_pert = compute_com(n_body_sim, i_pert_low, i_pert_high);
	Vector CoV_pert = compute_cov(n_body_sim, i_pert_low, i_pert_high);

	// Relative position and velocity of resolved body to center of perturbers
	Vector xr = CoM - CoM_pert;
	Vector vr = CoV - CoV_pert;
	outfile << std::setprecision(5) << xr.getX() << "\t" << xr.getY() << "\t"
			<< xr.getZ() << "\t" << vr.getX() << "\t" << vr.getY() << "\t"
			<< vr.getZ();

	// Relative position and velocity of particle with largest radius to resolved body
	Vector x_largest = { particles[part_largest_r].x,
			particles[part_largest_r].y, particles[part_largest_r].z };
	Vector v_largest = { particles[part_largest_r].vx,
			particles[part_largest_r].vy, particles[part_largest_r].vz };
	Vector xs = x_largest - CoM;
	Vector vs = v_largest - CoV;
	outfile << std::setprecision(5) << xs.getX() << "\t" << xs.getY() << "\t"
			<< xs.getZ() << "\t" << vs.getX() << "\t" << vs.getY() << "\t"
			<< vs.getZ() << "\t";

	// If there's a binary set of perturbers, print out secondary's position and velocity relative to primary
	if (num_perts == 2) {
		Vector x_prim = { particles[i_pert_low + 1].x,
				particles[i_pert_low + 1].y, particles[i_pert_low + 1].z };
		Vector v_prim = { particles[i_pert_low + 1].vx,
				particles[i_pert_low + 1].vy, particles[i_pert_low + 1].vz };
		Vector x_sec = { particles[i_pert_low].x, particles[i_pert_low].y,
				particles[i_pert_low].z };
		Vector v_sec = { particles[i_pert_low].vx, particles[i_pert_low].vy,
				particles[i_pert_low].vz };
		Vector x_rel = x_sec - x_prim;
		Vector v_rel = v_sec - v_prim;

		outfile << std::setprecision(5) << x_rel.getX() << "\t" << x_rel.getY()
				<< "\t" << x_rel.getZ() << "\t" << v_rel.getX() << "\t"
				<< v_rel.getY() << "\t" << v_rel.getZ() << "\t";
	} else {
		outfile << "0\t0\t0\t0\t0\t0\t";
	}

	// Spin angular momentum
	Vector L = measure_L(n_body_sim, i_low, i_high);
	outfile << std::setprecision(5) << L.getX() << "\t" << L.getY() << "\t"
			<< L.getZ() << "\t";

	// Moment of inertia
	Matrix inertia_mat = mom_inertia(n_body_sim, i_low, i_high);
	outfile << std::setprecision(5) << inertia_mat.getXX() << "\t"
			<< inertia_mat.getYY() << "\t" << inertia_mat.getZZ() << "\t"
			<< inertia_mat.getXY() << "\t" << inertia_mat.getYZ() << "\t"
			<< inertia_mat.getXZ() << "\n";
	outfile.close();
}

// Write out springs to a file
void write_springs(reb_simulation *const n_body_sim, string fileroot,
		int index) {
	// Set filename
	string filename = fileroot + "_" + zero_pad_int(6, index) + "_springs.txt";

	// Open output file to append (created if it doesn't exist)
	std::ofstream outfile(filename, std::ios::out | std::ios::app);

	// Write out each spring
	for (int i = 0; i < num_springs; i++) {
		outfile << springs[i].particle_1 << "\t" << springs[i].particle_2
				<< "\t" << std::setprecision(6) << springs[i].k << "\t"
				<< springs[i].rs0 << "\t" << springs[i].gamma << "\t"
				<< springs[i].k_heat << "\t" << strain(n_body_sim, springs[i])
				<< "\t" << spring_r(n_body_sim, springs[i]).len() << "\n";
	}
	outfile.close();
	std::cout << "write_springs: Printed " << num_springs << " springs to "
			<< filename << std::endl;
}

// Write out particles to file
void write_particles(reb_simulation *const n_body_sim, string fileroot,
		int index) {
	// Get particle info
	reb_particle *particles = n_body_sim->particles;

	// Set filename
	string filename = fileroot + "_" + zero_pad_int(6, index)
			+ "_particles.txt";

	// Open output file to append (created if it doesn't exist)
	std::ofstream outfile(filename, std::ios::out | std::ios::app);

	// Write out each particle
	for (int i = 0; i < n_body_sim->N; i++) {
		outfile << std::setprecision(6) << particles[i].x << "\t"
				<< particles[i].y << "\t" << particles[i].z << "\t"
				<< particles[i].vx << "\t" << particles[i].vy << "\t"
				<< particles[i].vz << "\t" << particles[i].r << "\t"
				<< particles[i].m << "\n";
	}
	outfile.close();
	std::cout << "write_particles: Printed " << n_body_sim->N
			<< " particles to " << filename;
}

// Write node temperature file
void write_nodes(reb_simulation *const n_body_sim, string filename) {
	// Get particle info
	reb_particle *particles = n_body_sim->particles;

	// Open output file to append (created if it doesn't exist)
	std::ofstream outfile(filename, std::ios::out | std::ios::app);

	// Print header first time only
	static bool first = true;
	if (first) {
		outfile << "# i\tx\ty\tz\tvx\tvy\tvz\tT\tcv\tsurf\tm\txrot\tyrot\n";
		first = false;
	}

	// Time
	outfile << "# t = " << std::setprecision(3) << n_body_sim->t << "\n";

	// Get indices for resolved body
	int i_low = 0;
	int i_high = n_body_sim->N - num_perts;

	// Get center of mass of resolved body
	Vector CoM = compute_com(n_body_sim, i_low, i_high);

	// Get position of most recent perturber
	int i_pert_1 = n_body_sim->N - 1;
	Vector x_1 = { particles[i_pert_1].x, particles[i_pert_1].y,
			particles[i_pert_1].z };

	// Get relative direction of perturber
	double theta = atan2(x_1.getY() - CoM.getY(), x_1.getX() - CoM.getX());
	Matrix rot_mat = getRotMatZ(theta);
	outfile << "# x1\ty1\tz1\txb\tyb\tzb\ttheta\n";
	outfile << std::setprecision(3) << x_1.getX() << "\t" << x_1.getY() << "\t"
			<< x_1.getZ() << "\t" << CoM.getX() << "\t" << CoM.getY() << "\t"
			<< CoM.getZ() << "\t" << theta << "\n";

	for (int i = i_low; i < i_high; i++) {
		Vector x = { particles[i].x, particles[i].y, particles[i].z };
		Vector x_rot = rot_mat * x;
		outfile << std::setprecision(4) << i << "\t";
		outfile << std::setprecision(10) << particles[i].x << "\t"
				<< particles[i].y << "\t" << particles[i].z << "\t"
				<< particles[i].vx << "\t" << particles[i].vy << "\t"
				<< particles[i].vz << "\t";
		outfile << std::setprecision(10) << nodes[i].temp << "\t" << nodes[i].cv
				<< "\t" << nodes[i].is_surf << "\t";
		outfile << std::setprecision(10) << x_rot.getX() << "\t" << x_rot.getY()
				<< "\n";
	}
	outfile.close();
}

// Write out stress file
void write_stresses(reb_simulation *const n_body_sim, string filename) {
	// Get particle info
	reb_particle *particles = n_body_sim->particles;

	// Open output file to append (created if it doesn't exist)
	std::ofstream outfile(filename, std::ios::out | std::ios::app);

	// Print header first time only
	static bool first = true;
	if (first) {
		outfile << "# i\tm\tx\ty\tz\teig1\teig2\teig3\tfailing\n";
		first = false;
	}

	// Time
	outfile << "# t = " << std::setprecision(2) << n_body_sim->t << "\n";

	// Get resolved body indices
	int i_low = 0;
	int i_high = n_body_sim->N - num_perts;

	// Get center of mass of resolved body
	Vector CoM = compute_com(n_body_sim, i_low, i_high);

	// Print stress of each particle
	for (int i = i_low; i < i_high; i++) {
		Vector x = { particles[i].x, particles[i].y, particles[i].z };
		Vector dx = x - CoM;
		outfile << i << "\t" << std::setprecision(4) << particles[i].m << "\t"
				<< std::setprecision(3) << dx.getX() << "\t" << dx.getY()
				<< "\t" << dx.getZ() << "\t";
		outfile << std::setprecision(3) << stresses[i].eigs[0] << "\t"
				<< stresses[i].eigs[1] << "\t" << stresses[i].eigs[2] << "\t";
		outfile << stresses[i].failing << "\n";
	}
	outfile.close();
}

// Write out heat file
// Takes num_timesteps that we've recorded power at
// powerfac is how much of heat to distribute at center of spring rather than at nodes
// also outputs rotated xy, after rotates assuming xy plane for orbit and toward pertuber
void write_heat(reb_simulation *const n_body_sim, string filename,
		int num_timesteps, double power_fac) {
	reb_particle *particles = n_body_sim->particles;

	// Get average power over our timesteps
	normalize_tot_power(num_timesteps);

	// Open output file to append (created if it doesn't exist)
	std::ofstream outfile(filename, std::ios::out | std::ios::app);

	// Print header first time only
	static bool first = true;
	if (first) {
		outfile << "# x\ty\tz\tdedt\txrot\tyrot\n";
		first = false;
	}

	// Time
	outfile << std::setprecision(2) << "# t = " << n_body_sim->t;

	// Get resolved body indices
	int i_low = 0;
	int i_high = n_body_sim->N - num_perts;

	// Get center of mass
	Vector CoM = compute_com(n_body_sim, i_low, i_high);

	// Index for primary perturber
	int i_prim = n_body_sim->N - 1;

	// Location of primary perturber
	Vector x_prim = { particles[i_prim].x, particles[i_prim].y,
			particles[i_prim].z };
	double theta = atan2(x_prim.getY() - CoM.getY(),
			x_prim.getX() - CoM.getX());

	// Print primary perturber info
	outfile << std::setprecision(3) << "# primary info: " << x_prim.getX()
			<< "\t" << x_prim.getY() << "\t" << x_prim.getZ() << "\t"
			<< CoM.getX() << "\t" << CoM.getY() << "\t" << CoM.getZ() << "\t"
			<< theta << "\n";

	// Get rotation matrix
	Matrix rot_mat = getRotMatZ(theta);

	// Print info for each spring (and connected nodes)
	for (int i = 0; i < num_springs; i++) {

		// Get spring midpoint location
		Vector x = spr_mid(n_body_sim, springs[i], CoM);

		// Rotate around center of body in xy plane so that +x is toward perturber, +y is direction of rotation of perturber w.r.t to body
		// (So -y is headwind direction for body and +y is tailwind surface on body)
		Vector x_rot = rot_mat * x;

		// Get average power for this spring
		double power = tot_power[i];

		// Write average power applied to center of spring
		outfile << i << "\t" << std::setprecision(3) << x.getX() << "\t"
				<< x.getY() << "\t" << x.getZ() << "\t" << power * power_fac
				<< "\t" << x_rot.getX() << "\t" << x_rot.getY() << "\n";

		// Write heat on nodes as well as center of spring
		if (power_fac < 1.0) {
			// Get nodes spring is connected to
			int part_1 = springs[i].particle_1;
			int part_2 = springs[i].particle_2;

			// Get displacements of each particle from center of mass
			Vector x_1 = { particles[part_1].x, particles[part_1].y,
					particles[part_1].z };
			Vector dx_1 = x_1 - CoM;

			Vector x_2 = { particles[part_2].x, particles[part_2].y,
					particles[part_2].z };
			Vector dx_2 = x_2 - CoM;

			// Write info for each particle
			outfile << i << "\t" << std::setprecision(3) << dx_1.getX() << "\t"
					<< dx_1.getY() << "\t" << dx_1.getZ() << "\t"
					<< power * (1.0 - power_fac) * 0.5 << "\t" << x_rot.getX()
					<< "\t" << x_rot.getY() << "\n";
			outfile << i << "\t" << std::setprecision(3) << dx_2.getX() << "\t"
					<< dx_2.getY() << "\t" << dx_2.getZ() << "\t"
					<< power * (1.0 - power_fac) * 0.5 << "\t" << x_rot.getX()
					<< "\t" << x_rot.getY() << "\n";
		}
	}
	outfile.close();

	// Reset total power to start sum again
	for (int i = 0; i < num_springs; i++) {
		tot_power[i] = 0;
	}
}

/****************/
/* Heat helpers */
/****************/

// Normalize total power by number of timesteps for output
void normalize_tot_power(double ndt) {
	for (int i = 0; i < num_springs; i++) {
		tot_power[i] /= ndt;
	}
}

/********************/
/* Filename helpers */
/********************/

// Make a heat filename dependent on interval at which to print
string heat_filename(reb_simulation *const n_body_sim, string root,
		double print_interval) {
	// Get integer file number
	int xd = (int) (n_body_sim->t / print_interval);

	// Return filename
	return root + "_" + zero_pad_int(6, xd) + "_heat.txt";
}

// Make a node filename dependent on interval at which to print
string node_filename(reb_simulation *const n_body_sim, string root,
		double print_interval) {
	// Get integer file number
	int xd = (int) (n_body_sim->t / print_interval);

	// Return filename
	return root + "_" + zero_pad_int(6, xd) + "_node.txt";
}

// Make a stress filename dependent on interval at which to print
string stress_filename(reb_simulation *const n_body_sim, string root,
		double print_interval) {
	// Get integer file number
	int xd = (int) (n_body_sim->t / print_interval);

	// Return filename
	return root + "_" + zero_pad_int(6, xd) + "_stress.txt";
}
