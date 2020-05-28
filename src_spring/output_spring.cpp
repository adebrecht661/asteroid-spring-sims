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
#include "shapes.h"
#include "orb.h"
#include "input_spring.h" // For padding function - nowhere good for misc functions
#include "output_spring.h"

using std::string;
using std::vector;
using std::abs;

extern vector<spring> springs;
extern int num_springs;
extern int num_perts;
extern vector<stress_tensor> stresses;

/********************/
/* Output functions */
/********************/

// Write out positions and velocities of surface particles
// Can write at multiple times
void write_surf_part(reb_simulation *const n_body_sim, int i_low, int i_high,
		string filename, vector<bool> is_surf) {
	// Set width of each field
	int width = 12;

	// Get particle info
	reb_particle *particles = n_body_sim->particles;

	// Open output file to append or truncate, depending on if first run or not
	// File is created if it doesn't exist
	std::ofstream outfile;

	// Get time of first print for constant offset for future prints
	static bool first = true;
	static double t_off = 0.0;
	if (first) {
		outfile.open(filename, std::ios::out | std::ios::trunc);
		first = false;
		t_off = n_body_sim->t;
		std::cout << std::setprecision(6) << "print_surf: t_off = " << t_off
				<< "\n";

		// Print header only once
		outfile << "#" << std::setw(width - 1) << "t" << std::setw(width) << "x"
				<< std::setw(width) << "y" << std::setw(width) << "z"
				<< std::setw(width) << "vx" << std::setw(width) << "vy"
				<< std::setw(width) << "vz" << std::setw(width) << "ax"
				<< std::setw(width) << "ay" << std::setw(width) << "az"
				<< std::setw(width) << "r" << std::setw(width) << "m\n";
	} else {
		outfile.open(filename, std::ios::out | std::ios::app);
	}

	// Write actual current time
	outfile << std::setprecision(8) << "# t = " << n_body_sim->t << " #\n";

	int num_printed = 0;
	// Print tab-separated particle info at time t - t_off
	for (int i = 0; i < n_body_sim->N; i++) {
		if (is_surf[i]) {
			outfile << std::setprecision(6) << std::setw(width)
					<< n_body_sim->t - t_off << std::setw(width)
					<< particles[i].x << std::setw(width) << particles[i].y
					<< std::setw(width) << particles[i].z << std::setw(width)
					<< particles[i].vx << std::setw(width) << particles[i].vy
					<< std::setw(width) << particles[i].vz << std::setw(width)
					<< particles[i].ax << std::setw(width) << particles[i].ay
					<< std::setw(width) << particles[i].az << std::setw(width)
					<< particles[i].r << std::setw(width) << particles[i].m
					<< "\n";
			num_printed++;
		}
	}

	// Close file and report
	outfile.close();
	std::cout << "print_surf: Printed " << num_printed << " particles to file "
			<< filename << std::endl;
}

// Write out information for an extended body with indices [i_low, i_high)
void write_resolved_with_E(reb_simulation *const n_body_sim, int i_low,
		int i_high, string filename, double dEdt_ave) {
	// Set width of each field
	int width = 14;

	// Open output file to append (created if it doesn't exist)
	std::ofstream outfile;

	// Print header once
	static bool first = true;
	if (first) {
		first = false;
		outfile.open(filename, std::ios::out | std::ios::trunc);
		outfile << "#" << std::setw(width - 1) << "t" << std::setw(width) << "x"
				<< std::setw(width) << "y" << std::setw(width) << "z"
				<< std::setw(width) << "vx" << std::setw(width) << "vy"
				<< std::setw(width) << "vz" << std::setw(width) << "omx"
				<< std::setw(width) << "omy" << std::setw(width) << "omz"
				<< std::setw(width) << "llx" << std::setw(width) << "lly"
				<< std::setw(width) << "llz" << std::setw(width) << "Ixx"
				<< std::setw(width) << "Iyy" << std::setw(width) << "Izz"
				<< std::setw(width) << "Ixy" << std::setw(width) << "Iyz"
				<< std::setw(width) << "Ixz" << std::setw(width) << "KErot"
				<< std::setw(width) << "PEspr" << std::setw(width) << "PEgrav"
				<< std::setw(width) << "Etot" << std::setw(width) << "dEdtnow"
				<< std::setw(width) << "dEdtave" << "\n";
	} else {
		outfile.open(filename, std::ios::out | std::ios::app);
	}

	// Time
	outfile << std::setprecision(3) << std::setw(width) << n_body_sim->t;

	// Center of mass location and velocity
	Vector CoM = compute_com(n_body_sim, i_low, i_high);
	Vector CoV = compute_cov(n_body_sim, i_low, i_high);
	outfile << std::setprecision(5) << std::setw(width) << CoM.getX()
			<< std::setw(width) << CoM.getY() << std::setw(width) << CoM.getZ();
	outfile << std::setprecision(5) << std::setw(width) << CoV.getX()
			<< std::setw(width) << CoV.getY() << std::setw(width) << CoV.getZ();

	// Spin angular momentum
	Vector omega = body_spin(n_body_sim, i_low, i_high);
	Vector L = measure_L(n_body_sim, i_low, i_high);
	outfile << std::setprecision(4) << std::setw(width) << omega.getX()
			<< std::setw(width) << omega.getY() << std::setw(width)
			<< omega.getZ();
	outfile << std::setprecision(4) << std::setw(width) << L.getX()
			<< std::setw(width) << L.getY() << std::setw(width) << L.getZ();

	// Moment of inertia
	Matrix inertia_mat = mom_inertia(n_body_sim, i_low, i_high);
	outfile << std::setprecision(5) << std::setw(width) << inertia_mat.getXX()
			<< std::setw(width) << inertia_mat.getYY() << std::setw(width)
			<< inertia_mat.getZZ() << std::setw(width) << inertia_mat.getXY()
			<< std::setw(width) << inertia_mat.getYZ() << std::setw(width)
			<< inertia_mat.getXZ();

	// Non-translational kinetic energy
	double KErot = compute_rot_kin(n_body_sim, i_low, i_high);
	outfile << std::setprecision(5) << std::setw(width) << KErot;

	// Spring potential energy
	double SPE = spring_potential_energy(n_body_sim);
	outfile << std::setprecision(5) << std::setw(width) << SPE;

	// Gravitational potential energy
	double GPE = grav_potential_energy(n_body_sim, i_low, i_high);
	outfile << std::setprecision(10) << std::setw(width) << GPE;

	// Total internal energy
	outfile << std::setprecision(10) << std::setw(width) << (KErot + SPE + GPE);

	// Total power dissipated at this time by spring damping
	outfile << std::setprecision(10) << std::setw(width)
			<< dEdt_total(n_body_sim);

	// Average power dissipation by spring damping
	outfile << std::setprecision(10) << std::setw(width) << dEdt_ave << "\n";
	outfile.close();
}

// Write out information for an extended body with indices [i_low, i_high)
// No energy terms
void write_resolved_no_E(reb_simulation *const n_body_sim, int i_low,
		int i_high, string filename) {
	// Set width of each field
	int width = 12;

	// Open output file to append (created if it doesn't exist)
	std::ofstream outfile;

	// Print header once
	static bool first = true;
	if (first) {
		first = false;
		outfile.open(filename, std::ios::out | std::ios::trunc);
		outfile << "#" << std::setw(width - 1) << "t" << std::setw(width) << "x"
				<< std::setw(width) << "y" << std::setw(width) << "z"
				<< std::setw(width) << "vx" << std::setw(width) << "vy"
				<< std::setw(width) << "vz" << std::setw(width) << "omx"
				<< std::setw(width) << "omy" << std::setw(width) << "omz"
				<< std::setw(width) << "llx" << std::setw(width) << "lly"
				<< std::setw(width) << "llz" << std::setw(width) << "Ixx"
				<< std::setw(width) << "Iyy" << std::setw(width) << "Izz"
				<< std::setw(width) << "Ixy" << std::setw(width) << "Iyz"
				<< std::setw(width) << "Ixz" << "\n";
	} else {
		outfile.open(filename, std::ios::out | std::ios::app);
	}

	// Time
	outfile << std::setprecision(3) << std::setw(width) << n_body_sim->t;

	// Center of mass location and velocity
	Vector CoM = compute_com(n_body_sim, i_low, i_high);
	Vector CoV = compute_cov(n_body_sim, i_low, i_high);
	outfile << std::setprecision(5) << std::setw(width) << CoM.getX()
			<< std::setw(width) << CoM.getY() << std::setw(width) << CoM.getZ();
	outfile << std::setprecision(5) << std::setw(width) << CoV.getX()
			<< std::setw(width) << CoV.getY() << std::setw(width) << CoV.getZ();

	// Spin angular momentum
	Vector omega = body_spin(n_body_sim, i_low, i_high);
	Vector L = measure_L(n_body_sim, i_low, i_high);
	outfile << std::setprecision(4) << std::setw(width) << omega.getX()
			<< std::setw(width) << omega.getY() << std::setw(width)
			<< omega.getZ();
	outfile << std::setprecision(4) << std::setw(width) << L.getX()
			<< std::setw(width) << L.getY() << std::setw(width) << L.getZ();

	// Moment of inertia
	Matrix inertia_mat = mom_inertia(n_body_sim, i_low, i_high);
	outfile << std::setprecision(5) << std::setw(width) << inertia_mat.getXX()
			<< std::setw(width) << inertia_mat.getYY() << std::setw(width)
			<< inertia_mat.getZZ() << std::setw(width) << inertia_mat.getXY()
			<< std::setw(width) << inertia_mat.getYZ() << std::setw(width)
			<< inertia_mat.getXZ() << "\n";
	outfile.close();
}

// Write out information for two specific nodes in an extended body with indices [i_low, i_high)
// The first node has largest x value at initial time
// The second node has largest z value at initial time
// Same nodes are printed each time after they are chosen
void write_resolved_2nodes(reb_simulation *const n_body_sim, int i_low,
		int i_high, string filename) {
	// Set width of each field
	int width = 12;

	// Get particle info
	reb_particle *particles = n_body_sim->particles;

	// Indices of particles with largest z and x values
	static int i_x = 0;
	static int i_z = 0;

	// Open output file to append (created if it doesn't exist)
	std::ofstream outfile;

	// Write header once
	static bool first = true;
	if (first) {
		first = false;
		outfile.open(filename, std::ios::out | std::ios::trunc);
		outfile << "#" << std::setw(width - 1) << "#" << std::setw(width * 3)
				<< "xnode" << std::setw(width * 3) << "#"
				<< std::setw(width * 3) << "znode" << std::setw(width * 3)
				<< "#" << "\n";
		outfile << "#" << std::setw(width - 1) << "t" << std::setw(width) << "x"
				<< std::setw(width) << "y" << std::setw(width) << "z"
				<< std::setw(width) << "vx" << std::setw(width) << "vy"
				<< std::setw(width) << "vz" << std::setw(width) << "x"
				<< std::setw(width) << "y" << std::setw(width) << "z"
				<< std::setw(width) << "vx" << std::setw(width) << "vy"
				<< std::setw(width) << "vz" << "\n";

		// Get particles with largest x and z values (points should be well outside the system)
		i_x = nearest_to_point(n_body_sim, i_low, i_high, Vector( { 10.0, 0.0,
				0.0 }));
		i_z = nearest_to_point(n_body_sim, i_low, i_high, Vector( { 0.0, 0.0,
				10.0 }));
	} else {
		outfile.open(filename, std::ios::out | std::ios::app);
	}

	// Time
	outfile << std::setprecision(3) << std::setw(width) << n_body_sim->t;

	// x particle location and velocity
	outfile << std::setprecision(5) << std::setw(width) << particles[i_x].x
			<< std::setw(width) << particles[i_x].y << std::setw(width)
			<< particles[i_x].z << std::setw(width) << particles[i_x].vx
			<< std::setw(width) << particles[i_x].vy << std::setw(width)
			<< particles[i_x].vz;

	// z particle location and velocity
	outfile << std::setprecision(5) << std::setw(width) << particles[i_z].x
			<< std::setw(width) << particles[i_z].y << std::setw(width)
			<< particles[i_z].z << std::setw(width) << particles[i_z].vx
			<< std::setw(width) << particles[i_z].vy << std::setw(width)
			<< particles[i_z].vz << "\n";
	outfile.close();
}

// Write out information for a point mass with particle index i_pert
// pert_num is count of perturbing particle
void write_pt_mass(reb_simulation *const n_body_sim, int i_pert, int pert_num,
		string filename) {
	// Set width of each field
	int width = 12;

	// Get particle info
	reb_particle *particles = n_body_sim->particles;

	// Record whether it's the first run for each particle
	const int max_perts = 100;
	static bool first[max_perts];

	// Open output file to append (created if it doesn't exist)
	std::ofstream outfile;

	// Initialize boolean array
	static bool first_run = true;
	if (first_run) {
		outfile.open(filename, std::ios::out | std::ios::trunc);
		for (int i = 0; i < num_perts; i++) {
			first[i] = true;
		}
	} else {
		outfile.open(filename, std::ios::out | std::ios::app);
	}

	// Write header once per particle
	if (first[pert_num]) {
		first[pert_num] = false;
		outfile << "#" << std::setw(width - 1) << "t" << std::setw(width) << "x"
				<< std::setw(width) << "y" << std::setw(width) << "z"
				<< std::setw(width) << "vx" << std::setw(width) << "vy"
				<< std::setw(width) << "vz" << std::setw(width) << "m" << "\n";
	}

	// Time
	outfile << std::setprecision(3) << std::setw(width) << n_body_sim->t;

	// Particle info
	outfile << std::setprecision(5) << std::setw(width) << particles[i_pert].x
			<< std::setw(width) << particles[i_pert].y << std::setw(width)
			<< particles[i_pert].z << std::setw(width) << particles[i_pert].vx
			<< std::setw(width) << particles[i_pert].vy << std::setw(width)
			<< particles[i_pert].vz;
	outfile << std::setprecision(3) << std::setw(width) << particles[i_pert].m
			<< "\n";
	outfile.close();
}

// Write information about resolved body orbiting one or more point particles
void write_resolved_orb(reb_simulation *const n_body_sim, string filename) {
	// Set width of each field
	int width = 12;

	// Open output file to append (created if it doesn't exist)
	std::ofstream outfile;

	// Write header once
	static bool first = true;
	if (first) {
		first = false;
		outfile.open(filename, std::ios::out | std::ios::trunc);
		outfile << "#" << std::setw(width - 1) << "t" << std::setw(width) << "a"
				<< std::setw(width) << "n" << std::setw(width) << "e"
				<< std::setw(width) << "i" << std::setw(width) << "omx"
				<< std::setw(width) << "omy" << std::setw(width) << "omz"
				<< std::setw(width) << "A" << std::setw(width) << "B"
				<< std::setw(width) << "C" << std::setw(width) << "E"
				<< std::setw(width) << "lx" << std::setw(width) << "ly"
				<< std::setw(width) << "lz" << std::setw(width) << "ang"
				<< std::setw(width) << "lox" << std::setw(width) << "loy"
				<< std::setw(width) << "loz" << std::setw(width) << "Ixx"
				<< std::setw(width) << "Iyy" << std::setw(width) << "Izz"
				<< std::setw(width) << "Ixy" << std::setw(width) << "Iyz"
				<< std::setw(width) << "Ixz" << "\n";
	} else {
		outfile.open(filename, std::ios::out | std::ios::app);
	}

	// Get indices of resolved body
	int i_low = 0;
	int i_high = n_body_sim->N - num_perts;

	// Time
	outfile << std::setprecision(3) << std::setw(width) << n_body_sim->t;

	// Orbital properties of resolved body around perturbers
	double a, mean_motion, ecc, incl, Lorb;
	compute_orb(n_body_sim, i_low, i_high, &a, &mean_motion, &ecc, &incl,
			&Lorb);
	outfile << std::setprecision(6) << std::setw(width) << a;
	outfile << std::setprecision(4) << std::setw(width) << mean_motion
			<< std::setw(width) << ecc << std::setw(width) << incl;

	// Spin
	Vector omega = body_spin(n_body_sim, i_low, i_high);
	outfile << std::setprecision(4) << std::setw(width) << omega.getX()
			<< std::setw(width) << omega.getY() << std::setw(width)
			<< omega.getZ();

	// Get moment of inertia eigenvalues
	Matrix inertia_mat = mom_inertia(n_body_sim, i_low, i_high);
	double eigs[3];
	eigenvalues(inertia_mat, eigs);
	outfile << std::setprecision(3) << std::setw(width) << eigs[0]
			<< std::setw(width) << eigs[1] << std::setw(width) << eigs[2];

	// Young's modulus
	double r_min = 0.0, r_max = 0.5;
	outfile << std::setprecision(3) << std::setw(width)
			<< Young_mesh(n_body_sim, i_low, i_high, r_min, r_max);

	// Spin angular momentum
	Vector L = measure_L(n_body_sim, i_low, i_high);
	outfile << std::setprecision(5) << std::setw(width) << L.getX()
			<< std::setw(width) << L.getY() << std::setw(width) << L.getZ();

	// Get orbital angular momentum
	Vector L_o = compute_Lorb(n_body_sim, i_low, i_high);

	// Angle between spin and orbit
	outfile << std::setprecision(5) << std::setw(width)
			<< acos(dot(L, L_o) / (L.len() * L_o.len()));

	// Output orbital angular momentum
	outfile << std::setprecision(5) << std::setw(width) << L_o.getX()
			<< std::setw(width) << L_o.getY() << std::setw(width) << L_o.getZ();

	// Print moment of inertia
	outfile << std::setprecision(5) << std::setw(width) << inertia_mat.getXX()
			<< std::setw(width) << inertia_mat.getYY() << std::setw(width)
			<< inertia_mat.getZZ() << std::setw(width) << inertia_mat.getXY()
			<< std::setw(width) << inertia_mat.getYZ() << std::setw(width)
			<< inertia_mat.getXZ() << "\n";
	outfile.close();
}

// Write out resolved body and point particle info
void write_resolved_bin(reb_simulation *const n_body_sim, string filename) {
	// Set width of each field
	int width = 12;

	// Get particle info
	reb_particle *particles = n_body_sim->particles;

	// Get particle with largest radius
	static int part_largest_r = -1;

	// Open output file to append (created if it doesn't exist)
	std::ofstream outfile;

	// Write header once
	static bool first = true;
	if (first) {
		first = false;
		outfile.open(filename, std::ios::out | std::ios::trunc);
		outfile << "#" << std::setw(width - 1) << "#" << std::setw(width * 3)
				<< "resolved body to perts" << std::setw(width * 3) << "#"
				<< std::setw(width * 3) << "largest radius to resolved"
				<< std::setw(width * 3) << "#" << std::setw(width * 3)
				<< "1st pert to 2nd pert" << std::setw(width * 3) << "#"
				<< std::setw(width * 4) << "resolved body\n";
		outfile << "#" << std::setw(width - 1) << "t" << std::setw(width) << "x"
				<< std::setw(width) << "y" << std::setw(width) << "z"
				<< std::setw(width) << "vx" << std::setw(width) << "vy"
				<< std::setw(width) << "vz" << std::setw(width) << "x"
				<< std::setw(width) << "y" << std::setw(width) << "z"
				<< std::setw(width) << "vx" << std::setw(width) << "vy"
				<< std::setw(width) << "vz" << std::setw(width) << "x"
				<< std::setw(width) << "y" << std::setw(width) << "z"
				<< std::setw(width) << "vx" << std::setw(width) << "vy"
				<< std::setw(width) << "vz" << std::setw(width) << "lx"
				<< std::setw(width) << "ly" << std::setw(width) << "lz"
				<< std::setw(width) << "Ixx" << std::setw(width) << "Iyy"
				<< std::setw(width) << "Izz" << std::setw(width) << "Ixy"
				<< std::setw(width) << "Iyz" << std::setw(width) << "Ixz\n";

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
	} else {
		outfile.open(filename, std::ios::out | std::ios::app);
	}

	// Indices of resolved body
	int i_low = 0;
	int i_high = n_body_sim->N - num_perts;

	// Indices of perturbing masses
	int i_pert_low = n_body_sim->N - num_perts;
	int i_pert_high = n_body_sim->N;

	// Time
	outfile << std::setprecision(3) << std::setw(width) << n_body_sim->t;

	// Center of mass of resolved body location and velocity
	Vector CoM = compute_com(n_body_sim, i_low, i_high);
	Vector CoV = compute_cov(n_body_sim, i_low, i_high);

	// Center of mass of perturbing particles location and velocity
	Vector CoM_pert = compute_com(n_body_sim, i_pert_low, i_pert_high);
	Vector CoV_pert = compute_cov(n_body_sim, i_pert_low, i_pert_high);

	// Relative position and velocity of resolved body to center of perturbers
	Vector xr = CoM - CoM_pert;
	Vector vr = CoV - CoV_pert;
	outfile << std::setprecision(5) << std::setw(width) << xr.getX()
			<< std::setw(width) << xr.getY() << std::setw(width) << xr.getZ()
			<< std::setw(width) << vr.getX() << std::setw(width) << vr.getY()
			<< std::setw(width) << vr.getZ();

	// Relative position and velocity of particle with largest radius to resolved body
	Vector x_largest = { particles[part_largest_r].x,
			particles[part_largest_r].y, particles[part_largest_r].z };
	Vector v_largest = { particles[part_largest_r].vx,
			particles[part_largest_r].vy, particles[part_largest_r].vz };
	Vector xs = x_largest - CoM;
	Vector vs = v_largest - CoV;
	outfile << std::setprecision(5) << std::setw(width) << xs.getX()
			<< std::setw(width) << xs.getY() << std::setw(width) << xs.getZ()
			<< std::setw(width) << vs.getX() << std::setw(width) << vs.getY()
			<< std::setw(width) << vs.getZ();

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

		outfile << std::setprecision(5) << std::setw(width) << x_rel.getX()
				<< std::setw(width) << x_rel.getY() << std::setw(width)
				<< x_rel.getZ() << std::setw(width) << v_rel.getX()
				<< std::setw(width) << v_rel.getY() << std::setw(width)
				<< v_rel.getZ();
	} else {
		outfile << std::setw(width) << 0 << std::setw(width) << 0
				<< std::setw(width) << 0 << std::setw(width) << 0
				<< std::setw(width) << 0 << std::setw(width) << 0;
	}

	// Spin angular momentum
	Vector L = measure_L(n_body_sim, i_low, i_high);
	outfile << std::setprecision(5) << std::setw(width) << L.getX()
			<< std::setw(width) << L.getY() << std::setw(width) << L.getZ();

	// Moment of inertia
	Matrix inertia_mat = mom_inertia(n_body_sim, i_low, i_high);
	outfile << std::setprecision(5) << std::setw(width) << inertia_mat.getXX()
			<< std::setw(width) << inertia_mat.getYY() << std::setw(width)
			<< inertia_mat.getZZ() << std::setw(width) << inertia_mat.getXY()
			<< std::setw(width) << inertia_mat.getYZ() << std::setw(width)
			<< inertia_mat.getXZ() << "\n";
	outfile.close();
}

// Write out springs to a file
void write_springs(reb_simulation *const n_body_sim, string fileroot,
		int index) {
	// Set width of each field
	int width = 12;

	// Set filename
	string filename = fileroot + "_" + zero_pad_int(6, index) + "_springs.txt";

	// Open output file to append (created if it doesn't exist)
	std::ofstream outfile;
	static bool first = true;
	if (first) {
		first = false;
		outfile.open(filename, std::ios::out | std::ios::trunc);
	} else {
		outfile.open(filename, std::ios::out | std::ios::app);
	}

	// Write out each spring
	for (int i = 0; i < num_springs; i++) {
		outfile << std::setw(width) << springs[i].particle_1 << std::setw(width)
				<< springs[i].particle_2 << std::setw(width)
				<< std::setprecision(6) << springs[i].k << std::setw(width)
				<< springs[i].rs0 << std::setw(width) << springs[i].gamma
				<< std::setw(width) << strain(n_body_sim, springs[i])
				<< std::setw(width) << spring_r(n_body_sim, springs[i]).len()
				<< "\n";
	}
	outfile.close();
	std::cout << "write_springs: Printed " << num_springs << " springs to "
			<< filename << std::endl;
}

// Write out particles to file
void write_particles(reb_simulation *const n_body_sim, string fileroot,
		int index) {
	// Set width of each field
	int width = 12;

	// Get particle info
	reb_particle *particles = n_body_sim->particles;

	// Set filename
	string filename = fileroot + "_" + zero_pad_int(6, index)
			+ "_particles.txt";

	// Open output file to append (created if it doesn't exist)
	std::ofstream outfile;
	static bool first = true;
	if (first) {
		first = false;
		outfile.open(filename, std::ios::out | std::ios::trunc);
	} else {
		outfile.open(filename, std::ios::out | std::ios::app);
	}

	// Write out each particle
	for (int i = 0; i < n_body_sim->N; i++) {
		outfile << std::setprecision(6) << std::setw(width) << particles[i].x
				<< std::setw(width) << particles[i].y << std::setw(width)
				<< particles[i].z << std::setw(width) << particles[i].vx
				<< std::setw(width) << particles[i].vy << std::setw(width)
				<< particles[i].vz << std::setw(width) << particles[i].r
				<< std::setw(width) << particles[i].m << "\n";
	}
	outfile.close();
	std::cout << "write_particles: Printed " << n_body_sim->N
			<< " particles to " << filename << std::endl;
}

// Write out stress file
void write_stresses(reb_simulation *const n_body_sim, string filename) {
	// Set width of each field
	int width = 12;

	// Get particle info
	reb_particle *particles = n_body_sim->particles;

	// Open output file to append (created if it doesn't exist)
	std::ofstream outfile;

	// Print header first time only
	static bool first = true;
	if (first) {
		first = false;
		outfile.open(filename, std::ios::out | std::ios::trunc);
		outfile << "#" << std::setw(width - 1) << "i" << std::setw(width) << "m"
				<< std::setw(width) << "x" << std::setw(width) << "y"
				<< std::setw(width) << "z" << std::setw(width) << "eig1"
				<< std::setw(width) << "eig2" << std::setw(width) << "eig3"
				<< std::setw(width) << "failing" << "\n";
	} else {
		outfile.open(filename, std::ios::out | std::ios::app);
	}

	// Time
	outfile << "# t = " << std::setprecision(2) << n_body_sim->t << " #\n";

	// Get resolved body indices
	int i_low = 0;
	int i_high = n_body_sim->N - num_perts;

	// Get center of mass of resolved body
	Vector CoM = compute_com(n_body_sim, i_low, i_high);

	// Print stress of each particle
	for (int i = i_low; i < i_high; i++) {
		Vector x = { particles[i].x, particles[i].y, particles[i].z };
		Vector dx = x - CoM;
		outfile << std::setw(width) << i << std::setw(width)
				<< std::setprecision(4) << particles[i].m << std::setw(width)
				<< std::setprecision(3) << dx.getX() << std::setw(width)
				<< dx.getY() << std::setw(width) << dx.getZ();
		outfile << std::setprecision(3) << std::setw(width)
				<< stresses[i].eigs[0] << std::setw(width)
				<< stresses[i].eigs[1] << std::setw(width)
				<< stresses[i].eigs[2];
		outfile << std::setw(width) << stresses[i].failing << "\n";
	}
	outfile.close();
}

// Print doubles to file and standard out
void print_run_double(double quantity, string label, std::ofstream *outfile) {
	// Set precision based on size of quantity
	if ((abs(log10(quantity)) > 4)) {
		std::cout << label << std::setprecision(4) << " " << quantity << "\n";
		*outfile << label << std::setprecision(4) << " " << quantity << "\n";
	} else {
		std::cout << label << std::setprecision(3) << " " << quantity << "\n";
		*outfile << label << std::setprecision(3) << " " << quantity << "\n";
	}
}

/********************/
/* Filename helpers */
/********************/

// Make a stress filename dependent on interval at which to print
string stress_filename(reb_simulation *const n_body_sim, string root,
		double print_interval) {
	// Get integer file number
	int xd = (int) (n_body_sim->t / print_interval);

	// Return filename
	return root + "_" + zero_pad_int(6, xd) + "_stress.txt";
}
