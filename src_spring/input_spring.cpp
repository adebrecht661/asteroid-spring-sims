#ifdef __cplusplus
# 	ifdef __GNUC__
#		define restrict __restrict__
#	else
#		define restrict
#	endif
#endif

/*
 * input.cpp
 *
 *  Created on: Apr 15, 2020
 *      Author: alex
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
#include "springs.h"
#include "shapes.h" // For mindist - nowhere good for misc functions
#include "input_spring.h"

using namespace libconfig;
using std::string;
using std::vector;

extern int num_springs;		// Total number of springs

// Global scales
extern double mass_scale, time_scale, length_scale, temp_scale, omega_scale,
		vel_scale, p_scale, L_scale, a_scale, F_scale, E_scale, dEdt_scale,
		P_scale;

/******************/
/* Input routines */
/******************/

// Read scales from problem.cfg and output to scales file
void read_scales(Config *cfg) {
	// Read and output scale info
	if (!(cfg->lookupValue("mass_scale", mass_scale)
			&& cfg->lookupValue("time_scale", time_scale)
			&& cfg->lookupValue("length_scale", length_scale)
			&& cfg->lookupValue("temp_scale", temp_scale))) {
		throw "Error reading scales from problem.cfg. Exiting.";
	}
	omega_scale = 1.0 / time_scale;
	vel_scale = length_scale / time_scale;
	p_scale = mass_scale * vel_scale;
	L_scale = length_scale * p_scale;
	a_scale = vel_scale / time_scale;
	F_scale = a_scale * mass_scale;
	E_scale = F_scale * length_scale;
	dEdt_scale = E_scale / time_scale;
	P_scale = F_scale / pow(length_scale, 2.0);

	std::ofstream scalefile("scales.txt", std::ios::out | std::ios::trunc);
	scalefile << "Multiply values by these scales to get physical values:\n"
			<< "Mass scale: " << mass_scale << "\nTime scale: " << time_scale
			<< "\nLength scale: " << length_scale << "\nTemperature scale: "
			<< temp_scale << "\nAngular velocity scale: " << omega_scale
			<< "\nVelocity scale: " << vel_scale << "\nMomentum scale: "
			<< p_scale << "\nAngular momentum scale: " << L_scale
			<< "\nAcceleration scale: " << a_scale << "\nForce scale: "
			<< F_scale << "\nEnergy scale: " << E_scale << "\nPower scale: "
			<< dEdt_scale << "\nPressure scale: " << P_scale;
}

// Read springs in from specified file
void read_springs(string fileroot, int index) {
	// Set filename
	string filename = fileroot + "_" + zero_pad_int(6, index) + "_springs.txt";

	// Read from file
	std::cout << "Reading in springs from " << filename << std::endl;
	std::ifstream spring_file(filename, std::ios::in);
	spring spr;
	string line;

	// Each line is a different spring
	while (std::getline(spring_file, line)) {
		std::istringstream input_stream(line);

		if (!(input_stream >> spr.particle_1 >> spr.particle_2 >> spr.k
				>> spr.rs0 >> spr.gamma >> spr.k_heat))
			throw "I/O error in read_springs.";

		// Add spring
		try {
			add_spring_helper(spr);
		} catch (char *str) {
			std::cerr << str << "Exiting." << std::endl;
			exit(1);
		}
	}
	spring_file.close();
	std::cout << "Read " << num_springs << "springs." << std::endl;
}

// Read particles in from specified file
void read_particles(reb_simulation *const n_body_sim, string fileroot,
		int index) {
	// Set filename
	string filename = fileroot + "_" + zero_pad_int(6, index)
			+ "_particles.txt";

	// Read from file
	std::cout << "Reading in particles from " << filename << std::endl;
	std::ifstream particle_file(filename, std::ios::in);
	reb_particle pt;
	pt.ax = 0.0;
	pt.ay = 0.0;
	pt.az = 0.0;
	string line;

	// Each line is a different particle
	while (std::getline(particle_file, line)) {
		std::istringstream input_stream(line);

		if (!(input_stream >> pt.x >> pt.y >> pt.z >> pt.vx >> pt.vy >> pt.vz
				>> pt.r >> pt.m))
			throw "I/O error in read_particles.";

		// Add particle
		reb_add(n_body_sim, pt);
	}
	particle_file.close();
	std::cout << "Read " << n_body_sim->N << " particles." << std::endl;
}

// Read in a vertex (shape) file from filename
// File is ASCII in form "char x y z", where char is v if line defines a vertex and x, y, z are vertex positions in units of km
void read_vertices(reb_simulation *n_body_sim, string filename) {
	// Get current number of particles
	int i_low = n_body_sim->N;

	// Read vertex into particle structure
	reb_particle pt;
	std::ifstream vertex_file(filename, std::ios::in);
	pt.ax = 0.0;
	pt.ay = 0.0;
	pt.az = 0.0;
	pt.vx = 0.0;
	pt.vy = 0.0;
	pt.vz = 0.0;
	pt.m = 1.0; // arbitrary mass
	pt.r = 1.0; // arbitrary radius
	char type;
	string line;

	// Each line is (possibly) a different vertex
	while (std::getline(vertex_file, line)) {
		std::istringstream input_stream(line);

		if (!(input_stream >> type >> pt.x >> pt.y >> pt.z))
			throw "I/O error in read_vertices.";

		// If line is a vertex, add it
		if (type == 'v') {
			reb_add(n_body_sim, pt);
		}
	}
	vertex_file.close();

	// Get new number of particles
	int i_high = n_body_sim->N;
	std::cout << "Read " << i_high - i_low << " vertices from " << filename
			<< std::endl;

	// Get smallest distance between particles
	double min_dist = mindist(n_body_sim, i_low, i_high);

	// Adjust radius of each vertex particle
	for (int i = i_low; i < i_high; i++) {
		n_body_sim->particles[i].r = min_dist / 4.0;
	}
}

/***********/
/* Helpers */
/***********/

// Helper function to pad ints
string zero_pad_int(int width, int num) {
	std::ostringstream ss;
	ss << std::setw(width) << std::setfill('0') << num;
	return ss.str();
}
