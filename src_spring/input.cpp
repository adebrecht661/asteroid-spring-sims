/*
 * input.cpp
 *
 *  Created on: Apr 15, 2020
 *      Author: alex
 */

#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
extern "C" {
#include "rebound.h"
}
#include "springs.h"
#include "shapes.h" // For mindist - nowhere good for misc functions
#include "input.h"
using std::string;

extern int num_springs; // Total number of springs

/******************/
/* Input routines */
/******************/

// Read springs in from specified file
void read_springs(string fileroot, int index) {
	// Set filename
	string filename = fileroot + "_" + zero_pad_int(6, index) + "_springs.txt";

	// Read from file
	std::cout << "Reading in springs from " << filename << std::endl;
	std::ifstream spring_file(filename, std::ios::in);
	spring spr;
	while (spring_file.good()) {
		spring_file >> spr.particle_1;
		spring_file >> spr.particle_2;
		spring_file >> spr.k;
		spring_file >> spr.rs0;
		spring_file >> spr.gamma;
		spring_file >> spr.k_heat;
		// Add spring
		add_spring_helper(spr);
	}
	spring_file.close();
	std::cout << "Read " << num_springs << "springs." << std::endl;
}

// Read particles in from specified file
void read_particles(struct reb_simulation* const n_body_sim, string fileroot, int index) {
	// Set filename
	string filename = fileroot + "_" + zero_pad_int(6, index) + "_particles.txt";

	// Read from file
	std::cout << "Reading in particles from " << filename << std::endl;
	std::ifstream particle_file(filename, std::ios::in);
	struct reb_particle pt;
	pt.ax = 0.0;
	pt.ay = 0.0;
	pt.az = 0.0;
	while (particle_file.good()) {
		particle_file >> pt.x;
		particle_file >> pt.y;
		particle_file >> pt.z;
		particle_file >> pt.vx;
		particle_file >> pt.vy;
		particle_file >> pt.vz;
		particle_file >> pt.r;
		particle_file >> pt.m;
		// Add particle
		reb_add(n_body_sim, pt);
	}
	particle_file.close();
	std::cout << "Read " << n_body_sim->N << " particles." << std::endl;
}

// Read in a vertex (shape) file from filename
// File is ASCII in form "char x y z", where char is v if line defines a vertex and x, y, z are vertex positions in units of km
void read_vertex_file(struct reb_simulation* n_body_sim, string filename) {
	// Get current number of particles
	int i_low = n_body_sim->N;

	// Read vertex into particle structure
	struct reb_particle pt;
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
	while (vertex_file.good()) {
		vertex_file >> type;
		vertex_file >> pt.x;
		vertex_file >> pt.y;
		vertex_file >> pt.z;
		// If line is a vertex, add it
		if (type == 'v') {
			reb_add(n_body_sim, pt);
		}
	}
	vertex_file.close();

	// Get new number of particles
	int i_high = n_body_sim->N;
	std::cout << "Read " << i_high - i_low << " vertices from " << filename << std::endl;

	// Get smallest distance between particles
	double min_dist = mindist(n_body_sim, i_low, i_high);
	// Adjust radius of each vertex particle
	////// Why do we keep these particles???????
	for (int i = i_low; i < i_high; i++) {
		n_body_sim->particles[i].r = min_dist / 4.0;
	}
}

/***********/
/* Helpers */
/***********/

// Helper function to pad ints
string zero_pad_int(int width, int num)
{
    std::ostringstream ss;
    ss << std::setw( width ) << std::setfill( '0' ) << num;
    return ss.str();
}
