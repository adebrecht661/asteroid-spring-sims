/*
 * input.h
 *
 * Input routines to read in springs, particles, vertices from specified files
 *
 *  Created on: Apr 15, 2020
 *      Author: alex
 */

#ifndef SRC_SPRING_INPUT_H_
#define SRC_SPRING_INPUT_H_

#include <string>
using std::string;

// Read springs from fileroot_%06d_springs.txt
void read_springs(string fileroot, int index);
// Read particles from fileroot_%06d_particles.txt
void read_particles(struct reb_simulation* const n_body_sim, string fileroot, int index);
// Read vertices from filename
void read_vertex_file(struct reb_simulation* n_body_sim, string filename);
// Pad a number to width with zeros
string zero_pad_int(int width, int num);

#endif /* SRC_SPRING_INPUT_H_ */
