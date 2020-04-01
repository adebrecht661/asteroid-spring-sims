/*
 * stress.h
 *
 *  Created on: Mar 26, 2020
 *      Author: alex
 */

#ifndef SRC_SPRING_STRESS_H_
#define SRC_SPRING_STRESS_H_

// Stress tensor for each node
typedef struct stress_tensor {
	// Stress tensor
	double sigxx;
	double sigyy;
	double sigzz;
	double sigxy;
	double sigyz;
	double sigxz;

	// Eigenvalues of stress tensor, from largest to smallest
	double eig1;
	double eig2;
	double eig3;

	double maxF;  // max Force of a spring
	int s_index;  // index of spring giving max Force
	bool fail;  // is there material failure
} stress_tensor;

void spring_force_one(struct reb_simulation *const r, int i, double *Fx,
		double *Fy, double *Fz, double *Lx, double *Ly, double *Lz);

#endif /* SRC_SPRING_STRESS_H_ */
