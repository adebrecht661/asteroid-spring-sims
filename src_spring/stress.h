/*
 * stress.h
 *
 *  Created on: Mar 26, 2020
 *      Author: alex
 */

#ifndef SRC_SPRING_STRESS_H_
#define SRC_SPRING_STRESS_H_

void spring_force_one(struct reb_simulation *const r, int i, double *Fx,
		double *Fy, double *Fz, double *Lx, double *Ly, double *Lz);

#endif /* SRC_SPRING_STRESS_H_ */
