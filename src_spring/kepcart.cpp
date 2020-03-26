/*
 keplerian conversion routines
 */

/* routines:
 keplerian: returns orbital elements from cartesian coords
 cartesian: returns cartesian from orbital elements
 ecc_ano:   solve Kepler's equation eccentric orbits
 ecc_anohyp:   solve Kepler's equation hyperbolic orbits
 */

#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include "kepcart.h"

#define SIGN(a) ((a) < 0 ? -1 : 1)

/* given PhaseState returns orbital elements */
void cart_to_kep(double GM, PhaseState state, OrbitalElements *orb_el) {
	double L_x, L_y, L_z, Lsqd, L;
	double r, vsqd, rdotv, rdot, ecostrueanom, esintrueanom, cosnode, sinnode;
	double rcosu, rsinu, u, trueanom, eccanom;

	/* find direction of angular momentum vector */
	L_x = state.y * state.vz - state.z * state.vy;
	L_y = state.z * state.vx - state.x * state.vz;
	L_z = state.x * state.vy - state.y * state.vx;
	Lsqd = L_x * L_x + L_y * L_y + L_z * L_z;
	L = sqrt(Lsqd);

	r = sqrt(state.x * state.x + state.y * state.y + state.z * state.z);
	vsqd = state.vx * state.vx + state.vy * state.vy + state.vz * state.vz;

	rdotv = state.x * state.vx + state.y * state.vy + state.z * state.vz;
	rdot = rdotv / r;

	orb_el->i = acos(L_z / L);

	if (L_x != 0.0 || L_y != 0.0) {
		orb_el->longnode = atan2(L_x, -L_y);
	} else
		orb_el->longnode = 0.0;

	orb_el->a = 1.0 / (2.0 / r - vsqd / GM); // could be negative

	ecostrueanom = Lsqd / (GM * r) - 1.0;
	esintrueanom = rdot * L / GM;
	orb_el->e = sqrt(ecostrueanom * ecostrueanom + esintrueanom * esintrueanom);

	if (esintrueanom != 0.0 || ecostrueanom != 0.0) {
		trueanom = atan2(esintrueanom, ecostrueanom);
	} else
		trueanom = 0.0;

	cosnode = cos(orb_el->longnode);
	sinnode = sin(orb_el->longnode);

	/* u is the argument of latitude */
	rcosu = state.x * cosnode + state.y * sinnode;
	rsinu = (state.y * cosnode - state.x * sinnode) / cos(orb_el->i);
	// potential divide by zero here!!!!!!!!

	if (rsinu != 0.0 || rcosu != 0.0) {
		u = atan2(rsinu, rcosu);
	} else
		u = 0.0;

	orb_el->argperi = u - trueanom;

	double foo = sqrt(fabs(1.0 - orb_el->e) / (1.0 + orb_el->e));
	if (orb_el->e < 1.0) {
		eccanom = 2.0 * atan(foo * tan(trueanom / 2.0));
		orb_el->meananom = eccanom - orb_el->e * sin(eccanom);
	} else {
		eccanom = 2.0 * atanh(foo * tan(trueanom / 2.0));
		orb_el->meananom = orb_el->e * sinh(eccanom) - eccanom;
	}
	if (orb_el->argperi > M_PI)
		orb_el->argperi -= 2.0 * M_PI;
	if (orb_el->argperi < -M_PI)
		orb_el->argperi += 2.0 * M_PI;
}

/* given orbital elements return PhaseState */
// hyperbolic orbits have e>1, not sure about convention for a positive or not
void kep_to_cart(double GM, OrbitalElements orb_el, PhaseState *state) {
	double meanmotion, cosE, sinE, foo;
	double x, y, z, vx, vy, vz;
	double xp, yp, zp, xdp, ydp, zdp;
	double cosw, sinw, cosi, sini, cosnode, sinnode;
	double E0, rovera;
	double a = orb_el.a;
	double e = orb_el.e;
	double i = orb_el.i;
	double longnode = orb_el.longnode;
	double argperi = orb_el.argperi;
	double meananom = orb_el.meananom;
	/* double E1, E2, den; */

	/* compute eccentric anomaly */

	if (e < 1)
		E0 = eccentric_anomaly(e, meananom);
	else
		E0 = eccentric_anomaly_hyperbolic(e, meananom);

	//  E0 = kepler(e,meananom); // also works

	if (e < 1.0) {
		cosE = cos(E0);
		sinE = sin(E0);
	} else {
		cosE = cosh(E0);
		sinE = sinh(E0);
	}
	a = fabs(a);
	meanmotion = sqrt(GM / (a * a * a));
	foo = sqrt(fabs(1.0 - e * e));
	/* compute unrotated positions and velocities */
	rovera = (1.0 - e * cosE);
	if (e > 1.0)
		rovera *= -1.0;
	x = a * (cosE - e);
	y = foo * a * sinE;
	z = 0.0;
	vx = -a * meanmotion * sinE / rovera;
	vy = foo * a * meanmotion * cosE / rovera;
	vz = 0.0;
	if (e > 1.0)
		x *= -1.0;

	/* rotate by argument of perihelion in orbit plane*/
	cosw = cos(argperi);
	sinw = sin(argperi);
	xp = x * cosw - y * sinw;
	yp = x * sinw + y * cosw;
	zp = z;
	xdp = vx * cosw - vy * sinw;
	ydp = vx * sinw + vy * cosw;
	zdp = vz;

	/* rotate by inclination about x axis */
	cosi = cos(i);
	sini = sin(i);
	x = xp;
	y = yp * cosi - zp * sini;
	z = yp * sini + zp * cosi;
	vx = xdp;
	vy = ydp * cosi - zdp * sini;
	vz = ydp * sini + zdp * cosi;

	/* rotate by longitude of node about z axis */
	cosnode = cos(longnode);
	sinnode = sin(longnode);
	state->x = x * cosnode - y * sinnode;
	state->y = x * sinnode + y * cosnode;
	state->z = z;
	state->vx = vx * cosnode - vy * sinnode;
	state->vy = vx * sinnode + vy * cosnode;
	state->vz = vz;
}

/* Helpers for coordinate conversion */

/* ----------------Solve Kepler's equation
 iterate to get an estimate for eccentric anomaly (u) given the mean anomaly (l).
 Appears to be accurate to level specified, I checked this
 and it works for u,lambda in all quadrants and at high eccentricity
 */
#define PREC_ecc_ano 1e-16  /* no reason that this must be very accurate in code at present */
double eccentric_anomaly(double ecc, double mean_anom) {
	double delta_u, u_guess, l_guess; // mean_anom_guess from the current guess of ecc anom
	delta_u = 1.0;
	u_guess = mean_anom + ecc * sin(mean_anom)
			+ 0.5 * pow(ecc, 2) * sin(2.0 * mean_anom);
	// also see M+D equation 2.55
	/* supposed to be good to second order in e, from Brouwer+Clemence
	 u0 is first guess */
	int counter = 0;
	while (fabs(delta_u) > PREC_ecc_ano) {
		l_guess = u_guess - ecc * sin(u_guess);
		delta_u = (mean_anom - l_guess) / (1.0 - ecc * cos(u_guess));
		u_guess += delta_u; /* this gives a better guess */
		counter++;
		if (counter > 10000)
			break;
		// equation 2.58 from M+D
	}
	return u_guess;
}

/* ----------------Solve Kepler's equation in hyperbolic case
 iterate to get an estimate for eccentric anomaly (u) given the mean anomaly (l).
 */
double eccentric_anomaly_hyperbolic(double ecc, double mean_anom) {
	double delta_u, u_guess, f_guess;
	delta_u = 1.0;
	u_guess = log(2.0 * mean_anom / ecc + 1.8); //danby guess
	int counter = 0;
	while (fabs(delta_u) > PREC_ecc_ano) {
		f_guess = u_guess - ecc * sinh(u_guess);
		delta_u = (-mean_anom - f_guess) / (1.0 - ecc * cosh(u_guess));
		u_guess += delta_u;
		counter++;
		if (counter > 10000)
			break;
	}
	return u_guess;
}
