/*
 keplerian conversion routines
 */

/* routines:
 keplerian: returns orbital elements from cartesian coords
 cartesian: returns cartesian from orbital elements
 ecc_ano:   solve Kepler's equation eccentric orbits
 ecc_anohyp:   solve Kepler's equation hyperbolic orbits
 */

#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include "kepcart.h"

#define SIGN(a) ((a) < 0 ? -1 : 1)

/* given PhaseState returns orbital elements */
void cart_to_kep(double GM, PhaseState state, OrbitalElements *orbel) {
	double rxv_x, rxv_y, rxv_z, hs, h;
	double r, vs, rdotv, rdot, ecostrueanom, esintrueanom, cosnode, sinnode;
	double rcosu, rsinu, u, trueanom, eccanom;

	/* find direction of angular momentum vector */
	rxv_x = state.y * state.vz - state.z * state.vy;
	rxv_y = state.z * state.vx - state.x * state.vz;
	rxv_z = state.x * state.vy - state.y * state.vx;
	hs = rxv_x * rxv_x + rxv_y * rxv_y + rxv_z * rxv_z;
	h = sqrt(hs);

	r = sqrt(state.x * state.x + state.y * state.y + state.z * state.z);
	vs = state.vx * state.vx + state.vy * state.vy + state.vz * state.vz;

	rdotv = state.x * state.vx + state.y * state.vy + state.z * state.vz;
	rdot = rdotv / r;

	orbel->i = acos(rxv_z / h);

	if (rxv_x != 0.0 || rxv_y != 0.0) {
		orbel->longnode = atan2(rxv_x, -rxv_y);
	} else
		orbel->longnode = 0.0;

	orbel->a = 1.0 / (2.0 / r - vs / GM); // could be negative

	ecostrueanom = hs / (GM * r) - 1.0;
	esintrueanom = rdot * h / GM;
	orbel->e = sqrt(ecostrueanom * ecostrueanom + esintrueanom * esintrueanom);

	if (esintrueanom != 0.0 || ecostrueanom != 0.0) {
		trueanom = atan2(esintrueanom, ecostrueanom);
	} else
		trueanom = 0.0;

	cosnode = cos(orbel->longnode);
	sinnode = sin(orbel->longnode);

	/* u is the argument of latitude */
	rcosu = state.x * cosnode + state.y * sinnode;
	rsinu = (state.y * cosnode - state.x * sinnode) / cos(orbel->i);
	// potential divide by zero here!!!!!!!!

	if (rsinu != 0.0 || rcosu != 0.0) {
		u = atan2(rsinu, rcosu);
	} else
		u = 0.0;

	orbel->argperi = u - trueanom;

	double foo = sqrt(fabs(1.0 - orbel->e) / (1.0 + orbel->e));
	if (orbel->e < 1.0) {
		eccanom = 2.0 * atan(foo * tan(trueanom / 2.0));
		orbel->meananom = eccanom - orbel->e * sin(eccanom);
		//     if (orbel->meananom> M_PI) orbel->meananom-= 2.0*M_PI;
		//     if (orbel->meananom< -M_PI) orbel->meananom+= 2.0*M_PI;
		// only shift M if elliptic orbit
	} else {
		eccanom = 2.0 * atanh(foo * tan(trueanom / 2.0));
		orbel->meananom = orbel->e * sinh(eccanom) - eccanom;
	}
	if (orbel->argperi > M_PI)
		orbel->argperi -= 2.0 * M_PI;
	if (orbel->argperi < -M_PI)
		orbel->argperi += 2.0 * M_PI;
}

/* given orbital elements return PhaseState */
// hyperbolic orbits have e>1, not sure about convention for a positive or not
void kep_to_cart(double GM, OrbitalElements orb_el, PhaseState *state) {
	double meanmotion, cosE, sinE, foo;
	double x, y, z, xd, yd, zd;
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
	xd = -a * meanmotion * sinE / rovera;
	yd = foo * a * meanmotion * cosE / rovera;
	zd = 0.0;
	if (e > 1.0)
		x *= -1.0;

	/* rotate by argument of perihelion in orbit plane*/
	cosw = cos(argperi);
	sinw = sin(argperi);
	xp = x * cosw - y * sinw;
	yp = x * sinw + y * cosw;
	zp = z;
	xdp = xd * cosw - yd * sinw;
	ydp = xd * sinw + yd * cosw;
	zdp = zd;

	/* rotate by inclination about x axis */
	cosi = cos(i);
	sini = sin(i);
	x = xp;
	y = yp * cosi - zp * sini;
	z = yp * sini + zp * cosi;
	xd = xdp;
	yd = ydp * cosi - zdp * sini;
	zd = ydp * sini + zdp * cosi;

	/* rotate by longitude of node about z axis */
	cosnode = cos(longnode);
	sinnode = sin(longnode);
	state->x = x * cosnode - y * sinnode;
	state->y = x * sinnode + y * cosnode;
	state->z = z;
	state->vx = xd * cosnode - yd * sinnode;
	state->vy = xd * sinnode + yd * cosnode;
	state->vz = zd;
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
