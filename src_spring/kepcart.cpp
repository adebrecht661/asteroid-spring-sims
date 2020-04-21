/*
 * Keplerian conversion routines
 */

#include <cmath>
#include "matrix_math.h"
#include "kepcart.h"

// Returns orbital elements given the phase space state
OrbitalElements cart_to_kep(double GM, PhaseState state) {

	// Get angular momentum
	Vector pos = state.x;
	Vector r_hat = pos / pos.len();
	Vector vel = state.v;
	Vector L = cross(pos, vel);

	// Init return values
	OrbitalElements orb_el;
	// Calculate orbital inclination
	orb_el.i = acos(L.getZ() / L.len());

	// Calculate longitude of ascending node
	if (L.getX() != 0.0 || L.getY() != 0.0) {
		orb_el.long_asc_node = atan2(L.getX(), -L.getY());
	} else {
		orb_el.long_asc_node = 0.0;
	}

	// Calculate semimajor axis
	// Negative if orbit is hyperbolic
	orb_el.a = 1.0 / (2.0 / pos.len() - pow(vel.len(), 2.0) / GM);

	// Calculate eccentricity
	double e_cos_true_anom = pow(L.len(), 2.0) / (GM * pos.len()) - 1.0;
	double e_sin_true_anom = dot(r_hat, vel) * L.len() / GM;
	orb_el.e = sqrt(
			e_cos_true_anom * e_cos_true_anom
					+ e_sin_true_anom * e_sin_true_anom);

	// Calculate true anomaly
	double true_anom;
	if (e_sin_true_anom != 0.0 || e_cos_true_anom != 0.0) {
		true_anom = atan2(e_sin_true_anom, e_cos_true_anom);
	} else {
		true_anom = 0.0;
	}

	// Calculate argument of latitude
	double cos_node = cos(orb_el.long_asc_node);
	double sin_node = sin(orb_el.long_asc_node);

	double rcosu = pos.getX() * cos_node + pos.getY() * sin_node;
	if (cos(orb_el.i) == 0.0)
		throw "You're about to divide by zero.";
	double rsinu = (pos.getY() * cos_node - pos.getX() * sin_node)
			/ cos(orb_el.i);

	double u;
	if (rsinu != 0.0 || rcosu != 0.0) {
		u = atan2(rsinu, rcosu);
	} else {
		u = 0.0;
	}

	// Calculate and shift argument of periapsis
	orb_el.arg_peri = u - true_anom;
	if (orb_el.arg_peri > M_PI)
		orb_el.arg_peri -= 2.0 * M_PI;
	if (orb_el.arg_peri < -M_PI)
		orb_el.arg_peri += 2.0 * M_PI;

	// Calculate mean anomaly (and eccentric anomaly)
	double foo = sqrt(abs(1.0 - orb_el.e) / (1.0 + orb_el.e));
	double ecc_anom;
	if (orb_el.e < 1.0) {
		ecc_anom = 2.0 * atan(foo * tan(true_anom / 2.0));
		orb_el.mean_anom = ecc_anom - orb_el.e * sin(ecc_anom);
	} else {
		ecc_anom = 2.0 * atanh(foo * tan(true_anom / 2.0));
		orb_el.mean_anom = orb_el.e * sinh(ecc_anom) - ecc_anom;
	}

	return orb_el;
}

/* given orbital elements return PhaseState */
// hyperbolic orbits have e>1, not sure about convention for a positive or not
PhaseState kep_to_cart(double GM, OrbitalElements orb_el) {

	// Calculate eccentric anomaly
	double e = orb_el.e, ecc_anom;
	if (e < 1) {
		ecc_anom = eccentric_anomaly(e, orb_el.mean_anom);
	} else {
		ecc_anom = eccentric_anomaly_hyperbolic(e, orb_el.mean_anom);
	}

	// Calculate unrotated position and velocity
	double cos_ecc, sin_ecc;
	if (e < 1.0) {
		cos_ecc = cos(ecc_anom);
		sin_ecc = sin(ecc_anom);
	} else {
		cos_ecc = cosh(ecc_anom);
		sin_ecc = sinh(ecc_anom);
	}
	double a = abs(orb_el.a);
	double mean_motion = sqrt(GM / pow(a, 3.0));
	double foo = sqrt(abs(1.0 - pow(e, 2.0)));
	double r_over_a = (1.0 - e * cos_ecc);
	if (e > 1.0)
		r_over_a *= -1.0;
	double x = a * (cos_ecc - e);
	if (e > 1.0)
		x *= -1.0;
	Vector pos = { x, foo * a * sin_ecc, 0.0 };
	Vector vel = { -a * mean_motion * sin_ecc / r_over_a, foo * a * mean_motion
			* cos_ecc / r_over_a, 0.0 };

	// Rotate about z axis by argument of periapsis
	Matrix rot_peri = getRotMatZ(orb_el.arg_peri);
	pos = rot_peri * pos;
	vel = rot_peri * vel;

	// Rotate about x axis by inclination
	Matrix rot_inc = getRotMatX(orb_el.i);
	pos = rot_inc * pos;
	vel = rot_inc * vel;

	// Rotate about z axis by longitude of ascending node
	Matrix rot_long = getRotMatZ(orb_el.long_asc_node);
	pos = rot_long * pos;
	vel = rot_long * vel;

	// Set and return state vector
	PhaseState state;
	state.x = pos;
	state.v = vel;

	return state;
}

/************/
/* Helpers  */
/************/

// Solver precision
#define PREC_ecc_ano 1e-16

// Solve Kepler's equation
// Iterate to get an estimate for eccentric anomaly (u) given the mean anomaly (l).
// Works for values in all quadrants and at high eccentricity
double eccentric_anomaly(double ecc, double mean_anom) {
	// Initialize change, guess for eccentric anomaly
	double delta_u = 1.0;
	double u_guess = mean_anom + ecc * sin(mean_anom)
			+ 0.5 * pow(ecc, 2) * sin(2.0 * mean_anom);

	// Iterate to find eccentric anomaly
	// From Brouwer+Clemence, also see M+D equations 2.55, 2.58
	// Good to second order in e
	int counter = 0;
	while (abs(delta_u) > PREC_ecc_ano && counter < 10000) {

		// Current estimate for mean anomaly
		double l_guess = u_guess - ecc * sin(u_guess);

		// Adjust eccentric anomaly based on current error in mean anomaly
		delta_u = (mean_anom - l_guess) / (1.0 - ecc * cos(u_guess));

		// Update guess at eccentric anomaly
		u_guess += delta_u;

		// Increment counter
		counter++;
	}

	// Return best guess for eccentric anomaly
	return u_guess;
}

// Solve Kepler's equation in hyperbolic case
// Iterate to get an estimate for eccentric anomaly (u) given the mean anomaly (l)
double eccentric_anomaly_hyperbolic(double ecc, double mean_anom) {
	// Initialize change, guess for eccentric anomaly
	double delta_u = 1.0;
	double u_guess = log(2.0 * mean_anom / ecc + 1.8); // Danby guess

	// Iterate to find eccentric anomaly
	int counter = 0;
	while (abs(delta_u) > PREC_ecc_ano && counter < 10000) {

		// Current guess for mean anomaly
		double l_guess = u_guess - ecc * sinh(u_guess);

		// Adjust eccentric anomaly based on current error in mean anomaly
		delta_u = (-mean_anom - l_guess) / (1.0 - ecc * cosh(u_guess));

		// Update guess at eccentric anomaly
		u_guess += delta_u;

		// Increment counter
		counter++;
	}

	// Return best guess for eccentric anomaly
	return u_guess;
}
