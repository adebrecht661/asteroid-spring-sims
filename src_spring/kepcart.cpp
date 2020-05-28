/*
 * Keplerian conversion routines
 */

#include <cmath>
#include "matrix_math.h"
#include "kepcart.h"

using std::abs;

// Returns orbital elements given the phase space state
OrbitalElements cart_to_kep(double GM, PhaseState state) {

	// Get angular momentum
	Vector pos = state.x;
	Vector r_hat = pos / pos.len();
	Vector vel = state.v;
	Vector h = cross(pos, vel);

	// Init return values
	OrbitalElements orb_el;

	// Calculate semimajor axis
	// Negative if orbit is hyperbolic
	orb_el.a = 1.0 / (2.0 / pos.len() - pow(vel.len(), 2.0) / GM);

	// Calculate eccentricity
	double p = pow(h.len(), 2.) / GM;
	orb_el.e = sqrt(1. - p / orb_el.a);

	// Calculate orbital inclination
	orb_el.i = acos(h.getZ() / h.len());

	// Calculate longitude of ascending node
	if (h.getX() != 0.0 || h.getY() != 0.0) {
		orb_el.long_asc_node = atan2(h.getX(), -h.getY());
	} else {
		orb_el.long_asc_node = 0.0;
	}

	// Calculate eccentricity
	double e_cos_true_anom = pow(h.len(), 2.0) / (GM * pos.len()) - 1.0;
	double e_sin_true_anom = dot(r_hat, vel) * h.len() / GM;

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

	// Calculate mean anomaly using series expansion
	orb_el.mean_anom = true_anom - 2 * orb_el.e * sin(true_anom)
			+ (0.75 * pow(orb_el.e, 2.) + 0.125 * pow(orb_el.e, 4.))
					* sin(2 * true_anom)
			- 1. / 3. * pow(orb_el.e, 3.) * sin(3 * true_anom)
			+ 5. / 32. * pow(orb_el.e, 4.) * sin(4 * true_anom);

	return orb_el;
}

// Given orbital elements, return PhaseState
// Hyperbolic orbits have e > 1
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

	PhaseState state;

	// Rotate about z axis by argument of periapsis, x by inclination, z by longitude of ascending node
	Matrix rot_peri = getRotMatZ(orb_el.arg_peri);
	Matrix rot_inc = getRotMatX(orb_el.i);
	Matrix rot_long = getRotMatZ(orb_el.long_asc_node);

	state.x = rot_long * rot_inc * rot_peri * pos;
	state.v = rot_long * rot_inc * rot_peri * vel;

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

	if (ecc < 0) {
		throw "Eccentricity cannot be negative.";
	}
	if (ecc > 1) {
		throw "You're in a hyperbolic orbit. Use eccentric_anomaly_hyperbolic.";
	}

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

	if (ecc < 0) {
		throw "Eccentricity cannot be negative.";
	}
	if (ecc < 1) {
		throw "You're in an elliptic orbit. Use eccentric_anomaly.";
	}

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

/*************/
/* Operators */
/*************/

// Equality
bool operator==(OrbitalElements lhs, OrbitalElements rhs) {
	return lhs.a == rhs.a && lhs.e == rhs.e && lhs.i == rhs.i
			&& lhs.long_asc_node == rhs.long_asc_node
			&& lhs.arg_peri == rhs.arg_peri
			&& lhs.mean_anom == rhs.long_asc_node;
}

// Stream output
std::ostream& operator<<(std::ostream &os, const OrbitalElements &orb_el) {
	os << "Semi-major axis: " << orb_el.a << "\nEccentricity: " << orb_el.e
			<< "\nInclination: " << orb_el.i
			<< "\nLongitude of ascending node: " << orb_el.long_asc_node
			<< "\nArgument of periapsis: " << orb_el.arg_peri
			<< "\nMean anomaly: " << orb_el.mean_anom << std::endl;
	return os;
}

// Equality
bool operator==(PhaseState lhs, PhaseState rhs) {
	return lhs.x == rhs.x && lhs.v == rhs.v;
}

// Stream output
std::ostream& operator<<(std::ostream &os, const PhaseState &state) {
	os << "Position: " << state.x << "\nVelocity: " << state.v << std::endl;
	return os;
}
