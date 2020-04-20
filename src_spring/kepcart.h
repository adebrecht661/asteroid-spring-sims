#ifndef _KEPCART_
#define _KEPCART_

/*
 * Phase space and orbital element structures
 */

// Orbital state in phase space (position, velocity)
typedef struct PhaseState {
	double x, y, z, v_x, v_y, v_z;
} PhaseState;

// Standard Keplerian elements:
// Semimajor axis, eccentricity, inclination, longitude of ascending node, argument of periapsis, mean anomaly
typedef struct OrbitalElements {
	double a, e, i, long_asc_node, arg_peri, mean_anom;
} OrbitalElements;

/**************************************/
/* Keplerian to Cartesian conversions */
/**************************************/

// Cartesian to Keplerian
void cart_to_kep(double GM, PhaseState state, OrbitalElements *orbel);
// Keplerian to Cartesian
void kep_to_cart(double GM, OrbitalElements orbel, PhaseState *state);
// Helper for conversion -- Get eccentric anomaly from mean anomaly
double eccentric_anomaly(double eccentricity, double mean_anomaly); // solves kepler's eqn
// ... in hyperbolic case
double eccentric_anomaly_hyperbolic(double eccentricity, double mean_anomaly); // solves kepler's eqn

#endif
