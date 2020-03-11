#ifndef _KEPCART_
#define _KEPCART_

/* Phase space and orbital element structures */

// Orbital state in phase space (position, velocity)
typedef struct PhaseState {
	double x, y, z, vx, vy, vz;
} PhaseState;

// Standard Keplerian elements:
// Semimajor axis, eccentricity, inclination, longitude of ascending node, argument of periapsis, mean anomaly
typedef struct OrbitalElements {
	double a, e, i, longnode, argperi, meananom;
} OrbitalElements;

/* Physical constants
 (may want to move them) */

#define GG 6.6732e-8      // gravitational constant in cgs (cm^3/g/s^2)
#define CC 2.997924562e10 // speed of light, cm/s
#define Msol 1.989e33     // mass of sun,  g
#define AU 1.49597892e13  // astronomical unit, cm

/* Defined in kepcart.c */

/* Keplerian to Cartesian conversions */

// Cartesian to Keplerian
void cart_to_kep(double GM, PhaseState state, OrbitalElements *orbel);
// Keplerian to Cartesian
void kep_to_cart(double GM, OrbitalElements orbel, PhaseState *state);
// Helper for conversion -- Get eccentric anomaly from mean anomaly
double eccentric_anomaly(double eccentricity, double mean_anomaly); // solves kepler's eqn
// ... in hyperbolic case
double eccentric_anomaly_hyperbolic(double eccentricity, double mean_anomaly); // solves kepler's eqn

#endif
