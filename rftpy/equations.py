"""Equations for space calculations."""

import math


def linear_eccentricity(a: float, b: float) -> float:
    """
    Calculate the linear eccentricity coefficient.
    :param a: semi-major axis [km]
    :param b: semi-minor axis [km]
    :return: linear eccentricity [km]
    """
    return math.sqrt(a**2 - b**2)


def numerical_eccentricity(a: float, b: float) -> float:
    """
    Calculate the numerical eccentricity coefficient.
    :param a: semi-major axis [km]
    :param b: semi-minor axis [km]
    :return: numerical eccentricity [1]
    """
    c: float = linear_eccentricity(a, b)
    return c / a


def semi_latus_rectum(a: float, b: float) -> float:
    """
    Calculate the Semi-Latus Rectum.
    :param a: semi-major axis [km]
    :param b: semi-minor axis [km]
    :return: Semi-Latus Rectum [km]
    """
    return b**2 / a


def r_nu(p: float, e: float, v: float) -> float:
    """
    Calculate the radius of curvature.
    :param p: Semi-Latus Rectum [km]
    :param e: numerical eccentricity [1]
    :param v: angle [rad]
    :return: radius of curvature (r_nu) [km]
    """
    return p / (1 + e * math.cos(v))


def r_perigee(a: float, e: float) -> float:
    """
    Calculate the radius at perigee.
    :param a: semi-major axis [km]
    :param e: numerical eccentricity [1]
    :return: distance at perigee [km]
    """
    return a * (1 - e)


def r_apogee(a: float, e: float) -> float:
    """
    Calculate the radius at apogee.
    :param a: semi-major axis [km]
    :param e: numerical eccentricity [1]
    :return: distance at apogee [km]
    """
    return a * (1 + e)


def orbit_interval(a: float, mu: float) -> float:
    """
    Calculate the interval of an orbiting body.
    :param a: semi-major axis [km]
    :param mu: standard gravitational parameter [km^3/s^2]
    :return: interval of an orbiting body [s]
    """
    return 2 * math.pi * math.sqrt(a**3 / mu)


def interval2a(interval: float, mu: float) -> float:
    """
    Calculate the semi-major axis by the orbit interval of an orbiting body.
    :param interval: interval of an orbiting body [s]
    :param mu: standard gravitational parameter [km^3/s^2]
    :return: semi-major axis [km]
    """
    return ((interval / (2 * math.pi)) ** 2 * mu) ** (1/3)


def velocity2a(v: float, mu: float, r: float) -> float:
    """
    Calculate the semi-major axis by the velocity of an orbiting body.
    :param v: velocity of an orbiting body [km/s]
    :param mu: standard gravitational parameter [km^3/s^2]
    :param r: distance of an orbiting body [km]
    :return: semi-major axis [km]
    """
    return -1 / (v**2 / mu - (2/r))


def vis_viva(mu: float, r: float, a: float) -> float:
    """
    Vis-Viva-Equation: Calculate the speed of an orbiting body.
    :param mu: standard gravitational parameter [km^3/s^2]
    :param r: distance of an orbiting body [km]
    :param a: semi-major axis [km]
    :return: speed of an orbiting body [km/s]
    """
    return math.sqrt(mu * (2/r - 1/a))


def red_vis_viva(mu: float, r: float) -> float:
    """
    Reduced Vis-Viva-Equation: Calculate the speed of an orbiting body.
    :param mu: standard gravitational parameter [km^3/s^2]
    :param r: distance of an orbiting body [km]
    :return: speed of an orbiting body [km/s]
    """
    return math.sqrt(mu / r)


def v_cosmic1(mu: float, r: float) -> float:
    """
    Calculate the 1st cosmic velocity (Reduced Vis-Viva-Equation).
    :param mu: standard gravitational parameter [km^3/s^2]
    :param r: distance of an orbiting body [km]
    :return: v_c (1st cosmic velocity) [km/s]
    """
    return red_vis_viva(mu, r)


def v_cosmic2(mu: float, r: float) -> float:
    """
    Calculate the 2nd cosmic velocity.
    :param mu: standard gravitational parameter [km^3/s^2]
    :param r: distance of an orbiting body [km]
    :return: v_esc (2nd cosmic velocity) [km/s]
    """
    return math.sqrt(2 * (mu / r))


def angular_momentum(r: float, v: float) -> float:
    """
    Calculate the angular momentum.
    :param r: distance of an orbiting body [km]
    :param v: speed of an orbiting body [km/s]
    :return: orbital angular momentum [km^2/s]
    """
    return r * v


def spc_angular_momentum(h: float, m: float) -> float:
    """
    Calculate the specific angular momentum of a satellite.
    :param h: orbital angular momentum [km^2/s]
    :param m: satellite mass [kg]
    :return: specific angular momentum [km^2/kgs]
    """
    return h / m


def drive_req(v1: float, v2: float, di:float) -> float:
    """
    Calculate the required drive velocity.
    :param v1: speed of an orbiting body [km/s]
    :param v2: speed of an orbiting body [km/s]
    :param di: change of inclination (i2 - i1) [deg]
    :return: change of velocity [km/s]
    """
    return math.sqrt(v1**2 + v2**2 - (2 * v1 * v2 * math.cos(math.radians(di))))


def n(mu: float, a: float) -> float:
    """
    Calculate the average motion.
    :param mu: standard gravitational parameter [km^3/s^2]
    :param a: semi-major axis [km]
    :return: n (average motion) [1/s]
    """
    return math.sqrt(mu / a**3)


def domega(n: float, j: float, r_e: float, a: float, i: float, e: float) -> float:
    """
    Calculate the node drift.
    :param n: n (average motion) [1/s]
    :param j: geopotential-coefficient [1]
    :param r_e: planet radius [km]
    :param a: semi-major axis [km]
    :param i: iclination [deg]
    :param e: numerical eccentricity [1]
    :return: domega (node drift) [1/s]
    """
    return - (3/2) * n * j * (r_e / a)**2 * (math.cos(math.radians(i)) / (1 - e**2)**2)


def ziolkowski(c: float, m1: float, m2: float) -> float:
    """
    Ziolkowski rocket equation.
    :param c: Exit speed [km/s]
    :param m1: Starting mass [kg]
    :param m2: End mass [kg]
    :return: change of velocity [km/s]
    """
    return c * math.log((m1 / m2), math.e)


def mass_flow(m: float, dt: float) -> float:
    """
    Calculate the mass flow.
    :param m: mass [kg]
    :param dt: time passed [s]
    :return: mass flow [kg/s]
    """
    return m / dt


def c_e(F: float, dotm: float) -> float:
    """
    Calculate the nozzle exit speed.
    :param F: boost power [N]
    :param dotm: mass flow [kg/s]
    :return: exit speed [m/s]
    """
    return F / dotm


def boost(c_e: float, dotm: float) -> float:
    """
    Calculate the rocket boost.
    :param c_e: exit speed [m/s]
    :param dotm: mass flow [kg/s]
    :return: boost power [N]
    """
    return c_e * dotm


def grav_force(M: float, m: float, r: float) -> float:
    """
    Newtons law: Calculate the gravitational force.
    :param M: planet mass [kg]
    :param m: satellite mass [kg]
    :param r: radius [km]
    :return: gravitational force [N]
    """
    G: float = 6.67430e-11  # Gravitational constant [m^3/kgs^-2]
    return -G * (M * m) / r ** 2


def di(v1: float, dv: float,) -> float:
    """
    Calculate the change of inclination.
    :param v1: velocity of the circular orbit [km/s]
    :param dv: change of velocity [km/s]
    :return: change of inclination [deg]
    """
    return math.asin(dv / (2 * v1)) * 2


def v_0_lat(r: float, lat: float) -> float:
    """
    Calculate the starting velocity depending on the latitude.
    :param r: planet radius [km]
    :param lat: latitude [deg]
    :return: starting velocity [km/s]
    """
    return ((2 * math.pi * r) / (24 * 3600)) * math.cos(math.radians(lat))


def alpha(r2_p: float, v2_p: float, r1: float, v2: float) -> float:
    """
    Calculate the angle between an circular and elliptical orbit.
    :param r2_p: radius at perigee [km]
    :param v2_p: velocity at perigee [km/s]
    :param r1: radius of circular orbit [km]
    :param v2: velocity on circular orbit [km/s]
    :return: angle [deg]
    """
    return math.degrees(math.acos(
        (r2_p * v2_p) / (r1 * v2)
    ))
