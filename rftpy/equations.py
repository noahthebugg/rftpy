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
    c = linear_eccentricity(a, b)
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
    :return: radius at perigee [km]
    """
    return a * (1 - e)

def r_apogee(a: float, e: float) -> float:
    """
    Calculate the radius at apogee.
    :param a: semi-major axis [km]
    :param e: numerical eccentricity [1]
    :return: radius at apogee [km]
    """
    return a * (1 + e)

def orbit_interval(a: float, mu: float) -> float:
    """
    Calculate the orbit interval.
    :param a: semi-major axis [km]
    :param mu: gravitational constant [km^3/s^2]
    :return: orbit interval [s]
    """
    return 2 * math.pi * math.sqrt(a**3 / mu)

def interval2a(interval: float, mu: float) -> float:
    """
    Calculate the semi-major axis by orbit interval.
    :param interval: orbit interval [s]
    :param mu: gravitational constant [km^3/s^2]
    :return: semi-major axis [km]
    """
    return ((interval / (2 * math.pi)) ** 2 * mu) ** (1/3)

def vis_viva(mu: float, r: float, a: float) -> float:
    """
    Vis-Viva-Equation: Calculate the orbital velocity.
    :param mu: gravitational constant [km^3/s^2]
    :param r: orbital radius [km]
    :param a: semi-major axis [km]
    :return: orbital velocity [km/s]
    """
    return math.sqrt(mu * (2/r - 1/a))

def red_vis_viva(mu: float, r: float) -> float:
    """
    Reduced Vis-Viva-Equation: Calculate the orbital velocity.
    :param mu: gravitational constant [km^3/s^2]
    :param r: orbital radius [km]
    :return: orbital velocity [km/s]
    """
    return math.sqrt(mu / r)

def v_cosmic1(mu: float, r: float) -> float:
    """
    Calculate the 1st cosmic velocity (Reduced Vis-Viva-Equation).
    :param mu: gravitational constant [km^3/s^2]
    :param r: orbital radius [km]
    :return: v_c (1st cosmic velocity) [km/s]
    """
    return red_vis_viva(mu, r)

def v_cosmic2(mu: float, r: float) -> float:
    """
    Calculate the 2nd cosmic velocity.
    :param mu: gravitational constant [km^3/s^2]
    :param r: orbital radius [km]
    :return: v_esc (2nd cosmic velocity) [km/s]
    """
    return math.sqrt(2 * (mu / r))

def angular_momentum(r: float, v: float) -> float:
    """
    Calculate the angular momentum.
    :param r: orbital radius [km]
    :param v: orbital velocity [km/s]
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