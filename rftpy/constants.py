"""Constants for space calculations."""

def constant(f):
    def fset(self, value):
        """
        :param self: constant
        :param value: value to be set
        :return: TypeError
        """
        raise TypeError
    def fget(self):
        """
        :param self: constant
        :return: value of constant
        """
        return f(self)
    return property(fget, fset)


class _Space(object):
    """
    General space constants
    """
    @constant
    def G(self) -> float:
        """
        :return: gravitational constant [m^3/kgs^-2]
        """
        return 6.67430e-11


class _Earth(object):
    """
    Earth constants
    """
    @constant
    def radius(self) -> float:
        """
        :return: radius [km]
        """
        return 6371.0
    @constant
    def mu(self) -> float:
        """
        :return: standard gravitational parameter [km^3/s^2]
        """
        return 398600.4415
    @constant
    def g0(self) -> float:
        """
        :return: [m/s^2]
        """
        return 9.81


class _Moon(object):
    """
    Moon constants
    """
    @constant
    def radius(self) -> float:
        """
        :return: radius [km]
        """
        return 1738.1
    @constant
    def mu(self) -> float:
        """
        :return: standard gravitational parameter [km^3/s^2]
        """
        return 4900.0


# Returning private values as readable values
Space = _Space()
Earth = _Earth()
Moon = _Moon()
