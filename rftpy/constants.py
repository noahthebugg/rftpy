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
        :return: gravitational constant [km^3/s^2]
        """
        return 398600.44

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
        :return: gravitational constant [km^3/s^2]
        """
        return 4900.0

"""
Returning private values as readable values
"""
Earth = _Earth()
Moon = _Moon()
