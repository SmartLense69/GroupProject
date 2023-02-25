"""
*Physical constants*
To prevent confusion with units, each constant is its own class.
To get the constant in the corresponding unit,
checkout the fields of each constant.

Example:
    Speed of light - C
    Speed of light in meter per seconds - C.mps
"""
import numpy as np


class Sys:

    threads = 16
    laneEmdenRunTime = 50


class C:
    """
    Speed of light.
    Available in units:

    * meter per second (MtrPSec)
    * centimeter per second (cmtrPsec)
    """

    mtrPsec = 2.99792458e8
    cmtrPsec = 2.99792458e10


class G:
    """
    Gravitational constant G.
    Available in units:

    * meter^3 per kg second^2 (Mtr3PKgSec2)

    """
    Mtr3PKgSec2 = 6.67384e-11

    # TODO: Determine units of this constant.
    whatUnitHuh = 6.67384e-8


def resetN():
    Var.n = 1.5


class Var:
    """
    Variables like K (equation of state coefficient) or n.
    Includes:

    * K
    * n
    * m0
    * rho0
    """
    K = 1e13
    n = 1.5
    m = 1
    rho = 1.5e11


class RKF:

    A = np.asarray([
        [None, None, None, None, None, None, None],
        [None, None, None, None, None, None, None],
        [None, (1 / 4), None, None, None, None, None],
        [None, (3 / 32), (9 / 32), None, None, None, None],
        [None, (1932 / 2197), (-7200 / 2197), (7296 / 2197), None, None, None],
        [None, (439 / 216), (-8), (3680 / 513), (-845 / 4104), None, None],
        [None, (-8 / 27), (2), (-3544 / 2565), (1859 / 4104), (-11 / 40), None],
    ], dtype=np.float64)
    """
    "a" Coefficients for Runge-Kutta-Fehlberg.
    Index the array like the index mentioned in the math scripts, i.e.
    first element starts at 1.
    This was made to make it easier to copy the algorithm.
    """

    B0 = np.asarray(
        [None, (35 / 384), 0, (500 / 1113), (125 / 192), (-2187 / 6784), (11 / 84), 0], dtype=np.float64)

    B = np.asarray(
        [None, (16 / 135), 0, (6656 / 12825), (28561 / 56430), (-9 / 50), (2 / 55)], dtype=np.float64)
    """
    "b" Coefficients for Runge-Kutta-Fehlberg.
    """

    Bm = np.asarray(
        [None, (5179 / 57600), 0, (7517 / 16695), (393 / 640), (-92097 / 339200), (187 / 2100), (1 / 40)],
        dtype=np.float64)
    """
    "b*" Coefficients for Runge-Kutta-Fehlberg.
    """

    C0 = np.asarray(
        [None, 0, (1 / 5), (3 / 10), (4 / 5), (8 / 9), 1, 1], dtype=np.float64)

    C = np.asarray(
        [None, 0, (1 / 4), (3 / 8), (12 / 13), 1, (1 / 2)], dtype=np.float64)
    """
    "c" Coefficients for Runge-Kutta-Fehlberg.
    """

class RK4:

    A = np.array([[None, None, None, None],
                     [None, None, None, None],
                     [None, 1 / 2, None, None],
                     [None, 0, 1 / 2, None],
                     [None, 0, 0, 1]], dtype=np.float64)
    B = np.array([None, 1 / 6, 1 / 3, 1 / 3, 1 / 6], dtype=np.float64)
    C = np.array([None, 0, 1 / 2, 1 / 2, 1], dtype=np.float64)


class RK2:

    C = np.asarray([None, 1, 0.5], dtype=np.float64)
