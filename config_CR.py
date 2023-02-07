"""
*Physical constants*
To prevent confusion with units, each constant is its own class.
To get the constant in the corresponding unit,
checkout the fields of each constant.

Example:
    Speed of light - C
    Speed of light in meter per seconds - C.mps
"""


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
