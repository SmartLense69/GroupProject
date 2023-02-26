"""
    Describes behaviour about stars.
    Includes polytropic and neutron stars, as well as white dwarfs
"""

# abc stands for abstract base classes
# This enables abstract classes and derivations

from typing import Tuple
import abc
import CR_config as cf
import numpy as np
import CR_diffsolver as df
import CR_exceptions as ex

# TODO: [Marc] Add string parsing for CR_star objects.


class Star(abc.ABC):
    """
        A star is defined by an equation of state (EOS), which should return density and pressure separately.
        The class has the relativistic and the non-relativistic pressure function included, (i.e. TOV and Hydro)
        as well as the mass equation, which is same for all types of stars.
        A function to get the mass-radius and the density-radius is provided for all child classes.
    """

    density: float

    @property
    def pressureEOS(self, *args):
        """
            Equation of state, that should return a pressure.
            Overwrite to the equation of state for the child class.
            :param args: inputList, preferably involving density.
            :return: Pressure.
        """

        raise NotImplementedError

    @property
    def densityEOS(self, *args):
        """
            Equation of state, that should return a density.
            Overwrite to the equation of state for the child class.
            :param args: inputList, preferably involving pressure.
            :return: Density.
        """

        raise NotImplementedError

    def massEquation(self, inputList: np.ndarray):
        """
            The mass equation of the star, being ``4 * np.pi * r ** 2 * self.densityEOS(P)``
            :param inputList: [pressure, radius]
            :return:
        """
        P = inputList[0]
        r = inputList[1]
        return 4 * np.pi * r ** 2 * self.densityEOS(P)

    def nonRelativePressure(self, inputList: np.ndarray):
        """
        Non-relativistic pressure equation, also known as the hydrostatic equation.
        :param inputList: [mass, pressure, radius]
        :return: pressure
        """
        m = inputList[0]
        P = inputList[1]
        r = inputList[2]
        return - ((cf.G.whatUnitHuh * m * self.densityEOS(P)) / r ** 2)

    def relativisticPressure(self, inputList: np.ndarray):
        """
            Relativistic pressure equation, also known as the Tolman–Oppenheimer (TOV) equation.
            :param inputList: [mass, pressure, radius]
            :return: pressure
        """

        m = inputList[0]
        P = inputList[1]
        r = inputList[2]
        return - ((cf.G.whatUnitHuh * m * self.densityEOS(P)) / r ** 2) \
            * (1 + P / (self.densityEOS(P) * (cf.C.cmtrPsec ** 2))) \
            * (1 + (4 * np.pi * (r ** 3) * P) / (m * (cf.C.cmtrPsec ** 2))) \
            * ((1 - (2 * cf.G.whatUnitHuh * m) / (r * (cf.C.cmtrPsec ** 2))) ** (-1))

    def getDensityRadius(self, rhoH: float = 1e4, stopTime: float = 2e10, method: str = "rk4",
                         pressure: str = "Relativistic", verbose: bool = True) -> Tuple[np.ndarray, np.ndarray]:
        """
        Get the density-radius data for the star.
        :param rhoH: Density step. Do not increase blindly.
        :param stopTime: How far should the differential equation solver solve for, in terms of density.
        :param method: The numerical method that should be used to solve the neutron star solution.
        :param pressure: The type of pressure function. Choose TOV or Hydro with ["Relativistic", "Non-Relativistic"]
        :param verbose: If additional information should be printed, like current calculation interation.
        :return: Radius and density numpy array. Returns 0, 0 if an invalid value was encountered in the calculations.
        """

        if verbose:
            print("Called with n = {0}".format(cf.Var.n))

        diffM = df.DifferentialEquation(["p", "r"], "m", self.massEquation, cf.Var.m, 1)

        if pressure == "Relativistic":
            diffP = df.DifferentialEquation(["m", "p", "r"], "p", self.relativisticPressure,
                                            self.pressureEOS(self.density), 2)
        elif pressure == "Non-Relativistic":
            diffP = df.DifferentialEquation(["m", "p", "r"], "p", self.nonRelativePressure,
                                            self.pressureEOS(self.density), 2)
        else:
            raise ex.InvalidPressureMethod(pressure)

        diffEqs = df.DifferentialEquationSystem([diffM, diffP])
        diffS = df.DifferentialSolver(diffEqs, rhoH, stopTime=stopTime)

        rMod = diffS.varDict.get("r")
        rMod[0] = 1
        diffS.varDict.update({"r": rMod})
        diffS.addThreshold({"p": 1e-10})

        if method == "rk2":
            diffS.rk2()
        elif method == "rk4":
            diffS.rk4()
        elif method == "rkf":
            diffS.rkf()
        elif method == "euler":
            diffS.euler()
        else:
            raise ex.InvalidNumericalMethod(method)

        if len(diffS.varDict.get("r")) != 0 and len(diffS.varDict.get("m")) != 0:

            r = diffS.varDict.get("r")/1e5
            rho = self.densityEOS(diffS.varDict.get("p"))

        else:
            r = 0
            rho = 0

        return r, rho

    def getMassRadius(self, rhoH: float = 1e5, stopTime: float = 2e10, method: str = "rk4",
                      pressure: str = "Relativistic", verbose: bool = False):
        """
        Get the mass-radius data for the star.
        :param rhoH: Density step. Do not increase blindly.
        :param stopTime: How far should the differential equation solver solve for, in terms of density.
        :param method: The numerical method that should be used to solve the neutron star solution.
        :param pressure: The type of pressure function. Choose TOV or Hydro with ["Relativistic", "Non-Relativistic"]
        :param verbose: If additional information should be printed, like current calculation interation.
        :return: Radius and mass numpy array. Returns 0, 0 if an invalid value was encountered in the calculations.
        """

        dataValues = np.zeros(2)

        diffM = df.DifferentialEquation(["p", "r"], "m", self.massEquation, cf.Var.m, 1)

        if pressure == "Relativistic":
            diffP = df.DifferentialEquation(["m", "p", "r"], "p", self.relativisticPressure,
                                            self.pressureEOS(self.density), 2)
        elif pressure == "Non-Relativistic":
            diffP = df.DifferentialEquation(["m", "p", "r"], "p", self.nonRelativePressure,
                                            self.pressureEOS(self.density), 2)
        else:
            raise ex.InvalidPressureMethod(pressure)

        diffEqs = df.DifferentialEquationSystem([diffM, diffP])
        diffS = df.DifferentialSolver(diffEqs, rhoH, stopTime=stopTime)

        rMod = diffS.varDict.get("r")
        rMod[0] = 1
        diffS.varDict.update({"r": rMod})
        diffS.addThreshold({"p": 1e-5})

        if method == "rk4":
            diffS.rk4()
        elif method == "rk2":
            diffS.rk2()
        elif method == "euler":
            diffS.euler()
        elif method == "rkf":
            diffS.rkf()
        else:
            raise ex.InvalidNumericalMethod(wrongMethod=method)

        if len(diffS.varDict.get("r")) != 0 and len(diffS.varDict.get("m")) != 0:

            r = diffS.varDict.get("r")[-1] / 1e5
            m = diffS.varDict.get("m")[-1] / 2e33

            print("{0}: Calculation with rho = {1}\n"
                  "gave rise to r = {2} and m = {3}\n".format(method, self.density, r, m))

            dataValues[0] = r
            dataValues[1] = m

        else:
            print("Calculation at rho = {0} skipped!".format(self.density))
            return 0, 0

        if verbose:
            print("For {0} pressure and rho = {1}: r = {2}, m = {3}"
                  .format(pressure, self.density, dataValues[0], dataValues[1]))

        return dataValues[0], dataValues[1]

    def getEOSData(self, rhoMin=1, rhoMax=1e14, rhoNum=100):
        """
        Get equation of state data.
        :param rhoMin: The minium range for density.
        :param rhoMax: The maximum range for density.
        :param rhoNum: The number of points in that range. (Generated by np.linspace)
        :return:
        """
        densityValues = np.linspace(rhoMin, rhoMax, rhoNum)
        return densityValues, self.pressureEOS(densityValues)


class Polytropic(Star):
    """
        A polytropic star is a star, where its expansion is merely limited by gravitational and hydrostatic force.
    """

    def pressureEOS(self, inputList: np.ndarray | float):
        rho = inputList
        return cf.Var.K * (rho ** (1 + (1 / cf.Var.n)))

    def densityEOS(self, inputList: np.ndarray | float):
        P = inputList
        return (np.abs(P) / cf.Var.K) ** (cf.Var.n / (cf.Var.n + 1))

    def __init__(self, density):
        self.density = density


class WhiteDwarf(Star):
    """
        # TODO: [Olivia] Was it Fermi-Gas that was different to a polytropic star?
        A white dwarf is just a really dense polytropic star.
    """
    def pressureEOS(self, inputList: np.ndarray | float):
        return self.pressureEOS(inputList)

    def densityEOS(self, inputList: np.ndarray | float):
        return self.densityEOS(inputList)

    def __init__(self, density, cubicSplinePressure, cubicSplineDensity):
        self.density = density
        self.pressureEOS = cubicSplinePressure
        self.densityEOS = cubicSplineDensity


class NeutronStar(Star):
    """
        Neutron stars are even denser.
        Ich bin vielleicht dicht, aber Goethe is Dichter.
    """

    def pressureEOS(self, inputList: np.ndarray | float):
        return self.pressureEOS(inputList)

    def densityEOS(self, inputList: np.ndarray | float):
        return self.densityEOS(inputList)

    def __init__(self, density, cubicSplinePressure, cubicSplineDensity):
        self.density = density
        self.pressureEOS = cubicSplinePressure
        self.densityEOS = cubicSplineDensity
