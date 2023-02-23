import abc

from matplotlib import pyplot as plt
from scipy.interpolate import CubicSpline

import CR_config as cf
import numpy as np
import CR_diffsolver as df
import CR_exceptions as ex


class Star(abc.ABC):
    density: float

    @property
    def pressureEOS(self, *args):
        raise NotImplementedError

    @property
    def densityEOS(self, *args):
        raise NotImplementedError

    def massEquation(self, inputList: np.ndarray):
        P = inputList[0]
        r = inputList[1]
        return 4 * np.pi * r ** 2 * self.densityEOS(P)

    def nonRelativePressure(self, inputList: np.ndarray):
        m = inputList[0]
        P = inputList[1]
        r = inputList[2]
        return - ((cf.G.whatUnitHuh * m * self.densityEOS(P)) / r ** 2)

    def relativisticPressure(self, inputList: np.ndarray):
        m = inputList[0]
        P = inputList[1]
        r = inputList[2]
        return - ((cf.G.whatUnitHuh * m * self.densityEOS(P)) / r ** 2) \
            * (1 + P / (self.densityEOS(P) * (cf.C.cmtrPsec ** 2))) \
            * (1 + (4 * np.pi * (r ** 3) * P) / (m * (cf.C.cmtrPsec ** 2))) \
            * ((1 - (2 * cf.G.whatUnitHuh * m) / (r * (cf.C.cmtrPsec ** 2))) ** (-1))

    def getDensityRadius(self, rhoH=1e4, stopTime=2e10, method="rk4", pressure="Relativistic"):

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

        if method == "rk4":
            diffS.rk4()
        elif method == "euler":
            diffS.euler()
        else:
            raise ex.InvalidNumericalMethod(method)

        if len(diffS.varDict.get("r")) != 0 and len(diffS.varDict.get("m")) != 0:

            r = diffS.varDict.get("r")/1e5
            rho = self.densityEOS(diffS.varDict.get("p"))

        else:
            return 0, 0

        return r, rho

    def getMassRadius(self, rhoH=1e5, stopTime=2e10, method="rk4", pressure="Relativistic"):

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
        elif method == "euler":
            diffS.euler()
        else:
            raise ex.InvalidNumericalMethod(method)

        if len(diffS.varDict.get("r")) != 0 and len(diffS.varDict.get("m")) != 0:

            r = diffS.varDict.get("r")[-1] / 1e5
            m = diffS.varDict.get("m")[-1] / 2e33

            print("{0}: Calculation with rho = {1}\n"
                  "gave rise to r = {2} and m = {3}\n".format(method, self.density, r, m))

            dataValues[0] = r
            dataValues[1] = m

        else:
            print("Calculation at rho = {0} skipped!".format(self.density))

        return dataValues[0], dataValues[1]

    def getEOSData(self, rhoMin=1, rhoMax=1e14, rhoNum=100):
        densityValues = np.linspace(rhoMin, rhoMax, rhoNum)
        return densityValues, self.pressureEOS(densityValues)


class Polytropic(Star):

    def pressureEOS(self, inputList: np.ndarray | float):
        rho = inputList
        return cf.Var.K * (rho ** (1 + (1 / cf.Var.n)))

    def densityEOS(self, inputList: np.ndarray | float):
        P = inputList
        return (np.abs(P) / cf.Var.K) ** (cf.Var.n / (cf.Var.n + 1))

    def __init__(self, density):
        self.density = density


class WhiteDwarf(Star):

    def pressureEOS(self, inputList: np.ndarray | float):
        return self.pressureEOS(inputList)

    def densityEOS(self, inputList: np.ndarray | float):
        return self.densityEOS(inputList)

    def __init__(self, density, cubicSplinePressure, cubicSplineDensity):
        self.density = density
        self.pressureEOS = cubicSplinePressure
        self.densityEOS = cubicSplineDensity


class NeutronStar(Star):

    def pressureEOS(self, inputList: np.ndarray | float):
        return self.pressureEOS(inputList)

    def densityEOS(self, inputList: np.ndarray | float):
        return self.densityEOS(inputList)


    def __init__(self, density, cubicSplinePressure, cubicSplineDensity):
        self.density = density
        self.pressureEOS = cubicSplinePressure
        self.densityEOS = cubicSplineDensity


