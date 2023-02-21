import abc
import CR_config as cf
import numpy as np


class Star(abc.ABC):

    density: float

    @property
    def getDensity(self):
        raise NotImplementedError

    @property
    def pressureEOS(self, *args):
        raise NotImplementedError

    @property
    def densityEOS(self, *args):
        raise NotImplementedError

    @property
    def massEquation(self, *args):
        raise NotImplementedError

    @property
    def nonRelativePressure(self, *args):
        raise NotImplementedError

    @property
    def relativePressure(self, *args):
        raise NotImplementedError


class Polytropic(Star):

    @property
    def getDensity(self):
        return self.getDensity

    def pressureEOS(self, inputList: np.ndarray):
        rho = inputList
        return cf.Var.K * (rho ** (1 + 1 / cf.Var.n))

    def densityEOS(self, inputList: np.ndarray):
        P = inputList
        return (np.abs(P) / cf.Var.K) ** (cf.Var.n / (cf.Var.n + 1))

    def massEquation(self, inputList: np.ndarray):
        P = inputList[0]
        r = inputList[1]
        return 4 * np.pi * r ** 2 * self.densityEOS(P)

    def nonRelativePressure(self, inputList: np.ndarray):
        m = inputList[0]
        P = inputList[1]
        r = inputList[2]
        return - ((cf.G.whatUnitHuh * m * self.densityEOS(P)) / r ** 2)

    def relativePressure(self, inputList: np.ndarray):
        m = inputList[0]
        P = inputList[1]
        r = inputList[2]
        return - ((cf.G.whatUnitHuh * m * self.densityEOS(P)) / r ** 2) \
            * (1 + P / (self.densityEOS(P) * cf.C.cmtrPsec ** 2)) \
            * (1 + (4 * np.pi * r ** 3 * P) / (m * cf.C.cmtrPsec ** 2)) \
            * (1 - (2 * cf.G.whatUnitHuh * m) / (r * cf.C.cmtrPsec ** 2)) ** (-1)

    def __init__(self, density):
        self.density = density


