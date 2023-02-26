import numpy as np
import CR_diffsolver as df
import CR_exceptions as ex
import CR_config as cf


class LaneEmden:
    """
        This class handles behaviour about the Lane-Emden equation.
        It includes
    """

    def _phi(self, listInput: np.ndarray):
        phi0 = listInput[0]
        theta0 = listInput[1]
        xi = listInput[2]
        if xi == 0:
            return xi
        else:
            return ((-2 / xi) * phi0) - theta0 ** self.n

    def _theta(self, phi: np.ndarray):
        return phi[0]

    def getSolution(self, xiH=0.0005, stopTime=cf.Sys.laneEmdenRunTime, method="rk4"):
        diffEq1 = df.DifferentialEquation(["phi", "theta", "xi"], "phi", self._phi, 0, 2)
        diffEq2 = df.DifferentialEquation(["phi"], "theta", self._theta, 1, 0)
        differentialSystem = df.DifferentialEquationSystem([diffEq1, diffEq2])
        differentialSolver = df.DifferentialSolver(differentialSystem, xiH, stopTime, stepBegin=0)
        differentialSolver.addThreshold({"theta": -0.5})

        if method == "rk4":
            differentialSolver.rk4()
        elif method == "rk2":
            differentialSolver.rk2()
        elif method == "rkf":
            differentialSolver.rkf()
        elif method == "euler":
            differentialSolver.euler()
        else:
            raise ex.InvalidNumericalMethod(method)

        return differentialSolver.varDict.get("xi"), differentialSolver.varDict.get("theta")

    def getAnalyticalSolution(self, num=1000, stopTime=cf.Sys.laneEmdenRunTime):
        match self.n:
            case 0.0:
                xi = np.linspace(0, stopTime, num)
                val = -(1 / 6) * (xi ** 2) + 1
                return xi, val
            case 1.0:
                xi = np.linspace(0, stopTime, num)
                val = np.sin(xi) / xi
                return xi, val
            case 5.0:
                xi = np.linspace(0, stopTime, num)
                val = 1 / (np.sqrt(1 + ((xi ** 2) / 3)))
                return xi, val
            case _:
                raise ex.NoAnalyticalSolution

    def __init__(self, inputN: int):
        self.n = inputN
