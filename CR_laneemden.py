import numpy as np
import CR_diffsolver as df
import CR_exceptions as ex
import CR_config as cf
from typing import Tuple

class LaneEmden:
    """
    This class handles behaviour about the Lane-Emden equation.
    """

    def _phi(self, listInput: np.ndarray) -> float:
        """
        Function for gradient phi.

        Is of form: ``((-2 / xi) * phi0) - theta0^n``

        :param listInput: [phi, theta, xi]
        :return: Gradient of phi.
        """
        phi0 = listInput[0]
        theta0 = listInput[1]
        xi = listInput[2]

        # Added xi == 0 as a failsafe.
        if xi == 0:
            return xi
        else:
            return ((-2 / xi) * phi0) - theta0 ** self.n

    @staticmethod
    def _theta(phi: np.ndarray):
        """
        Function for theta gradient.

        Is of form: ``dtheta/dx = phi``

        :param phi:
        :return:
        """
        return phi[0]

    def getSolution(self, xiH=0.0005, stopTime=cf.Sys.laneEmdenRunTime, method="rk4") -> Tuple[np.ndarray, np.ndarray]:
        """
        Gets the solution for Lane-Emden.
        :param xiH: The step size of xi.
        :param stopTime: How far xi should be running to
        :param method: The numerical method to solve it for.
        :return: xi and theta.
        """
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
        """
        Gets the analytical solution.
        :param num:
        :param stopTime:
        :return:
        """
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

    def __init__(self, inputN: float):
        self.n = inputN
