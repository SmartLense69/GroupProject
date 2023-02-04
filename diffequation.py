
class DifferentialEquation:
    """
    An object to represent a differential equation.
    Made to simplify the input of the RKF algorithm.
    """

    def __init__(self, decoupledFunction, initialCondition):
        """
        Constructor of the Differential Equation Class.
        :param decoupledFunction: The differential equation,
        represented by a pointer to a function
        :param initialCondition: The initial condition for t = 0
        """
        self.func: callable(any) = decoupledFunction
        self.iniC: float = initialCondition
