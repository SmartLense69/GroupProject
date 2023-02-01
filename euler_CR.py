import numpy as np
import diffequation as de

# Raise exceptions instead of runtime warnings.
# Introduced so debugging could be easier.
# np.seterr(all='raise')


def euler2Diff(diff1: de.DifferentialEquation, diff2: de.DifferentialEquation,
                  stepSize: float, stop: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Solves two sets of decoupled differential equations (DE) from
    a 2nd order DEs.
    Example: Let the first differential euqation be x=x'(t) and the second  y=y'(x)
    :param diff1: First decoupled DE.
    :param diff2: Second decoupled DE.
    :param stepSize: Step size of the algorithm. Increase for more accurate results
    within the end of interval defined by parameter `stop`.
    :param stop: To which time the algorithm should run.
    For better results at the end of a solution, increase `stop`.
    :return: Tuple, with t, x and y as solution
    """

    # Create arrays for time, x and y
    _tValues = np.arange(0, stop+stepSize, stepSize)
    _xValues = np.arange(0, stop+stepSize, stepSize)
    _yValues = np.arange(0, stop+stepSize, stepSize)

    # Defining inital conditions
    _yValues[0] = diff2.iniC
    _xValues[0] = diff1.iniC

    # Go through the algorithm
    for i in range(1, len(_tValues)):
        _xValues[i] = _xValues[i-1] + stepSize*diff2.func(_tValues[i-1], _yValues[i-1], _xValues[i-1])
        _yValues[i] = _yValues[i-1] + stepSize*diff2.func(_xValues[i-1])

    return _tValues, _xValues, _yValues


def euler1Diff(diff1: de.DifferentialEquation, stepSize: float, stop: float) -> tuple[np.ndarray, np.ndarray]:
    """
    Solves a differential equations (DE) from a 1st order DE.
    Example: x'(t)=t
    :param diff1: DE
    :param stepSize: Step size of the algorithm. Increase for more accurate results
    within the end of interval defined by parameter `stop`.
    :param stop: To which time the algorithm should run.
    For better results at the end of a solution, increase `stop`.
    :return: Tuple, with t and x as solution
    """

    # Create arrays for time, x and y
    _tValues = np.arange(0, stop + stepSize, stepSize)
    _xValues = np.arange(0, stop + stepSize, stepSize)

    # Defining inital conditions
    _xValues[0] = diff1.iniC

    # Go through the algorithm
    for i in range(1, len(_tValues)):
        _xValues[i] = _xValues[i - 1] + stepSize * diff1.func(_tValues[i - 1])

    return _tValues, _xValues


def eulerOLD(stop, h, n, func1, ic1, func2=None, ic2=None):
    """
    Function solve the Lane-Emden equation for a given n

    :param h: step size
    :param n: polytropic index
    :param func1: second derivative given as a first derivative with another equation
    :param ic1: initial condition for the second derivative
    :param func2: first derivative given as a constant function
    :param ic2: initial condition for first derivative
    :return: array of xi values and their corresponding theta values
    """

    xivalues = np.arange(0, stop+h, h)
    steps = len(xivalues)
    thetasol = np.arange(0, stop+h, h)
    phisol = np.arange(0, stop+h, h)
    # defining inital conditions
    thetasol[0] = ic2
    phisol[0] = ic1

    for i in range(1, steps):
        if func2 is not None:
            phisol[i] = phisol[i-1] + h*func1(xivalues[i-1], n, thetasol[i-1], phisol[i-1])
        thetasol[i] = thetasol[i-1] + h*func2(phisol[i-1])

    return xivalues, thetasol