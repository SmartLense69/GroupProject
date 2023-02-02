import numpy as np
import diffequation as de

# Raise exceptions instead of runtime warnings.
# Introduced so debugging could be easier.
# np.seterr(all='raise')


def euler2Diff(diff1: de.DifferentialEquation, diff2: de.DifferentialEquation,
               stepSize: float, stop: float, yLim: float = None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Solves two sets of decoupled differential equations (DE) from
    a 2nd order DEs.
    Example: Let the first differential euqation be x=x'(t) and the second  y=y'(x)
    :param diff1: First decoupled DE.
    :param diff2: Second decoupled DE.
    :param stepSize: Step size of the algorithm. Increase for more accurate results within the end of interval.
    :param stop: To which time the algorithm should run. For better results at the end of a solution, increase `stop`.
    :param yLim: When should the function stop if the solution is smaller than ylim.
    :return: Tuple, with t, x and y as solution
    """

    # Create arrays for time, x and y
    _tValues = np.arange(0, stop+stepSize, stepSize)
    _xValues = np.arange(0, stop+stepSize, stepSize)
    _yValues = np.arange(0, stop+stepSize, stepSize)

    # Defining initial conditions
    _xValues[0] = diff1.iniC
    _yValues[0] = diff2.iniC

    # Instead of asking if a ylim is None every loop,
    # do different loops depending on the optional parameters
    if yLim is None:
        for i in range(1, len(_tValues)):
            _xValues[i] = _xValues[i - 1] + stepSize * diff1.func(_tValues[i - 1], _yValues[i - 1], _xValues[i - 1])
            _yValues[i] = _yValues[i - 1] + stepSize * diff2.func(_xValues[i])
    else:
        for i in range(1, len(_tValues)):
            _xValues[i] = _xValues[i - 1] + stepSize * diff1.func(_tValues[i - 1], _yValues[i - 1], _xValues[i - 1])
            _yValues[i] = _yValues[i - 1] + stepSize * diff2.func(_xValues[i])
            if _yValues[i] < yLim:
                _tValues = _tValues[:i]
                _xValues = _xValues[:i]
                _yValues = _yValues[:i]
                break

    return _tValues, _xValues, _yValues


def eulerTOV(diff1: de.DifferentialEquation, diff2: de.DifferentialEquation,
               stepSize: float, stop: float, yLim: float = None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Solves two sets of decoupled differential equations (DE) from
    a 2nd order DEs.
    Example: Let the first differential euqation be x=x'(t) and the second  y=y'(x)
    :param diff1: First decoupled DE.
    :param diff2: Second decoupled DE.
    :param stepSize: Step size of the algorithm. Increase for more accurate results within the end of interval.
    :param stop: To which time the algorithm should run. For better results at the end of a solution, increase `stop`.
    :param yLim: When should the function stop if the solution is smaller than ylim.
    :return: Tuple, with t, x and y as solution
    """

    # Create arrays for time, x and y
    _tValues = np.arange(0, stop+stepSize, stepSize)
    _xValues = np.arange(0, stop+stepSize, stepSize)
    _yValues = np.arange(0, stop+stepSize, stepSize)

    # Defining initial conditions
    _xValues[0] = diff1.iniC
    _yValues[0] = diff2.iniC

    # Instead of asking if a ylim is None every loop,
    # do different loops depending on the optional parameters
    if yLim is None:
        for i in range(1, len(_tValues)):
            _xValues[i] = _xValues[i - 1] + stepSize * diff1.func(_tValues[i - 1])
            _yValues[i] = _yValues[i - 1] + stepSize * diff2.func(_tValues[i], _yValues[i - 1])
    else:
        for i in range(1, len(_tValues)):
            _xValues[i] = _xValues[i - 1] + stepSize * diff1.func(_tValues[i - 1])
            _yValues[i] = _yValues[i - 1] + stepSize * diff2.func(_tValues[i], _yValues[i - 1])
            if _yValues[i]/_yValues[0] < yLim:
                _tValues = _tValues[:i]
                _xValues = _xValues[:i]
                _yValues = _yValues[:i]
                break

    return _tValues, _xValues, _yValues


def eulerlimit(masscont, starfunc, statefunction, m0, rho0, G, c0, K, n, h, *arg):
    """
    :param masscont:
    :param starfunc:
    :param statefunc
    :param m0:
    :param rho0:
    :param G:
    :param c0:
    :param K:
    :param n:
    :param h:
    :param arg:
    :return:
    """

    # initialise the arrays to be used
    rvalues = np.arange(0, 1e8+h, h)
    # number of steps
    steps = np.shape(rvalues)[0]
    # solution arrays
    msol = np.zeros(steps)
    Psol = np.zeros(steps)

    # set the initial conditions
    msol[0] = m0
    Psol[0] = statefunction(K, n, rho=rho0)

    # Loop over radius intervals
    for i in range(1, steps):

        # find m using euler method
        msol[i] = msol[i-1] + h*masscont(rvalues[i], Psol[i-1], K, n, *arg)
        # find P using euelr method
        Psol[i] = Psol[i - 1] + h*starfunc(rvalues[i], Psol[i-1], msol[i], G, c0, K, n, *arg)

        if Psol[i]/Psol[0] < 1e-5:
            Psol = Psol[:i]
            rvalues = rvalues[:i]
            msol = msol[:i]
            break

    return (rvalues, Psol, msol)