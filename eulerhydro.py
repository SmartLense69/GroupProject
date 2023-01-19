import numpy as np
import matplotlib.pyplot as plt
import constants as consts

def rhograd(phi0):
    out = phi0
    return out

def phigrad(r, n, rho0, phi0, G, K):
    out = ((-((1/n)-1)/(rho0))*phi0**2)-((2/r)*rho0)-((rho0**(2/(1/n)))*G*4*np.pi / K*(1+(1/n)))
    return out

def euler1(h, n, func1, ic1, func2=None, ic2=None):
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

    rvalues = np.arange(0, 7e6+h, h)
    steps = len(rvalues)
    rhosol = np.arange(0, 7e6+h, h)
    phisol = np.arange(0, 7e6+h, h)
    # defining inital conditions
    rhosol[0] = ic2
    phisol[0] = ic1

    for i in range(1, steps):
        if func2 is not None:
            phisol[i] = phisol[i-1] + h*func1(rvalues[i-1], n, rhosol[i-1], phisol[i-1], consts.G, 1)
        rhosol[i] = rhosol[i-1] + h*func2(phisol[i-1])

    return rvalues, rhosol


result = euler1(10, 1, phigrad, 0, rhograd, 1e9)
plt.plot(result[0], result[1])
plt.xlabel('r')
plt.ylabel(r'$\rho$')
plt.show()