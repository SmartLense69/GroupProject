import numpy as np
import matplotlib.pyplot as plt

def phigrad(xi, n, theta0, phi0):
    if xi==0:
        out = 0
    else:
        out = ((-2/xi)*phi0) - theta0**n
    return out

def RK4(xi1, phi0, theta0, h):

    # constants
    c = np.array([None, 0, 1/2, 1/2, 1])
    a = np.array([None, None, None, None],
                 [None, None, None, None],
                 [None, 1/2, 0, 0],
                 [None, None, 1/2, 0],
                 [None, None, None, 1])
    b = np.array([None, 1/6, 1/3, 1/3, 1/6])

    # initialise the arrays to be used

    xi = np.arange(0, xi1 + h, h)
    # number of time steps
    n = np.shape(xi)[0]

    phi = np.zeros(n)
    theta = np.zeros(n)

    # set the initial conditions
    phi[0] = phi0
    theta[0] = theta0

    # Loop over time
    for i in range(1, n):

        # find phi using RK4 method

        # calculate the values k1 through k4
        k1 = h * func(theta[i -1], xi[i - 1] + c[1]*h, phi[i - 1])
        k2 = h * func(theta[i -1], xi[i - 1] + c[2]*h, phi[i - 1] + a[2, 1]*k1)
        k3 = h * func(theta[i -1], xi[i - 1] + c[3]*h, phi[i - 1] + a[3, 1]*k1 + a[3, 2]*k2)
        k4 = h * func(theta[i -1], xi[i - 1] + c[4]*h, phi[i - 1] + a[4, 1]*k1 + a[4, 2]*k2 + a[4, 3]*k3)

        # calculate the next value for phi
        phi[i] = phi[i - 1] + b1*k1 + b2*k2 + b3*k3 + b4*k4

        # find theta using simple Euler method
        theta[i] = theta[i - 1] + h * phi[i]

    return (xi, theta)

RK4(1, 0, 1, 0.01)