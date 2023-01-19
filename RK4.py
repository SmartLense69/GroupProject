import numpy as np
import matplotlib.pyplot as plt

def phigrad(xi, phi, n, theta):
    if xi==0:
        out = 0
    else:
        out = ((-2/xi)*phi) - theta**n
    return out

def RK4(func, xi1, phi0, theta0, n, h, *arg):

    # constants
    c = np.array([None, 0, 1/2, 1/2, 1])
    a = np.array([[None, None, None, None],
                 [None, None, None, None],
                 [None, 1/2, None, None],
                 [None, 0, 1/2, None],
                 [None, 0, 0, 1]])
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
        k1 = h * func(xi[i - 1] + c[1]*h, phi[i - 1], theta[i -1], n, *arg)
        k2 = h * func(xi[i - 1] + c[2]*h, phi[i - 1] + a[2, 1]*k1, theta[i -1], n, *arg)
        k3 = h * func(xi[i - 1] + c[3]*h, phi[i - 1] + a[3, 1]*k1 + a[3, 2]*k2, theta[i -1], n, *arg)
        k4 = h * func(xi[i - 1] + c[4]*h, phi[i - 1] + a[4, 1]*k1 + a[4, 2]*k2 + a[4, 3]*k3, theta[i -1], n, *arg)

        # calculate the next value for phi
        phi[i] = phi[i - 1] + b[1]*k1 + b[2]*k2 + b[3]*k3 + b[4]*k4

        # find theta using simple Euler method
        theta[i] = theta[i - 1] + h * phi[i]

    return (xi, theta)

# outputs for different values of n
xi0, theta0 = RK4(phigrad, 1, 0, 1, 0, 0.01)
xi1, theta1 = RK4(phigrad, 1, 0, 1, 1, 0.01)
xi5, theta5 = RK4(phigrad, 1, 0, 1, 5, 0.01)

# plot
plt.plot(xi0, theta0, label='n=0')
plt.plot(xi1, theta1, label='n=1')
plt.plot(xi5, theta5, label='n=5')
plt.title('Euler Lane-Emden')
plt.xlabel(r'$\xi$')
plt.ylabel(r'$\theta$')
plt.legend()
plt.show()
