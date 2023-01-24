import numpy as np
import matplotlib.pyplot as plt

def phigrad(xi, phi, n, theta):
    if xi==0:
        out = 0
    else:
        out = ((-2/xi)*phi) - theta**n
    return out

def rungekutta4lane(func, xi1, phi0, theta0, n, h, *arg):

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
    steps = np.shape(xi)[0]

    phi = np.zeros(steps)
    theta = np.zeros(steps)

    # set the initial conditions
    phi[0] = phi0
    theta[0] = theta0

    # Loop over time
    for i in range(1, steps):

        # find phi using RK4 method

        # calculate the values k1 through k4
        k1 = h * func(xi[i] + c[1]*h, phi[i - 1], n, theta[i -1], *arg)
        k2 = h * func(xi[i] + c[2]*h, phi[i - 1] + a[2, 1]*k1, n, theta[i -1], *arg)
        k3 = h * func(xi[i] + c[3]*h, phi[i - 1] + a[3, 1]*k1 + a[3, 2]*k2, n, theta[i -1], *arg)
        k4 = h * func(xi[i] + c[4]*h, phi[i - 1] + a[4, 1]*k1 + a[4, 2]*k2 + a[4, 3]*k3, n, theta[i -1], *arg)

        # calculate the next value for phi
        phi[i] = phi[i - 1] + b[1]*k1 + b[2]*k2 + b[3]*k3 + b[4]*k4

        # find theta using simple Euler method
        theta[i] = theta[i - 1] + h * phi[i]

    return (xi, theta)

# outputs for different values of n
n_list = [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]
colours = ['r', 'orange', 'yellow', 'lightgreen', 'g', 'cyan', 'b', 'purple', 'pink', 'k', 'gray']


# # plot configs
# plt.figure(figsize=(12, 8))
#
# for n, c in zip(n_list, colours):
#
#     # run RK algorithm
#     xi, theta = rungekutta4lane(phigrad, 2.7, 0, 1, n, 0.01)
#
#     # plot curve
#     plt.plot(xi, theta**n, c=c, label='n={0}'.format(n))
#
# # plot configs
# plt.title('RK4 Lane-Emden')
# plt.xlabel(r'$\xi$')
# plt.ylabel(r'$\theta$')
# plt.legend()
# plt.show()
