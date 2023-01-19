import numpy as np
import matplotlib.pyplot as plt

# this is so that the state functions can be varied
# to be implemented at a later point
def statefunc(P, K, n):
    return (P / K)**(n/(n+1))

def masscont(r, P, K, n):
    return 4 * np.pi * r**2 * statefunc(P, K, n)

def hydro(r, P, m, G, K, n):
    return - ( G * m * statefunc(P, K, n) / r**2 )

def tov(r, P, m, G, c, K, n):
    return - ( (G * m * statefunc(P, K, n)) / r**2 ) * \
             ( 1 + P / (statefunc(P, K, n) * c**2) ) * \
             ( 1 + (4 * np.pi * r**3 * P) / (m * c**2) ) * \
             ( 1 - (2 * G * m) / (r * c**2) )**(-1)

def rungekutta4tov(masscont, starfunc, r1, m0, P0, G, K, n, h, *arg):

    # constants
    c = np.array([None, 0, 1/2, 1/2, 1])
    a = np.array([[None, None, None, None],
                 [None, None, None, None],
                 [None, 1/2, None, None],
                 [None, 0, 1/2, None],
                 [None, 0, 0, 1]])
    b = np.array([None, 1/6, 1/3, 1/3, 1/6])

    # initialise the arrays to be used

    r = np.arange(0, r1 + h, h)
    # number of time steps
    steps = np.shape(r)[0]

    m = np.zeros(steps)
    P = np.zeros(steps)

    # set the initial conditions
    m[0] = m0
    P[0] = P0

    # Loop over time
    for i in range(1, steps):

        # find m using RK4 method
        # calculate the values k1 through k4
        k1 = h * masscont(r[i - 1] + c[1]*h, P[i - 1], K, n, *arg)
        k2 = h * masscont(r[i - 1] + c[2]*h, P[i - 1], K, n, *arg)
        k3 = h * masscont(r[i - 1] + c[3]*h, P[i - 1], K, n, *arg)
        k4 = h * masscont(r[i - 1] + c[4]*h, P[i - 1], K, n, *arg)

        # find next value for m
        m[i] = m[i - 1] + b[1]*k1 + b[2]*k2 + b[3]*k3 + b[4]*k4

        # find P using RK4 method
        # calculate the values k1 through k4
        k1 = h * starfunc(r[i - 1] + c[1] * h, P[i - 1], m[i], G, K, n, *arg)
        k2 = h * starfunc(r[i - 1] + c[2] * h, P[i - 1] + a[2, 1] * k1, m[i], G, K, n, *arg)
        k3 = h * starfunc(r[i - 1] + c[3] * h, P[i - 1] + a[3, 1] * k1 + a[3, 2] * k2, m[i], G, K, n, *arg)
        k4 = h * starfunc(r[i - 1] + c[4] * h, P[i - 1] + a[4, 1] * k1 + a[4, 2] * k2 + a[4, 3] * k3, m[i], G, K, n, *arg)

        # find next value for m
        P[i] = P[i - 1] + b[1] * k1 + b[2] * k2 + b[3] * k3 + b[4] * k4

    return (r, P)

# outputs for different values of n
n_list = [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]
colours = ['r', 'orange', 'yellow', 'lightgreen', 'g', 'cyan', 'b', 'purple', 'pink', 'k', 'gray']

# plot configs
plt.figure(figsize=(12, 8))

for n, c in zip(n_list, colours):

    # run RK4
    r, P = rungekutta4tov(phigrad, 3, 0, 1, n, 0.01)

    # plot curve
    plt.plot(r, P, c=c, label='n={0}'.format(n))

# plot configs
plt.title('Euler Lane-Emden')
plt.xlabel(r'$\xi$')
plt.ylabel(r'$\theta$')
plt.legend()
plt.show()