import numpy as np
import matplotlib.pyplot as plt


# this is so that the state functions can be varied
def statefunc(K, n, P=None, rho=None):
    if P is not None:
        return (np.abs(P) / K)**(n/(n+1))
    if rho != None:
        return K * rho**(1 + 1/n)


def masscont(r, P, K, n):
    return 4 * np.pi * r**2 * statefunc(K, n, P=P)


def hydro(r, P, m, G, c, K, n):
    return - ( G * m * statefunc(K, n, P=P) / r**2 )


def tov(r, P, m, G, c, K, n):
    return - ( (G * m * statefunc(K, n, P=P)) / r**2 ) * ( 1 + P / (statefunc(K, n, P=P) * c**2) ) * ( 1 + (4 * np.pi * r**3 * P) / (m * c**2) ) * ( 1 - (2 * G * m) / (r * c**2) )**(-1)


def rungekutta4(masscont, starfunc, r1, m0, rho0, G, c0, K, n, h, *arg):

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
    P[0] = statefunc(K, n, rho=rho0)

    # Loop over time
    for i in range(1, steps):

        # find m using RK4 method
        # calculate the values k1 through k4
        k1 = h * masscont(r[i] + c[1]*h, np.abs(P[i - 1]), K, n, *arg)
        k2 = h * masscont(r[i] + c[2]*h, np.abs(P[i - 1]), K, n, *arg)
        k3 = h * masscont(r[i] + c[3]*h, np.abs(P[i - 1]), K, n, *arg)
        k4 = h * masscont(r[i] + c[4]*h, np.abs(P[i - 1]), K, n, *arg)

        # find next value for m
        m[i] = m[i - 1] + b[1]*k1 + b[2]*k2 + b[3]*k3 + b[4]*k4

        # find P using RK4 method
        # calculate the values k1 through k4
        k1 = h * starfunc(r[i] + c[1] * h, np.abs(P[i - 1]), m[i], G, c0, K, n, *arg)
        k2 = h * starfunc(r[i] + c[2] * h, np.abs(P[i - 1]) + a[2, 1] * k1, m[i], G, c0, K, n, *arg)
        k3 = h * starfunc(r[i] + c[3] * h, np.abs(P[i - 1]) + a[3, 1] * k1 + a[3, 2] * k2, m[i], G, c0, K, n, *arg)
        k4 = h * starfunc(r[i] + c[4] * h, np.abs(P[i - 1]) + a[4, 1] * k1 + a[4, 2] * k2 + a[4, 3] * k3, m[i], G, c0, K, n, *arg)

        # find next value for m
        P[i] = np.abs(P[i - 1]) + b[1] * k1 + b[2] * k2 + b[3] * k3 + b[4] * k4

        if P[i]/P[0] < 1e-5:
            P = P[:i]
            r = r[:i]
            m = m[:i]
            break

    return (r, P, m)

# outputs for different values of n
n_list = [0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]
colours = ['r', 'orange', 'yellow', 'lightgreen', 'g', 'cyan', 'b', 'purple', 'pink', 'k']

# plot configs
plt.figure(figsize=(12, 8))

for n, c in zip(n_list, colours):

    # run RK4
    r, P, m = rungekutta4(masscont, tov, 2e8, 0, 1e9, 6.67e-8, 3e10, 5e11, n, 1e5)
    print('n:', n, 'R:', r[-1], 'M:', m[-1])
    print()

    # switch to density
    rho = statefunc(5e11, n, P=P)

    # plot curve
    plt.plot(r/np.max(r), rho/np.max(rho), c=c, label='n={0}'.format(n))

# plot configs
plt.title('RK4 hydrostatic')
plt.xlabel('r')
plt.ylabel(r'$\rho$')
plt.legend()
plt.show()