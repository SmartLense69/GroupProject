import numpy as np
import matplotlib.pyplot as plt


def phigrad(r, m, P, K=None, n=None):
    xi = r
    phi = m
    theta = P
    if xi == 0:
        out = 0
    else:
        out = -((2/xi)*phi + theta**n)
    return out


def thetagrad(r, m, P, G, c, K, n):
    phi = m
    return phi


# this is so that the state functions can be varied
def statefunc(K, n, P=None, rho=None):
    if P is not None:
        return (np.abs(P) / K)**(n/(n+1))
    if rho != None:
        return K * rho**(1 + 1/n)


def masscont(r, m, P=None, K=None, n=None):
    return 4 * np.pi * r**2 * statefunc(K, n, P=P)


def hydro(r, m, P, G, c, K, n):
    return - ( G * m * statefunc(K, n, P=P) / r**2 )


def tov(r, m, P, G, c, K, n):
    return - ( (G * m * statefunc(K, n, P=P)) / r**2 ) * \
        ( 1 + P / (statefunc(K, n, P=P) * c**2) ) * \
        ( 1 + (4 * np.pi * r**3 * P) / (m * c**2) ) * \
        ( 1 - (2 * G * m) / (r * c**2) )**(-1)


def rungekutta4(diffEq1, diffEq2, statefunc, t1, x0, y0, G, c0, K, n, h, *arg):

    # constants
    #c = np.array([None, 0, 1/2, 1/2, 1])
    #a = np.array([[None, None, None, None],
    #             [None, None, None, None],
    #             [None, 1/2, None, None],
    #             [None, 0, 1/2, None],
    #             [None, 0, 0, 1]])
    #b = np.array([None, 1/6, 1/3, 1/3, 1/6])

    # initialise the arrays to be used

    t = np.arange(1, t1 + h, h)
    # number of time steps
    steps = np.shape(t)[0]

    x = np.zeros(steps)
    y = np.zeros(steps)

    # set the initial conditions
    x[0] = x0
    if diffEq2 == thetagrad:
        y[0] = y0
    else:
        y[0] = statefunc(K, n, rho=y0)

    # Loop over time
    for i in range(1, steps):

        # midpoint method
        j1 = h * diffEq1(t[i - 1], x[i - 1], y[i - 1], K=K, n=n, *arg)
        k1 = h * diffEq2(t[i - 1], x[i - 1], y[i - 1], G, c0, K, n, *arg)

        j2 = h * diffEq1(t[i - 1] + 0.5 * h, x[i - 1] + 0.5 * j1, y[i - 1] + 0.5 * k1, K=K, n=n, *arg)
        k2 = h * diffEq2(t[i - 1] + 0.5 * h, x[i - 1] + 0.5 * j1, y[i - 1] + 0.5 * k1, G, c0, K, n, *arg)

        x[i] = x[i - 1] + j2
        y[i] = y[i - 1] + k2

        # # find m using RK4 method
        # # calculate the values k1, j1 through k4, j4
        # j1 = h * diffEq1(t[i - 1], x[i - 1], y[i - 1], K=K, n=n, *arg)
        # k1 = h * diffEq2(t[i - 1], x[i - 1], y[i - 1], G, c0, K, n, *arg)
        #
        # j2 = h * diffEq1(t[i - 1] + 0.5 * h, x[i - 1] + 0.5 * j1, y[i - 1] + 0.5 * k1, K=K, n=n, *arg)
        # k2 = h * diffEq2(t[i - 1] + 0.5 * h, x[i - 1] + 0.5 * j1, y[i - 1] + 0.5 * k1, G, c0, K, n, *arg)
        #
        # j3 = h * diffEq1(t[i - 1] + 0.5 * h, x[i - 1] + 0.5 * j2, y[i - 1] + 0.5 * k2, K=K, n=n, *arg)
        # k3 = h * diffEq2(t[i - 1] + 0.5 * h, x[i - 1] + 0.5 * j2, y[i - 1] + 0.5 * k2, G, c0, K, n, *arg)
        #
        # j4 = h * diffEq1(t[i - 1] + h, x[i - 1] + j3, y[i - 1] + k3, K=K, n=n, *arg)
        # k4 = h * diffEq2(t[i - 1] + h, x[i - 1] + j3, y[i - 1] + k3, G, c0, K, n, *arg)
        #
        # # find next value for m
        # x[i] = x[i - 1] + (j1 + 2 * j2 + 2 * j3 + j4) / 6
        #
        # # find next value for P
        # y[i] = y[i - 1] + (k1 + 2 * k2 + 2 * k3 + k4) / 6

        if y[i]/y[0] < 1e-5:
            t = t[:i]
            x = x[:i]
            y = y[:i]
            break

    return t, x, y

# outputs for different values of n
n_list = [0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]
colours = ['r', 'orange', 'yellow', 'lightgreen', 'g', 'cyan', 'b', 'purple', 'pink', 'k']

# plot configs
plt.figure(figsize=(12, 8))

for n, c in zip(n_list, colours):

    # run RK4
    #xi, phi, theta = rungekutta4(phigrad, thetagrad, statefunc, 35, 0, 1, 6.67e-8, 3e10, 1e13, n, 0.01) #6.67e-8, 3e10, 1e13, n, 0.01)

    # plot curve
    #plt.plot(xi, theta, c=c, label='n={0}'.format(n))

    # run RK4
    r, m, P = rungekutta4(masscont, tov, statefunc, 1e11, 1, 1.5e11, 6.67e-8, 3e10, 1e13, n, 1e4)
    print('n:', n, 'R:', r[-1], 'M:', m[-1])

    # switch to density
    rho = statefunc(1.5e11, n, P=P)

    # plot curve
    plt.plot(r / np.max(r), rho / np.max(rho), c=c, label='n={0}'.format(n))

# plot configs
#plt.hlines(0, 0, 35, color='k', linestyles='--')
plt.title('RK4 hydrostatic/tov')
plt.xlabel('r')
plt.ylabel(r'$\rho$')
#plt.ylim([-0.25, 1.1])
plt.legend()
plt.show()


"""
def rungekutta4(masscont, starfunc, statefunc, r1, m0, rho0, G, c0, K, n, h, *arg):

    # constants
    c = np.array([None, 0, 1/2, 1/2, 1])
    a = np.array([[None, None, None, None],
                 [None, None, None, None],
                 [None, 1/2, None, None],
                 [None, 0, 1/2, None],
                 [None, 0, 0, 1]])
    b = np.array([None, 1/6, 1/3, 1/3, 1/6])

    # initialise the arrays to be used

    r = np.arange(1, r1 + h, h)
    # number of time steps
    steps = np.shape(r)[0]

    m = np.zeros(steps)
    P = np.zeros(steps)

    # set the initial conditions
    m[0] = m0
    if starfunc != thetagrad:
        P[0] = statefunc(K, n, rho=rho0)
    else:
        P[0] = rho0

    # Loop over time
    for i in range(1, steps):

        # find m using RK4 method
        # calculate the values k1, j1 through k4, j4
        j1 = h * masscont(r[i - 1] + c[1] * h, np.abs(P[i - 1]), m[i - 1], K, n, *arg)
        k1 = h * starfunc(r[i - 1] + c[1] * h, np.abs(P[i - 1]), m[i - 1], G, c0, K, n, *arg)

        j2 = h * masscont(r[i - 1] + c[2] * h, np.abs(P[i - 1]) + a[2, 1] * k1, m[i - 1] + a[2, 1] * j1, K, n, *arg)
        k2 = h * starfunc(r[i - 1] + c[2] * h, np.abs(P[i - 1]) + a[2, 1] * k1, m[i - 1] + a[2, 1] * j1, G, c0, K, n, *arg)

        j3 = h * masscont(r[i - 1] + c[3] * h, np.abs(P[i - 1]) + a[3, 1] * k1 + a[3, 2] * k2,
                          m[i - 1] + a[3, 1] * j1 + a[3, 2] * j2, K, n, *arg)
        k3 = h * starfunc(r[i - 1] + c[3] * h, np.abs(P[i - 1]) + a[3, 1] * k1 + a[3, 2] * k2,
                          m[i - 1] + a[3, 1] * j1 + a[3, 2] * j2, G, c0, K, n, *arg)

        j4 = h * masscont(r[i - 1] + c[4] * h, np.abs(P[i - 1]) + a[4, 1] * k1 + a[4, 2] * k2 + a[4, 3] * k3,
                          m[i - 1] + a[4, 1] * j1 + a[4, 2] * j2 + a[4, 3] * j3, K, n, *arg)
        k4 = h * starfunc(r[i - 1] + c[4] * h, np.abs(P[i - 1]) + a[4, 1] * k1 + a[4, 2] * k2 + a[4, 3] * k3,
                          m[i - 1] + a[4, 1] * j1 + a[4, 2] * j2 + a[4, 3] * j3, G, c0, K, n, *arg)

        # find next value for m
        m[i] = m[i - 1] + b[1] * j1 + b[2] * j2 + b[3] * j3 + b[4] * j4

        # find next value for P
        P[i] = np.abs(P[i - 1]) + b[1] * k1 + b[2] * k2 + b[3] * k3 + b[4] * k4

        if P[i]/P[0] < 1e-5:
            P = P[:i]
            r = r[:i]
            m = m[:i]
            break

    return r, P, m



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

# run RK4
    xi, theta, phi = rungekutta4(phigrad, thetagrad, statefunc, 20, 0, 1, 6.67e-8, 3e10, 1e13, n, 0.01)
    #print('n:', n, 'R:', r[-1], 'M:', m[-1])

    # switch to density
    #rho = statefunc(5e11, n, P=P)

    # plot curve
    plt.plot(xi, theta, c=c, label='n={0}'.format(n))

# run RK4
    r, P, m = rungekutta4(masscont, tov, statefunc, 1e10, 1, 1e9, 6.67e-8, 3e10, 1e13, n, 1e4)
    print('n:', n, 'R:', r[-1], 'M:', m[-1])

    # switch to density
    rho = statefunc(5e11, n, P=P)

    # plot curve
    plt.plot(r / np.max(r), rho / np.max(rho), c=c, label='n={0}'.format(n))
"""