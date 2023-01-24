import numpy as np
import matplotlib.pyplot as plt

np.seterr(all='raise')


# this is so that the state functions can be varied
def statefunc(K, n, P=None, rho=None):
    if P != None:
        return np.power(P / K, np.abs(n/(n+1)))
    if rho != None:
        return K * np.power(rho,(1 + 1/n))


# mass continuity equation
def masscont(r, P, K, n):
    return 4 * np.pi * r**2 * statefunc(K, n, P=P)


# hydrstatic equilibrium equation
def hydro(r, P, m, G, c, K, n):
    return - ( G * m * statefunc(K, n, P=P) / r**2 )


# TOV equation
def tov(r, P, m, G, c, K, n):
    return - ( (G * m * statefunc(K, n, P=P)) / r**2 ) * ( 1 + P / (statefunc(K, n, P=P) * c**2) ) * ( 1 + (4 * np.pi * r**3 * P) / (m * c**2) ) * ( 1 - (2 * G * m) / (r * c**2) )**(-1)


def euler1(masscont, starfunc, r1, m0, rho0, G, c0, K, n, h, *arg):
    """

    :param masscont:
    :param starfunc:
    :param r1: maximum radius value
    :param m0:
    :param rho0:
    :param G:
    :param c0:
    :param K:
    :param n:
    :param h: stepsize
    :param arg:
    :return:
    """

    # initialise the arrays to be used

    rvalues = np.arange(0, r1 + h, h)
    # number of steps
    steps = np.shape(rvalues)[0]
    msol = np.zeros(steps)
    Psol = np.zeros(steps)

    # set the initial conditions
    msol[0] = m0
    Psol[0] = statefunc(K, n, rho=rho0)

    # Loop over radius intervals
    for i in range(1, steps):

        # find m using euler method
        msol[i] = msol[i-1] + h*masscont(rvalues[i], np.abs(Psol[i-1]), K, n, *arg)

        # find P using euler method
        Psol[i] = np.abs(Psol[i - 1]) + h*starfunc(rvalues[i], np.abs(Psol[i-1]), msol[i], G, c0, K, n, *arg)

    return (rvalues, np.abs(Psol))

# outputs for different values of n
n_list = np.asarray([0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5])
colours = ['r', 'orange', 'yellow', 'lightgreen', 'g', 'cyan', 'b', 'purple', 'pink', 'k']
# n_list = [1.5]
#colours = ['deeppink']

# plot configs
plt.figure(figsize=(12, 8))


# performing plots loop
for n, c in zip(n_list, colours):

    # run euler
    K = np.asarray([5e11])
    r, P = euler1(masscont, hydro, 1.5e8, 0, 1e9, 6.67e-8, 3e10, K, n, 1e5)

    print(len(r))

    # switch to density
    rho = (P / K)**np.abs(n/(n+1))

    # plot curve
    plt.plot(r/np.max(r), rho/np.max(rho), c=c, label='n={0}'.format(n)) #rho/np.max(rho)

# plot configs
plt.title('Euler hydrostatic')
plt.xlabel('r')
plt.ylabel(r'$\rho$')
plt.legend()
plt.show()