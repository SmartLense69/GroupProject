import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimize

# testing the rk4 algorithm by solving a very simple 2nd order ode: y'' = -y

# coupled equations:
# z' = -y, solve first
# y' = z, solve second

def gradz(in_y):
    return -in_y

def grady(in_z):
    return in_z

def rungekutta4(h, stop):

    # constants
    c = np.array([None, 0, 0.5, 0.5, 1])
    a = np.array([[None, None, None, None],
                 [None, None, None, None],
                 [None, 0.5, None, None],
                 [None, 0, 0.5, None],
                 [None, 0, 0, 1]])
    b = np.array([None, 1/6, 1/3, 1/3, 1/6])

    x = np.arange(0, stop+h, h)
    steps = len(x)
    y = np.zeros(steps)
    z = np.zeros(steps)

    # setting initial conditions
    y[0] = 0
    z[0] = 1

    for i in range(1, steps):
        # calculate RK4 coefficients
        j1 = h*gradz(y[i-1])
        k1 = h*grady(z[i-1])

        j2 = h*gradz(y[i-1] + k1*0.5)
        k2 = h*grady(z[i-1] + j1*0.5)

        j3 = h*gradz(y[i-1] + k2*0.5)
        k3 = h*grady(z[i-1] + j2*0.5)

        j4 = h*gradz(y[i-1] + k3)
        k4 = h*grady(z[i-1] + j3)

        z[i] = z[i-1] + (1/6)*(j1 + 2*j2 + 2*j3 + j4)
        y[i] = y[i-1] + (1/6)*(k1 + 2*k2 + 2*k3 + k4)

    return x, y

out = rungekutta4(0.0001, 10)

# analytical solution
def an(in_x):
    return np.sin(in_x)


plt.plot(out[0], out[1], label='num')
plt.plot(out[0], an(out[0]), label='ana')
plt.legend()
plt.show()

testh = np.array([5e-1, 1e-1, 5e-2, 1e-2, 5e-3, 1e-3, 5e-4])
errors = np.zeros(len(testh))
count = np.arange(0, len(testh), 1)
for i in count:
    # setting h for the iteration of the loop
    h = testh[i]
    outit = rungekutta4(h, 100)
    abserror = np.abs(an(outit[0][len(outit[0])-1]) - outit[1][len(outit[0])-1])
    errors[i] = abserror

def fitlinear(x, a, b):
    return a*x + b

fith_smooth = np.linspace(min(testh), max(testh), 1000)

params, params_covariance = optimize.curve_fit(fitlinear, np.log10(testh), np.log10(errors), p0=[1, 1])
errorsfit = fitlinear(np.log10(fith_smooth),params[0], params[1])
print("Error \u221d h^{0:.3g}".format(params[0]))

plt.scatter(testh, errors, label='Errors')
plt.plot(fith_smooth, 10**errorsfit)
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.xlabel('h')
plt.ylabel('Absolute Error')
plt.show()

