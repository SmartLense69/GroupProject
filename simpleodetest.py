import numpy as np
import matplotlib.pyplot as plt

# testing the rk4 algorithm by solving a very simple 2nd order ode: y'' = -y

# coupled equations:
# z' = -y, solve first
# y' = z, solve second

def zgrad(in_y):
    return -in_y

def ygrad(in_z):
    return in_z

def rungekutta4(h, stop):

    # constants
    c = np.array([None, 0, 1/2, 1/2, 1])
    a = np.array([[None, None, None, None],
                 [None, None, None, None],
                 [None, 1/2, None, None],
                 [None, 0, 1/2, None],
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
        j1 = h*zgrad(y[i-1])
        k1 = h*ygrad(z[i-1])

        j2 = h*zgrad(y[i-1] + k1*a[2,1])
        k2 = h*ygrad(z[i-1] + j1*a[2,1])

        j3 = h*zgrad(y[i-1] + k1*a[3,1] + k2*a[3,2])
        k3 = h*ygrad(z[i-1] + j1*a[3,1] + j2*a[3,2])

        j4 = h*zgrad(y[i-1] + k1*a[4,1] + k2*a[4,2] + k3*a[4,3])
        k4 = h*ygrad(z[i] + j1*a[4,1] + j2*a[4,2] + j3*a[4,3])

        z[i] = z[i-1] + b[1]*j1 + b[2]*j2 + b[3]*j3 + b[4]*j4
        y[i] = y[i-1] + b[1]*k1 + b[2]*k2 + b[3]*k3 + b[4]*k4

    return x, y

out = rungekutta4(0.0001, 10)

# analytical solution
def an(in_x):
    return np.sin(in_x)

plt.plot(out[0], out[1], label='num')
plt.plot(out[0], an(out[0]), label='ana')
plt.legend()
plt.show()

testh = [1, 5e-1, 1e-1, 5e-2, 1e-2, 5e-3, 1e-3, 5e-4, 1e-4]
errors = np.zeros(len(testh))
count = np.arange(0, len(testh), 1)
for i in count:
    # setting h for the iteration of the loop
    h = testh[i]
    outit = rungekutta4(h, 10)
    print(outit[1][len(outit[0])-1])
    print(an(outit[0][len(outit[0])-1]))
    abserror = np.abs( an(outit[0][len(outit[0])-1]) - outit[1][len(outit[0])-1] )
    errors[i] = abserror


# plt.scatter(testh, errors)
# plt.legend()
# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel('h')
# plt.ylabel('Absolute Error')
# plt.show()

