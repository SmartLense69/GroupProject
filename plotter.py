import matplotlib.pyplot as plt
import numpy as np

import rungekuttafehlberg as rkf
import euler as eu
import rungekutta4lane as rk4l


def plotBegin():
    plt.rcParams['font.size'] = '14'
    plt.title("Dimension-less density $\\theta$ over\ndimension-less radius $\\xi$")


def plotLaneEmden(euPlot, rk4Plot, rkfPlot, n):
    if euPlot:
        eu.plot(n)
    if rk4Plot:
        rk4l.plot(n)
    if rkfPlot:
        rkf.plot(n)


def plotLaneEmdenAnalytical(n, num=1000):
    match n:
        case 0:
            xi = np.linspace(0, 1, num)
            val = -(1 / 6) * (xi ** 2) + 1
            plt.plot(xi, val, label="Analytical solution for n=0")
        case 1:
            xi = np.linspace(0, 1, num)
            val = np.sin(xi) / xi
            plt.plot(xi, val, label="Analytical solution for n=1")
        case 5:
            xi = np.linspace(0, 1, num)
            val = 1 / (np.sqrt(1 + ((xi ** 2) / 3)))
            plt.plot(xi, val, label="Analytical solution for n=5")


def plotEnd(xMin=0, xMax=1):
    #plt.xlim(xMin, xMax)
    plt.grid()
    plt.legend()
    plt.show()
