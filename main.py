import massradius as mr
import plotter
import numpy as np

if __name__ == '__main__':
    # mr.plotMassRadius()

    n = 1
    plotEuler = True
    plotRK4 = False
    plotRKF = False

    plotter.plotBegin()
    plotter.plotLaneEmden(plotEuler, plotRK4, plotRKF, n)
    # plotter.plotLaneEmdenAnalytical(n)
    plotter.plotEnd()
