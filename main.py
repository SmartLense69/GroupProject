import plotter
import numpy as np

if __name__ == '__main__':

    n = 5
    plotEuler = True
    plotRK4 = True
    plotRKF = True

    plotter.plotBegin()
#    for n in np.arange(0, 5.5, 0.5):
    plotter.plotLaneEmden(plotEuler, plotRK4, plotRKF, n)
    plotter.plotLaneEmdenAnalytical(n)
    plotter.plotEnd()
