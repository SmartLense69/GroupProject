import plotter

if __name__ == '__main__':

    n = 5
    plotEuler = True
    plotRK4 = True
    plotRKF = True

    plotter.plotBegin()
    plotter.plotLaneEmden(plotEuler, plotRK4, plotRKF, n)
    plotter.plotLaneEmdenAnalytical(n)
    plotter.plotEnd()
