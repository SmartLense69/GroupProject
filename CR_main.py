import sys
import CR_cmd as cmd
import CR_plotter as plotter

if __name__ == '__main__':

    argv = sys.argv
    argc = len(argv)

    CMD = cmd.InputReader(argc, argv)
    CMD.execute()

    # plotter.plotTOVWhiteDwarfNeutronStar()

    # plotter.plotPolytropicRhoRadius()
    # plotter.plotPolytropic()
    # plotter.plotWhiteDwarfRange()
    # plotter.plotNeutronStarRange()
    # plotter.plotLaneEmden()
    # plotter.plotPolytropicEOS()
    # plotter.plotWhiteDwarfEOS()
    # plotter.plotNeutronStarEOS()
    # plotter.plotLaneEmdenAnalyticalAll()
