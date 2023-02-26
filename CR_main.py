import CR_plotter as plotter
import sys
import CR_cmd as cmd


if __name__ == '__main__':

    argv = sys.argv
    argc = len(argv)

    CMD = cmd.CommandLineInput(argc, argv)
    CMD.execute()

    # plotter.plotPolytropicRhoRadius()
    # plotter.plotPolytropic()
    # plotter.plotWhiteDwarfRange()
    # plotter.plotNeutronStarRange()
    # plotter.plotLaneEmden()
    # plotter.plotPolytropicEOS()
    # plotter.plotWhiteDwarfEOS()
    # plotter.plotNeutronStarEOS()
    # plotter.plotLaneEmdenAnalyticalAll()
