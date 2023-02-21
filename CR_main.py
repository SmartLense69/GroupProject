import sys
import numpy as np
import CR_exitcodes as ec
from CR_param import Param as P
from CR_input import isStringPartOf, subtractStringBack

argv: list[str] = sys.argv
argc: int = len(argv)
densities: list[float] = []
densitySteps: list[float] | np.ndarray = [-1, -1, -1]
method: str | None = None


def matchMissingDensityInfo(i: int):
    match i:
        case 0:
            print("Missing minimum density.")
            sys.exit(ec.MISSING_ARGUMENT)
        case 1:
            print("Missing maximum density.")
            sys.exit(ec.MISSING_ARGUMENT)
        case 2:
            print("Missing step density")
            sys.exit(ec.MISSING_ARGUMENT)


def checkDensityMinimum(arg: str):
    if isStringPartOf(arg, P.DENSITY_MIN[0]):
        strRho = subtractStringBack(arg, P.DENSITY_MIN[0])
        try:
            rho = float(strRho)
            if rho >= 0:
                if densitySteps[0] != -1:
                    densitySteps[0] = rho
                else:
                    print("Please mention the minimum density once.")
                    sys.exit(ec.INVALID_COMMAND_ARG)
            else:
                print("{0} is not a valid density. Densities must be real positive."
                      .format(density))
                sys.exit(ec.INVALID_ARGUMENT)
        except ValueError:
            print("Density '{0}' must be a number.".format(strRho))
            sys.exit(ec.INVALID_ARGUMENT)
    elif isStringPartOf(arg, P.DENSITY_MIN[1]):
        strRho = subtractStringBack(arg, P.DENSITY_MIN[1])
        try:
            rho = float(strRho)
            if rho >= 0:
                if densitySteps[0] != -1:
                    densitySteps[0] = rho
                else:
                    print("Please mention the minimum density once.")
                    sys.exit(ec.INVALID_COMMAND_ARG)
            else:
                print("{0} is not a valid density. Densities must be real positive."
                      .format(density))
                sys.exit(ec.INVALID_ARGUMENT)
        except ValueError:
            print("Minimum density '{0}' must be a number.".format(strRho))
            sys.exit(ec.INVALID_ARGUMENT)


def checkDensityStep(arg: str):
    if isStringPartOf(arg, P.DENSITY_STEP[0]):
        strRho = subtractStringBack(arg, P.DENSITY_STEP[0])
        try:
            rho = float(strRho)
            if rho >= 0:
                if densitySteps[2] != -1:
                    densitySteps[2] = rho
                else:
                    print("Please mention the step density once.")
                    sys.exit(ec.INVALID_COMMAND_ARG)
            else:
                print("{0} is not a valid density step. Densities must be real positive."
                      .format(density))
                sys.exit(ec.INVALID_ARGUMENT)
        except ValueError:
            print("Density '{0}' must be a number.".format(strRho))
            sys.exit(ec.INVALID_ARGUMENT)
    elif isStringPartOf(arg, P.DENSITY_MIN[1]):
        strRho = subtractStringBack(arg, P.DENSITY_MIN[1])
        try:
            rho = float(strRho)
            if rho >= 0:
                if densitySteps[2] != -1:
                    densitySteps[2] = rho
                else:
                    print("Please mention the step density once.")
                    sys.exit(ec.INVALID_COMMAND_ARG)
            else:
                print("{0} is not a valid density step. Densities must be real positive."
                      .format(density))
                sys.exit(ec.INVALID_ARGUMENT)
        except ValueError:
            print("Density step '{0}' must be a number.".format(strRho))
            sys.exit(ec.INVALID_ARGUMENT)


def checkDensityMaximum(arg: str) -> bool:
    if isStringPartOf(arg, P.DENSITY_MAX[0]):
        strRho = subtractStringBack(arg, P.DENSITY_MAX[0])
        try:
            rho = float(strRho)
            if rho >= 0:
                if densitySteps[1] != -1:
                    densitySteps[1] = rho
                else:
                    print("Please mention the max density once.")
                    sys.exit(ec.INVALID_COMMAND_ARG)
            else:
                print("{0} is not a valid max density. Densities must be real positive."
                      .format(density))
                sys.exit(ec.INVALID_ARGUMENT)
        except ValueError:
            print("Density '{0}' must be a number.".format(strRho))
            sys.exit(ec.INVALID_ARGUMENT)
    elif isStringPartOf(arg, P.DENSITY_MIN[1]):
        strRho = subtractStringBack(arg, P.DENSITY_MIN[1])
        try:
            rho = float(strRho)
            if rho >= 0:
                if densitySteps[2] != -1:
                    densitySteps[2] = rho
                    return True
                else:
                    print("Please mention the max density once.")
                    sys.exit(ec.INVALID_COMMAND_ARG)
            else:
                print("{0} is not a valid max density. Densities must be real positive."
                      .format(density))
                sys.exit(ec.INVALID_ARGUMENT)
        except ValueError:
            print("Density step '{0}' must be a number.".format(strRho))
            sys.exit(ec.INVALID_ARGUMENT)
    else:
        return False


if __name__ == '__main__':

    if argc - 1 == 0:

        print("No command line arguments have been given.\n"
              "Use '{0} -h' or '{0} --help' for more information."
              .format(sys.argv[0]))
        sys.exit(ec.NO_ARGUMENTS)

    else:

        if argv[1] == P.LANE_EMDEN[0] or \
                argv[1] == P.LANE_EMDEN[1]:

            if argc - 2 >= 1:

                if argv[2] == P.DENSITY[0] or \
                        argv[2] == P.DENSITY[1]:

                    if argc - 3 >= 1:
                        for index in range(3, argc):
                            try:
                                density = float(argv[index])
                                if density > 0:

                                    densities.append(density)

                                else:
                                    print("{0} is not a valid density. Densities must be real positive."
                                          .format(density))
                                    sys.exit(ec.INVALID_ARGUMENT)

                            except ValueError:
                                if argv[index] == P.RK4[0] or argv[index] == P.RK4[1]:
                                    if method is not None:
                                        print("Please provide the numerical method once.")
                                        sys.exit(ec.INVALID_COMMAND_ARG)
                                    if index == argc - 1:
                                        method = "rk4"
                                        print("Lane Emden was called with densities {0} via RK4".format(densities))
                                    else:
                                        print("Please end the command with the numerical method of your choice.")
                                        sys.exit(ec.INVALID_ARGUMENT)
                                elif argv[index] == P.EULER[0] or argv[index] == P.EULER[1]:
                                    if method is not None:
                                        print("Please provide the numerical method once.")
                                        sys.exit(ec.INVALID_COMMAND_ARG)
                                    if index == argc - 1:
                                        method = "euler"
                                        print("Lane Emden was called with densities {0} via Euler".format(densities))
                                    else:
                                        print("Please end the command with the numerical method of your choice.")
                                        sys.exit(ec.INVALID_ARGUMENT)
                                else:
                                    print("Density '{0}' must be a number.".format(argv[index]))
                                    sys.exit(ec.INVALID_ARGUMENT)
                        if method is None:
                            print("Please provide a numerical method.")
                            sys.exit(ec.MISSING_ARGUMENT)
                    else:
                        print("Please enter at least one valid real positive density.")
                        sys.exit(ec.MISSING_ARGUMENT)
                else:
                    for index in range(2, argc):
                        if checkDensityMinimum(argv[index]) is False and \
                                checkDensityMaximum(argv[index]) is False and \
                                checkDensityStep(argv[index]) and False:
                            if argv[index] == P.RK4[0] or argv[index] == P.RK4[1]:
                                if index == argc - 1:
                                    if method is None:
                                        method = "rk4"
                                    else:
                                        print("Please provide the numerical method once.")
                                        sys.exit(ec.INVALID_COMMAND_ARG)
                                else:
                                    print("Please end the command with the numerical method of your choice.")
                                    sys.exit(ec.INVALID_ARGUMENT)
                            elif argv[index] == P.EULER[0] or argv[index] == P.EULER[1]:
                                if index == argc - 1:
                                    if method is None:
                                        method = "euler"
                                    else:
                                        print("Please provide the numerical method once.")
                                        sys.exit(ec.INVALID_COMMAND_ARG)
                                else:
                                    print("Please end the command with the numerical method of your choice.")
                                    sys.exit(ec.INVALID_ARGUMENT)
                            else:
                                print("{0} is not a valid argument.\n"
                                      "Use '{1} -h' or '{1} --help' for more information."
                                      .format(argv[index], argv[0]))
                                sys.exit(ec.INVALID_COMMAND_ARG)

                    densitySteps = np.asarray(densitySteps)
                    for index, val in enumerate(densitySteps):
                        if val == -1:
                            matchMissingDensityInfo(index)
                    if method is None:
                        print("Please provide a numerical method.")
                        sys.exit(ec.MISSING_ARGUMENT)

            else:
                print("Please enter a parameter for density.")
                sys.exit(ec.MISSING_ARGUMENT)
        elif argv[1] == P.POLYTROPIC[0] or \
                argv[1] == P.POLYTROPIC[1]:

        else:
            print("{0} is not a valid argument.\n"
                  "Use '{1} -h' or '{1} --help' for more information."
                  .format(argv[1], argv[0]))
            sys.exit(ec.INVALID_COMMAND_ARG)