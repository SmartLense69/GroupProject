import sys
import CR_exitcodes as ec
from CR_param import Param as P

if __name__ == '__main__':
    argv: list[str] = sys.argv
    argc: int = len(argv)
    if argc - 1 == 0:

        print("No command line arguments have been given.\n" \
              "Use '{0} -h' or '{0} --help' for more information." \
              .format(sys.argv[0]))
        sys.exit(ec.NO_ARGUMENTS)

    else:

        if argv[1] is P.LANE_EMDEN[0] or \
                argv[1] is P.LANE_EMDEN[1]:

            if argc - 2 >= 1:

                if argv[2] is P.DENSITY[0] or \
                        argv[2] is P.DENSITY[1]:

                    if argc - 3 >= 1:
                        for index in range(3, argc):
                            try:
                                density = float(argv[index])
                                if density > 0:
                                    print("Lane Emden was called with {0}".format(density))
                                else:
                                    print("{0} is not a valid density. Densities must be real positive."
                                          .format(density))
                                    sys.exit(ec.INVALID_ARGUMENT)
                            except ValueError:
                                print("Argument '{0}' must be a number.".format(argv[index]))
                                sys.exit(ec.INVALID_ARGUMENT)
                    else:
                        print("Please enter at least one valid real positive density.")
                        sys.exit(ec.MISSING_ARGUMENT)

            else:
                print("Please enter a parameter for density.")
                sys.exit(ec.MISSING_ARGUMENT)
