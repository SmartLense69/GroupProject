class Param:
    LANE_EMDEN:     list[str] = ["-le", "--lane-emden"]
    POLYTROPIC:     list[str] = ["-p", "--polytropic"]
    WHITE_DWARF:    list[str] = ["-w", "--whiteDwarf"]
    NEUTRON_STAR:   list[str] = ["-ns", "--neutronStar"]
    RK2:            list[str] = ["-rk2", "--runge-kutta-2"]
    RK4:            list[str] = ["-rk4", "--runge-kutta-4"]
    RKF:            list[str] = ["-rkf", "--runge-kutta-fehlberg"]
    EULER:          list[str] = ["-e", "--euler"]
    DENSITY:        list[str] = ["-rho=", "--density"]
    DENSITY_MIN:    list[str] = ["-rm=", "--minDensity="]
    DENSITY_MAX:    list[str] = ["-rM=", "--maxDensity="]
    DENSITY_STEP:   list[str] = ["-rS=", "--stepDensity="]
    DENSITY_NUM:    list[str] = ["-rN=", "--numDensity="]
    HELP:           list[str] = ["-h", "--help"]
    ANALYTICAL:     list[str] = ["-a", "--analytical"]
    SHOW_EOS:       list[str] = ["-es", "--showEOS"]
    N_SINGLE: str = "-n="
    N_MIN: str = "-nMin="
    N_MAX: str = "-nMax="
    N_STEP:   str = "-nH="
