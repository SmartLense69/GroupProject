"""
    In this module, all exit codes are defined.
"""

SUCCESS:            int = 0
"""
    If the application had no issues running.
"""

NO_ARGUMENTS:       int = 1
"""
    If the user provided no arguments.
"""

INVALID_ARGUMENT:   int = 2
"""
    If the user provided an invalid argument.
"""

MISSING_ARGUMENT:   int = 3
"""
    If the user forgot to provide an argument.
"""

TOO_MANY_ARGUMENTS: int = 4
"""
    If the user specified too many arguments.
"""

INVALID_PARAMETER:  int = 5
"""
    If the user put an invalid parameter, like a negative density.
"""
