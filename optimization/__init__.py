"""This module contains the optimization algorithm and the settings 
    """

import math


def run(F_desired,P_max, settings):
    """ Algorithm to determine the thruster with the best specific impulse for the given thrust and power requirements

    Args:
        F_desired (N): Required thrust
        P_max (W): Maximum power consumption
        In settings {dict}:
            FDP (dict): Fixed design parameters
            bounds (dict of 2-tuples): The constraints on the design parameters to be optimized
            NR (dict of functions): Nusselt relations that are used in various phases in the channel (liquid, two-phase, dry-out,gas)
            PDR (dict of functions): Pressure drop relations that are used in various phases in the channel (friction: liquid, two-phase, gas; contraction)
            steps (dict of ints): Fidelity of the temperature-steps in the various sections of the channel
    """

    ## First determine ideal engine performance

    print(settings)

