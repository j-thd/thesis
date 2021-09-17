"""This module contains the optimization algorithm and the settings 
    """

import math
from basic import IRT
from basic import chamber
from models.optimization_1D import optim_P_total
from models import one_D as oneD
from scipy.optimize import minimize, Bounds

def run(F_desired,T_chamber, channel_amount, settings, x_guess):
    """ Algorithm to determine the thruster with the best specific impulse for the given thrust and power requirements

    Args:
        F_desired (N): Required thrust
        T_chamber (K): Chamber temperature determines mass flow and Isp and ideal power consumption
        channel_amount (-): Amount of channels. Technically, this must be optmized, but using built-in optimization really doesn't work well for integers
        In settings {dict}:
            FDP (dict): Fixed design parameters
            bounds (dict of 2-tuples): The constraints on the design parameters to be optimized
            NR (dict of functions): Nusselt relations that are used in various phases in the channel (liquid, two-phase, dry-out,gas)
            PDR (dict of functions): Pressure drop relations that are used in various phases in the channel (friction: liquid, two-phase, gas; contraction)
            steps (dict of ints): Fidelity of the temperature-steps in the various sections of the channel
    """
    print("\n  ----- OPTIMIZING THRUSTER FOR POWER CONSUMPTION -----")
    print("\n    HIGH-LEVEL INPUTS")
    print("--- Desired thrust:        {:3.1f} mN".format(F_desired*1e3))
    print("--- Chamber temperature:   {:4.1f} K".format(T_chamber))
    
    ## First determine ideal engine performance, and load settings for that
    T_inlet = settings['FDP']['T_inlet'] # [K] Inlet temperature (at heating channel inlet manifold)
    p_inlet = settings['FDP']['p_inlet'] # [Pa] Chamber pressure assumed equal to inlet pressure
    AR_exit = settings['FDP']['AR_exit'] # [-] Exit area ratio of nozzle
    p_back = settings['FDP']['p_back'] # [Pa]
    assert p_back == 0 # Nozzle correction does not work for back pressure
    fp = settings['FDP']['fp'] # [object] FluidProperties object containing thermodynamic parameters for propellant

    print("--- Inlet pressure:          {:2.1f} bar".format(p_inlet*1e-5))
    print("--- Exit area ratio          {:2.1f}".format(AR_exit))

    


    # Check if the inputs are correct
    assert(T_chamber > T_inlet)
    ep_ideal = IRT.engine_performance_from_F_and_T(
        F_desired = F_desired,
        p_chamber=p_inlet,
        T_chamber=T_chamber,
        AR_exit=AR_exit,
        p_back=p_back,
        fp=fp
    )

    # Calculate ideal power consumption
    m_dot = ep_ideal['m_dot'] # [kg/s] Mass flow through chamber
    Isp = ep_ideal['Isp'] # [s] Specific impulse
    P_ideal = chamber.ideal_power_consumption( # [W] Ideal power consumption
        mass_flow=m_dot, T_inlet=T_inlet,p_inlet=p_inlet, T_outlet=T_chamber,p_outlet=p_inlet, fp=fp)

    print("\n    IDEAL PERFORMANCE")
    print("--- Nozzle status:           {}".format(ep_ideal['nozzle_status']))
    print("--- Mass flow:               {:3.2f} mg/s".format(m_dot*1e6))
    print("--- Isp:                     {:3.1f} s".format(Isp))
    print("--- Power consumption:       {:2.1f} W".format(P_ideal))



    # Calculate pre-calculated values of heating channel that are constant through entire optimization
    prepared_values = oneD.full_homogenous_preparation(
        T_inlet=T_inlet,                            # [K] Inlet temperature of the channel
        T_outlet=T_chamber,                         # [K] Outlet of the channel is equal to chamber temperature in IRT
        m_dot=m_dot,                                # [kg/s] Mass flow
        p_ref=p_inlet,                              # [Pa] The pressure in the channel is assumed to be constant for heat flow calculation purposes
        steps_l=settings['steps']['steps_l'],       # [-] Fidelity of simulation in liquid section of channel
        steps_tp=settings['steps']['steps_tp'],     # [-] Fidelity of simulation in two-phase section of channel
        steps_g=settings['steps']['steps_g'],       # [-] Fidelity of simulation in gas section of channel
        fp=fp,                                      # [object] Object to access thermodynamic parameters of propellant with
    )

    ### Define function to optimize here, and load settings into fixed design parameters
    ### IMPORTANT NOTE: if the order of arguments of f() is changed, the order of bounds must be changed too. 
    # The bounds are defined a bit early here, so one does not forgot this after changing either.
    b = Bounds(
        [ # Lower bounds
            settings['bounds']['w_channel'][0],
            settings['bounds']['w_channel_spacing'][0],
            settings['bounds']['T_wall_superheat'][0] + T_chamber,
        ],
        [ # Higher bounds
            settings['bounds']['w_channel'][1],
            settings['bounds']['w_channel_spacing'][1],
            settings['bounds']['T_wall_superheat'][1] + T_chamber,
        ])
    # Initial guess (just make it the higher bound (arbitary)) NOTE: SAME ORDER APPLIES ABOUNDS AND AS f() below
    x0 = None # Set the first guess based on a guess being provided or not
    if (x_guess is None):
        x0 = [ 
            settings['bounds']['w_channel'][1],
            settings['bounds']['w_channel_spacing'][1],
            settings['bounds']['T_wall_superheat'][1] + T_chamber,
            ]
    else:
        x0 = x_guess
    # The function to optimize
    def f(x): #  NOTE: change bounds above if argument change
        w_channel = x[0]
        w_channel_spacing = x[1]
        T_wall = x[2]
        # Wall arguments need to be encapsulated in a dictionary so they can be conveniently passed to the section of code calcuating wall effects
        # and the actual bottom wall temperature
        wall_args = {
            'kappa_wall': settings['FDP']['kappa_wall'],                            # [W/(m*K)] Thermal conductivity of chip wall
            'h_channel': settings['FDP']['h_channel'],                              # [m] Channel depth, channel height
            'w_channel': w_channel,                                                 # [m] Channel width
            'w_channel_spacing': w_channel_spacing,                                 # [m] Spacing between channels/wall thickness
            'emissivity_chip_bottom': settings['FDP']['emissivity_chip_bottom'],    # [m] Emissivity of the bottom of the chip
        }
        # Exit manifold length is zero based on existing VLM design, so should be zero in settings
        assert settings['FDP']['l_exit_manifold'] == 0 # It is merely a remnant from previous code choices
        # The function that returns the eventual power output
        res_P = optim_P_total(                                                                      # CAPITALIZED COMMENTS ARE FOR OPTIMIZED VARIABLES
            channel_amount=channel_amount,                                                  # [-] AMOUNT OF CHANNELS
            w_channel=w_channel,                                                            # [m] CHANNEL WIDTH
            h_channel=settings['FDP']['h_channel'],                                         # [m] Channel depth, channel height
            inlet_manifold_length_factor=settings['FDP']['inlet_manifold_length_factor'],   # [-] Linear factor to determine inlet manifold length
            inlet_manifold_width_factor=settings['FDP']['inlet_manifold_width_factor'],     # [-] Linear factor to determine inlet manifold width
            l_exit_manifold=settings['FDP']['l_exit_manifold'],                             # [-] NOTE: old part of design. Is set to zero.
            w_channel_spacing=w_channel_spacing,                                            # [m] SPACING BETWEEN CHANNELS/WALL THICKNESS
            w_outer_margin=settings['FDP']['w_outer_margin'],                               # [m] Margin around edge of heating channels, or nozzle width
            T_wall=T_wall,                                                                  # [K] WALL TEMPERATURE
            p_ref=p_inlet,                                                                  # [Pa] Pressure through channel
            m_dot=m_dot,                                                                    # [kg/s] Heat flow through channel
            prepared_values=prepared_values,                                                # {dict} Dictionary with many parameters that are costly to calcuate but are cosntant throughout optimization
            Nusselt_relations=settings['Nusselt_relations'],                                # {dict} Heat transfer relations used for optimization
            pressure_drop_relations=settings['pressure_drop_relations'],                    # {dict} Pressure drops used for optimization
            convergent_half_angle=settings['FDP']['convergent_half_angle'],                 # [rad] Half-angle of convergent part of 2D-conical nozzle
            divergent_half_angle=settings['FDP']['divergent_half_angle'],                   # [rad] Half-angle of divergent part of 2D-conical nozzle
            F_desired=F_desired,                                                            # [N] Desired thrust. Still required to correct for pressure drop
            p_back=p_back,                                                                  # [Pa] Back pressure. NOTE: Nozzle correction depends on this being 0.
            AR_exit=AR_exit,                                                                # [-] Exit area ratio to correct nozzle after pressure drop
            emissivity_top=settings['FDP']['emissivity_chip_top'],                               # [-] Emissivity of top chip (black-body-radiation)
            fp=fp,                                                                          # [object] Object to access thermodynamic parameters of propellant with
            wall_args=wall_args                                                             # {dict} Arguments about wall design. See above
        )
        # Return the total power consumption. The variable of interest
        P_total = P_ideal + res_P['P_loss']
        #print(x)
        
        # Rescaling
        factor = 1e6

        if math.isnan(P_total):
            #print("A constraint violated.")
            # Punish the objective function by outputting high power consumptions for invalid solutions
            # The punishment must however not be constant, or it will converge at the boundary where it is invalid.
            # Since the punishment stems from the pressure drop being higher than the initial p_chamber value, the lower the negative final pressure is the higher the punishment
            return (P_ideal*(2+1*res_P['pressure_drop_punishment']))*factor
        else:
            #print("P_total: {:2.5f} W".format(P_total))
            return P_total*factor

    ## OPTIMIZATION ALGORITHM
    # First bounds must be set in the same order that arguments are given in f()
    # The order is: channel_amount, w_channel,w_channel_spacing,T_wall
    minimize_results = minimize(
        fun=f,
        x0=x0,
        bounds=b
        )

    optim_results = {
        'minimize_results': minimize_results,
        'w_channel': minimize_results.x[0],
        'w_channel_spacing': minimize_results.x[1],
        'T_wall_top': minimize_results.x[2],
        'P_total': minimize_results.fun/1e6,
    }
    print(optim_results)
    return optim_results

        


