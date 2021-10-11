# File to plot optimization results from saved files

import math
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt

import optimization
import optimization.settings
from basic import IRT
from basic import chamber
from models.optimization_1D import optim_P_total
from models import one_D as oneD
from thermo.prop import FluidProperties

def run():
    F_desired = 4e-3 # [mN] Thrust
    T_chamber = 1000 # [K] Chamber temperature
    file_name = "optimization_results/optimization_results-F{:1.0f}mN-{}K.npz".format(F_desired*1e3,T_chamber)
    optimization_settings = optimization.settings.settings_1D_rectangular_multichannel # For the analysis of the answer

    npzfile = open(file_name, "rb")
    data = np.load(npzfile)

    # The Isp and m_dot is the same for all cases for a given thrust F and temperature T. It was only saved as an array for convenience
    F_desired = data['F_desired'] # Actually read the fle just to check
    T_chamber = data['T_chamber']
    m_dot = data['m_dot'][0]
    Isp = data['Isp'][0]
    p_inlet = data['p_inlet'][0]

    print("\n------ PLOTTING RESULTS FOR:     ")
    print(" - F:                {:3.1f} mN".format(F_desired*1e3))
    print(" - m_dot:            {:3.3f} mg/s".format(m_dot*1e6))
    print(" - Isp:              {:3.1f} s".format(Isp))
    print(" - T_chamber:        {:3.1f} K".format(T_chamber))

    # Say where minimal power consumption was found.
    id = np.argmin(data['P_loss']) # Index of result with lowest power consumption
    P_min = data['P_loss'][id] # [W] Minimal power consumption
    channel_min = data['channel_amount_range'][id] # [-] Amount of channels
    w_channel_min = data['w_channel'][id]
    w_channel_spacing_min = data['w_channel_spacing'][id]
    T_wall_top_min = data['T_wall'][id]
    print('\n----- BEST RESULT:')
    print(" Power consumption:       {:2.3f} W".format(P_min))
    print(" Amount of channels:      {:2.0f}".format(channel_min))
    print(" Channel width:           {:3.4f} micron".format(w_channel_min*1e6) )
    print(" Channel spacing:         {:3.4f} micron".format(w_channel_spacing_min*1e6) )
    print(" Top wall temperature:    {:3.4f} K".format(T_wall_top_min))

    # Determine parameters at best result
    Re_l_min = data['Re_channel_l'][id]
    Re_tp_min = data['Re_channel_tp'][id]
    Re_g_min = data['Re_channel_g'][id]
    Re_throat_min = data['Re_throat_new'][id]

    l_total_min = data['l_total'][id]
    l_channel_l_min = data['l_channel_l'][id]
    l_channel_tp_min = data['l_channel_tp'][id]
    l_channel_g_min = data['l_channel_g'][id]

    print("\n---- PARAMETERS AT BEST RESULT:")
    print(" - Reynolds:")
    print("  * Liquid:    {:3.1f}".format(Re_l_min))
    print("  * Two-phase: {:3.1f}".format(Re_tp_min))
    print("  * Gas:       {:3.1f}".format(Re_g_min))
    print("  * Throat:    {:3.1f}".format(Re_throat_min))
    print(" - Channel section length:")
    print("  * Liquid:    {:3.1f} micron".format(l_channel_l_min*1e6))
    print("  * Two-phase: {:3.1f} micron".format(l_channel_tp_min*1e6))
    print("  * Gas:       {:3.1f} micron".format(l_channel_g_min*1e6))
    


    # Analyze best result for sensitivity.

    # Determine left/rightbounds of optimum
    sensitivity = 0.05
    lb_id = id
    # Find the first value that is higher and then stop
    while data['P_loss'][lb_id] < P_min*(1+sensitivity):
        lb_id -= 1
    
    lb = data['channel_amount_range'][lb_id]
    print("Lower bound optimum: {} channels".format(lb))

    ub_id = id
    # Find the first value that is higher and then stop
    while data['P_loss'][ub_id] < P_min*(1+sensitivity):
        ub_id += 1
    
    ub = data['channel_amount_range'][ub_id]
    print("Higher bound optimum: {} channels".format(ub))

    str_title = "({:1.1f} mN, $p_{{in}}=${:1.0f} bar, $T_c=${:3.0f} K)".format(F_desired*1e3, p_inlet*1e-5, T_chamber)
    optimum_bounds = {
        'min': channel_min,
        'lb': lb,
        'ub': ub,
    }
    # plotReynolds(
    #     channel_amount=data['channel_amount_range'],
    #     Re_l=data['Re_channel_l'],
    #     Re_tp=data['Re_channel_tp'],
    #     Re_g=data['Re_channel_g'],
    #     Re_throat=data['Re_throat_new'],
    #     str_title=str_title,
    #     optimum_bounds=optimum_bounds,
    #     sens=sensitivity
    # )

    # plotOptimumOutcome(
    #     channel_amount=data['channel_amount_range'],
    #     P_loss=data['P_loss'],
    #     w_channel=data['w_channel'],
    #     w_channel_spacing=data['w_channel_spacing'],
    #     T_wall=data['T_wall'],
    #     str_title=str_title,
    #     optimum_bounds=optimum_bounds,
    #     sens=sensitivity)

    # plotChannelCharacteristics(
    #     channel_amount=data['channel_amount_range'],
    #     pressure_drop=data['pressure_drop'],
    #     l_channel=data['l_total'],
    #     w_throat_new=data['w_throat_new'],
    #     w_total=data['w_total'],
    #     str_title=str_title,
    #     optimum_bounds=optimum_bounds,
    #     sens=sensitivity
    # )

    #plotPressureDropVsThroatReynolds(pressure_drop=data['pressure_drop'], Re_throat=data['Re_throat_new'], str_title=str_title, optimum_bounds=optimum_bounds, sens=sensitivity)

    #analyze_best_result(id=id,data=data,settings=optimization_settings, str_title=str_title)

    monte_carlo_guess(id=id, data=data,settings=optimization_settings, str_title=str_title)

    plt.show()

def plotChannelCharacteristics(channel_amount, pressure_drop, l_channel, w_throat_new, w_total, str_title, optimum_bounds, sens):
    fig, axs = plt.subplots(2,2)
    # Total power consumption
    axs[0][0].plot(channel_amount, pressure_drop*1e-5)
    axs[0][0].set_ylabel("$\Delta P$ [bar]")
    axs[0][0].grid()
    plotOptimumBounds(optimum_bounds,sens,axs=axs[0][0],labels=True)
    axs[0][1].plot(channel_amount, l_channel*1e3)
    axs[0][1].set_ylabel("$l_c$ [mm]")
    axs[0][1].grid()
    plotOptimumBounds(optimum_bounds,sens,axs=axs[0][1])
    axs[1][0].plot(channel_amount, w_throat_new*1e6)
    axs[1][0].set_ylabel("$w_t$ [$\\mu$m]")
    axs[1][0].grid()
    plotOptimumBounds(optimum_bounds,sens,axs=axs[1][0])
    axs[1][1].plot(channel_amount, w_total*1e3)
    axs[1][1].set_ylabel("$w_{{chip}}}$ [mm]")
    axs[1][1].grid()
    plotOptimumBounds(optimum_bounds,sens, axs=axs[1][1])
    fig.suptitle("Geometry "+str_title)
    axs[0][0].legend()
    plt.tight_layout(pad=0.5)

def plotOptimumOutcome(channel_amount, P_loss, w_channel, w_channel_spacing, T_wall, str_title, optimum_bounds, sens):
    fig, axs = plt.subplots(2,2)
    # Total power consumption
    axs[0][0].plot(channel_amount, P_loss)
    axs[0][0].set_ylabel("$P_{{loss}}$ [W]")
    axs[0][0].grid()
    plotOptimumBounds(optimum_bounds,sens,axs=axs[0][0],labels=True)
    # Channel width
    axs[0][1].plot(channel_amount, w_channel*1e6)
    axs[0][1].set_ylabel("$w_c$ [$\\mu$m]")
    axs[0][1].grid()
    plotOptimumBounds(optimum_bounds,sens,axs=axs[0][1])
    # Channel width spacing
    axs[1][0].plot(channel_amount, w_channel_spacing*1e6)
    axs[1][0].set_ylabel("$s_c$ [$\\mu$m]")
    axs[1][0].grid()
    plotOptimumBounds(optimum_bounds,sens,axs=axs[1][0])
    # Top wall temperature
    axs[1][1].plot(channel_amount, T_wall)
    axs[1][1].set_ylabel("$T_{{wall}}$ [K]")
    axs[1][1].grid()
    plotOptimumBounds(optimum_bounds,sens,axs=axs[1][1])
    fig.suptitle("Optimum outcome and parameters "+str_title)
    axs[0][0].legend()
    plt.tight_layout(pad=0.5)

def plotPressureDropVsThroatReynolds(pressure_drop, Re_throat, str_title, optimum_bounds, sens):
    sens_label = "+{:1.0f}% loss".format(sens*100)
    plt.figure()

    plt.plot(pressure_drop*1e-5,Re_throat, label="Throat")
    plt.axvline(1e-5*pressure_drop[optimum_bounds['min']], label="Optimum", c="red", linestyle="dashed")
    plt.axvline(1e-5*pressure_drop[optimum_bounds['lb']], label=sens_label, c="red", linestyle="dotted")
    plt.axvline(1e-5*pressure_drop[optimum_bounds['ub']],  c="red", linestyle="dotted")
    plt.xlabel("Pressure drop $\Delta p$ [bar]")
    plt.ylabel("Throat Reynolds $Re_{{throat,dh}}$ [-]")
    plt.title("Throat Reynolds vs. pressure drop " + str_title)
    plt.grid()
    plt.legend()
    plt.tight_layout()

    print("\n Throat Reynold's vs. pressure drop info:")
    Re_lb = Re_throat[optimum_bounds['lb']]
    Re_min = Re_throat[optimum_bounds['min']]
    Re_ub = Re_throat[optimum_bounds['ub']]
    print("Re_lb {:4.0f} ({:2.1f})%".format(Re_lb, (Re_lb-Re_min)/Re_min*100))
    print("Re_min {:4.0f}".format(Re_min))
    print("Re_ub {:4.0f} ({:2.1f})%".format(Re_ub, (Re_ub-Re_min)/Re_min*100))

def plotReynolds(channel_amount, Re_l, Re_tp, Re_g, Re_throat, str_title, optimum_bounds, sens):
    plt.figure()
    
    plt.plot(channel_amount, Re_l, label="Liquid")
    plt.plot(channel_amount, Re_tp, label="Two-Phase")
    plt.plot(channel_amount, Re_g, label="Gas")
    plt.plot(channel_amount, Re_throat, label="Throat",c='black')
    plt.xlabel("Number of channels $N_c$[-]")
    plt.ylabel("$Re_{{D,hydraulic}}$ [-]")
    plotOptimumBounds(optimum_bounds, sens)
    plt.yscale('log')
    plt.grid()
    plt.title("Reynolds at start of section/throat "+str_title)
    plt.legend()
    plt.tight_layout()

def plotOptimumBounds(optimum_bounds, sens, axs = None, labels=False):
    sens_label = "+{:1.0f}% loss".format(sens*100)
    if axs is None:
        plt.axvline(optimum_bounds['min'], label='Optimum', color='red', linestyle='dashed')
        plt.axvline(optimum_bounds['lb'], label=sens_label, color='red', linestyle='dotted')
        plt.axvline(optimum_bounds['ub'], color='red', linestyle='dotted')
    elif labels == True:
        axs.axvline(optimum_bounds['min'], label='Optimum', color='red', linestyle='dashed')
        axs.axvline(optimum_bounds['lb'], label=sens_label, color='red', linestyle='dotted')
        axs.axvline(optimum_bounds['ub'], color='red', linestyle='dotted')
    else:
        axs.axvline(optimum_bounds['min'], color='red', linestyle='dashed')
        axs.axvline(optimum_bounds['lb'], color='red', linestyle='dotted')
        axs.axvline(optimum_bounds['ub'],  color='red', linestyle='dotted')

def monte_carlo_guess(id, data, settings, str_title, attempts=250):
    """Check iof the outcome of the optimization is dependent on the guess

    Args:
        id (): [description]
        data ([type]): [description]
        settings ([type]): [description]
        attempts (int, optional): [description]. Defaults to 100.
    """
        # Vary the parameters other than channel_amount to see if it is locally an optimum
    channel_amount = data['channel_amount_range'][id]
    P_loss_min = data['P_loss'][id]
    w_channel_min = data['w_channel'][id]
    w_channel_spacing_min = data['w_channel_spacing'][id]
    T_wall_min = data['T_wall'][id]

    F_desired = data['F_desired'] # Actually read the fle just to check
    T_chamber = float(data['T_chamber'])

    # Generate random guesses as x_guess input
    np.random.seed(20211006) # Current date set as random seed for reproducibility
    # For w_channel
    w_channel_lb = settings['bounds']['w_channel'][0]
    w_channel_ub = settings['bounds']['w_channel'][1]
    w_channel_guesses = w_channel_lb + np.random.rand(attempts)*(w_channel_ub-w_channel_lb)
    # For w_channel_spacing
    w_channel_spacing_lb = settings['bounds']['w_channel_spacing'][0]
    w_channel_spacing_ub = settings['bounds']['w_channel_spacing'][1]
    w_channel_spacing_guesses = w_channel_spacing_lb + np.random.rand(attempts)*(w_channel_spacing_ub-w_channel_spacing_lb)
    # For T_wall
    T_wall_lb = settings['bounds']['T_wall_superheat'][0] + T_chamber
    T_wall_ub = settings['bounds']['T_wall_superheat'][1] + T_chamber
    T_wall_guesses = T_wall_lb + np.random.rand(attempts)*(T_wall_ub-T_wall_lb)

    x_guess_array = np.vstack((w_channel_guesses,w_channel_spacing_guesses,T_wall_guesses)).transpose()
    print(x_guess_array)
    print(x_guess_array[0])

    P_loss = np.zeros(attempts)
    w_channel = np.zeros_like(P_loss)
    w_channel_spacing = np.zeros_like(P_loss)
    T_wall = np.zeros_like(P_loss)

    iter_P = np.nditer(P_loss, flags=['c_index'])

    for iP in iter_P:
        print("\nRandomized guess: {}".format(iter_P.index))
        print("x_guess: {}".format(x_guess_array[iter_P.index]))
        res = optimization.run(
        F_desired=F_desired,
        T_chamber=T_chamber,
        channel_amount=channel_amount,
        settings=settings,
        x_guess=x_guess_array[iter_P.index]
        )

        P_loss[iter_P.index] = res['P_loss']
        w_channel[iter_P.index] = res['w_channel']
        w_channel_spacing[iter_P.index] = res['w_channel_spacing']
        T_wall[iter_P.index] = res['T_wall_top']


    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    sensitivity = 1e-6 # [-] Whether found solution with different guess is within this bounds
    threshold = P_loss_min*(1+sensitivity)
    P_loss_mask_correct = P_loss >= threshold
    print(P_loss_mask_correct)
    P_loss_mask_wrong = P_loss < threshold
    P_loss_correct = ma.array(P_loss, mask=P_loss_mask_correct)
    print(P_loss_correct)
    print("Mean + STD correct: {} ({})".format(np.mean(P_loss_correct),np.std(P_loss_correct)))
    P_loss_wrong = ma.array(P_loss, mask=P_loss_mask_wrong)
    print(P_loss_wrong)
    print("Mean + STD wrong: {} ({})".format(np.mean(P_loss_wrong),np.std(P_loss_wrong)))

    # Masked array with the correct guesses
    w_channel_correct = ma.array(w_channel_guesses, mask = P_loss_mask_correct)
    print(w_channel_correct)
    w_channel_spacing_correct = ma.array(w_channel_spacing_guesses, mask = P_loss_mask_correct)
    print(w_channel_spacing_correct)
    T_wall_correct = ma.array(T_wall_guesses, mask = P_loss_mask_correct)

    # Masked array with the correct results
    w_channel_correct_results = ma.array(w_channel, mask = P_loss_mask_correct)
    print(w_channel_correct_results)
    w_channel_spacing_correct_results = ma.array(w_channel_spacing, mask = P_loss_mask_correct)
    print(w_channel_spacing_correct_results)
    T_wall_correct_results = ma.array(T_wall, mask = P_loss_mask_correct)
    print(T_wall_correct_results)
    print("\n --- CORRECT RESULTS  ({})---".format(ma.count_masked(P_loss_wrong)))
    print(" P_loss:       MIN: {:1.8F} W       MAX: {:1.8F} W     ".format(np.min(P_loss_correct), np.max(P_loss_correct)))
    print(" w_channel:    MIN: {:3.6f} micron  MAX: {:3.6f} micron".format(np.min(w_channel_correct_results*1e6),np.max(w_channel_correct_results*1e6)))
    print(" spacing:      MIN: {:3.6f} micron  MAX: {:3.6f} micron".format(np.min(w_channel_spacing_correct_results*1e6),np.max(w_channel_spacing_correct_results*1e6)))
    print(" T_wall:       MIN: {:4.5f} K       MAX: {:3.6f} K     ".format(np.min(T_wall_correct_results),np.max(T_wall_correct_results)))

    # Masked array with the wrong guesses
    w_channel_wrong = ma.array(w_channel_guesses, mask = P_loss_mask_wrong)
    w_channel_spacing_wrong = ma.array(w_channel_spacing_guesses, mask = P_loss_mask_wrong)
    T_wall_wrong = ma.array(T_wall_guesses, mask = P_loss_mask_wrong)

    # Masked array with the wrong guesses
    w_channel_wrong_results = ma.array(w_channel, mask = P_loss_mask_wrong)
    w_channel_spacing_wrong_results = ma.array(w_channel_spacing, mask = P_loss_mask_wrong)
    T_wall_wrong_results = ma.array(T_wall, mask = P_loss_mask_wrong)

    print("\n --- WRONG RESULTS ({})---".format(ma.count_masked(P_loss_correct)))
    print(" P_loss:       MIN: {:1.8F} W       MAX: {:1.8F} W     ".format(np.min(P_loss_wrong), np.max(P_loss_wrong)))
    print(" w_channel:    MIN: {:3.6f} micron  MAX: {:3.6f} micron".format(np.min(w_channel_wrong_results*1e6),np.max(w_channel_wrong_results*1e6)))
    print(" spacing:      MIN: {:3.6f} micron  MAX: {:3.6f} micron".format(np.min(w_channel_spacing_wrong_results*1e6),np.max(w_channel_spacing_wrong_results*1e6)))
    print(" T_wall:       MIN: {:4.5f} K       MAX: {:3.6f} K     ".format(np.min(T_wall_wrong_results),np.max(T_wall_wrong_results)))


    
    ax.scatter(w_channel_correct*1e6, w_channel_spacing_correct*1e6, T_wall_correct, marker="o", color="blue", label="$P_{{loss}}$ within {:1.4f}%".format(sensitivity*100))
    ax.scatter(w_channel_wrong*1e6, w_channel_spacing_wrong*1e6, T_wall_wrong, marker="X", color="red", label="$P_{{loss}}$ outside {:1.4f}%".format(sensitivity*100))
    ax.set_xlabel("Channel width - $w_c$ [$\\mu$m]")
    ax.set_ylabel("Channel spacing - $s_c$ [$\\mu$m]")
    ax.set_zlabel("Top wall temperature - $T_{{wall,top}}$ [K]")
    ax.legend()
    fig.suptitle("Monte-Carlo Sensitivity Analysis on Starting Guesses\n"+str_title)
    

def analyze_best_result(id, data, settings, str_title, resolution=100 ):
    """
        Id is the number of the optimum
    """   
        # Vary the parameters other than channel_amount to see if it is locally an optimum
    channel_amount = data['channel_amount_range'][id]
    P_loss_min = data['P_loss'][id]
    w_channel_min = data['w_channel'][id]
    w_channel_spacing_min = data['w_channel_spacing'][id]
    T_wall_min = data['T_wall'][id]

    F_desired = data['F_desired'] # Actually read the fle just to check
    T_chamber = float(data['T_chamber'])

    w_channel_range = np.linspace(start=settings['bounds']['w_channel'][0],stop=settings['bounds']['w_channel'][1], num=resolution)
    w_channel_spacing_range = np.linspace(start=settings['bounds']['w_channel_spacing'][0],stop=settings['bounds']['w_channel_spacing'][1], num=resolution)
    T_wall_range = np.linspace(start=settings['bounds']['T_wall_superheat'][0]+T_chamber,stop=settings['bounds']['T_wall_superheat'][1]+T_chamber, num=resolution)

    P_loss_channel = np.zeros_like(w_channel_range)
    P_loss_channel_spacing = np.zeros_like(w_channel_range)
    P_loss_wall = np.zeros_like(w_channel_range)

    # Calculate the alternative outcomes
    iter_P = np.nditer(P_loss_channel, flags=['c_index'])
    for iP in iter_P:
        # Calculate results for alternative w_channel
        res_channel = alternative_power_loss(\
            F_desired=F_desired,
            T_chamber=T_chamber,
            channel_amount=channel_amount,
            settings=settings,
            w_channel=w_channel_range[iter_P.index],
            w_channel_spacing=w_channel_spacing_min,
            T_wall=T_wall_min)
        # Calculate result for alternative channel spacing
        res_channel_spacing = alternative_power_loss(\
            F_desired=F_desired,
            T_chamber=T_chamber,
            channel_amount=channel_amount,
            settings=settings,
            w_channel=w_channel_min,
            w_channel_spacing=w_channel_spacing_range[iter_P.index],
            T_wall=T_wall_min)
        # Calculate result for alternative T_wall
        res_wall = alternative_power_loss(\
            F_desired=F_desired,
            T_chamber=T_chamber,
            channel_amount=channel_amount,
            settings=settings,
            w_channel=w_channel_min,
            w_channel_spacing=w_channel_spacing_min,
            T_wall=T_wall_range[iter_P.index])

        # Store alternative power outcomes
        P_loss_channel[iter_P.index]=res_channel['P_loss']
        P_loss_channel_spacing[iter_P.index]=res_channel_spacing['P_loss']
        P_loss_wall[iter_P.index]=res_wall['P_loss']

    fig, axs = plt.subplots(2,2)
    # P_loss vs channel amount
    axs[0][0].plot(data['channel_amount_range'], data['P_loss'])
    axs[0][0].set_xlabel("Number of channels $N_c$ [-]")
    axs[0][0].set_ylabel("$P_{{loss}}$ [W]")
    axs[0][0].axvline(data['channel_amount_range'][id], label='Optimum', color='red', linestyle='dashed')
    axs[0][0].grid()
    # P_loss vs w_channel
    axs[0][1].plot(w_channel_range*1e6, P_loss_channel)
    axs[0][1].set_xlabel("Channel width $w_c$ [$\\mu$m]")
    #axs[0][1].set_ylabel("$P_{{loss}}$ [W]")
    axs[0][1].axvline(data['w_channel'][id]*1e6, label='Optimum', color='red', linestyle='dashed')
    axs[0][1].grid()
    # Total power consumption
    axs[1][0].plot(w_channel_spacing_range*1e6, P_loss_channel_spacing)
    axs[1][0].set_xlabel("Channel spacing $s_c$ [$\\mu$m]")
    axs[1][0].set_ylabel("$P_{{loss}}$ [W]")
    axs[1][0].axvline(data['w_channel_spacing'][id]*1e6, label='Optimum', color='red', linestyle='dashed')
    axs[1][0].grid()
    # Total power consumption
    axs[1][1].plot(T_wall_range, P_loss_wall)
    axs[1][1].set_xlabel("Top wall temperature $T_{{wall,top}}$ [K]")
    #axs[1][1].set_ylabel("$P_{{loss}}$ [W]")
    axs[1][1].axvline(data['T_wall'][id], label='Optimum', color='red', linestyle='dashed')
    axs[1][1].grid()
    fig.suptitle("Power loss around found optimum "+ str_title)

    fig.tight_layout()
    
def alternative_power_loss(F_desired,T_chamber, channel_amount, settings, w_channel, w_channel_spacing, T_wall):
    """ Copy of the algorithem, but exploring other options

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
    # print("\n  ----- OPTIMIZING THRUSTER FOR POWER CONSUMPTION -----")
    # print("\n    HIGH-LEVEL INPUTS")
    # print("--- Desired thrust:        {:3.1f} mN".format(F_desired*1e3))
    # print("--- Chamber temperature:   {:4.1f} K".format(T_chamber))
    
    T_chamber=float(T_chamber)

    ## First determine ideal engine performance, and load settings for that
    T_inlet = settings['FDP']['T_inlet'] # [K] Inlet temperature (at heating channel inlet manifold)
    p_inlet = settings['FDP']['p_inlet'] # [Pa] Chamber pressure assumed equal to inlet pressure
    AR_exit = settings['FDP']['AR_exit'] # [-] Exit area ratio of nozzle
    p_back = settings['FDP']['p_back'] # [Pa]
    assert p_back == 0 # Nozzle correction does not work for back pressure
    fp = settings['FDP']['fp'] # [object] FluidProperties object containing thermodynamic parameters for propellant

    # print("--- Inlet pressure:          {:2.1f} bar".format(p_inlet*1e-5))
    # print("--- Exit area ratio          {:2.1f}".format(AR_exit))

    


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

    # print("\n    IDEAL PERFORMANCE")
    # print("--- Nozzle status:           {}".format(ep_ideal['nozzle_status']))
    # print("--- Mass flow:               {:3.2f} mg/s".format(m_dot*1e6))
    # print("--- Isp:                     {:3.1f} s".format(Isp))
    # print("--- Power consumption:       {:2.1f} W".format(P_ideal))



    # Calculate pre-calculated values of heating channel that are constant through entire optimization
    prepared_values = oneD.full_homogenous_preparation(
        T_inlet=T_inlet,                            # [K] Inlet temperature of the channel
        T_outlet=T_chamber,                         # [K] Outlet of the channel is equal to chamber temperature in IRT
        m_dot=m_dot/channel_amount,                 # [kg/s] Mass flow in a single channel
        p_ref=p_inlet,                              # [Pa] The pressure in the channel is assumed to be constant for heat flow calculation purposes
        steps_l=settings['steps']['steps_l'],       # [-] Fidelity of simulation in liquid section of channel
        steps_tp=settings['steps']['steps_tp'],     # [-] Fidelity of simulation in two-phase section of channel
        steps_g=settings['steps']['steps_g'],       # [-] Fidelity of simulation in gas section of channel
        fp=fp,                                      # [object] Object to access thermodynamic parameters of propellant with
    )
    # Result for alternative input
    x_alt = [
        w_channel,
        w_channel_spacing,
        T_wall
        ]
    
    # The function to optimize
    def f(x, return_full_results=False): #  NOTE: change bounds above if argument change
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
            wall_args=wall_args,                                                            # {dict} Arguments about wall design. See aboves                          
        )
        # Return the total power consumption. The variable of interest
        P_total = P_ideal + res_P['P_loss']
        #print(x)

        if math.isnan(P_total):
            #print("A constraint violated.")
            # Punish the objective function by outputting high power consumptions for invalid solutions
            # The punishment must however not be constant, or it will converge at the boundary where it is invalid.
            # Since the punishment stems from the pressure drop being higher than the initial p_chamber value, the lower the negative final pressure is the higher the punishment
            return (P_ideal*(2+1*res_P['pressure_drop_punishment']))*settings['function_scaling']
        else:
            #print("P_total: {:2.5f} W".format(P_total))
            # Scipy minimize can't handle multiple returns for the objective function, so the extensive final results must be pulled out afterwards
            if return_full_results:
                # print(res_P)
                return res_P
            else:
                return P_total*settings['function_scaling']



    res_final = f(x_alt, return_full_results=True)

    optim_results = {
        'minimize_results': None,
        'w_channel': None,
        'w_channel_spacing': None,
        'T_wall_top': None,
        'P_total': P_ideal + res_final['P_loss'],
        'P_loss': res_final['P_loss'],
        'P_ideal': P_ideal,
        'l_channel': res_final['l_channel'],
        'T_wall_bottom': res_final['T_wall_bottom'],
        'w_total': res_final['w_total'],
        'pressure_drop': res_final['pressure_drop'],
        'Re_channel_l': res_final['Re_channel_l'],
        'Re_channel_tp': res_final['Re_channel_tp'],
        'Re_channel_g': res_final['Re_channel_g'],
        'M_channel_exit_after_dP': res_final['M_channel_exit_after_dP'],
        'Re_channel_exit_after_dP': res_final['Re_channel_exit_after_dP'],
        'm_dot': ep_ideal['m_dot'],
        'Isp': ep_ideal['Isp'],
        'p_inlet': p_inlet,
        'w_throat_new': res_final['w_throat_new'],
        'Re_throat_new': res_final['Re_throat_new'],
    }
    return optim_results
    

if __name__ == "__main__":
    run()