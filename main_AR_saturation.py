# File to plot what the minimum temperature is for each area ratio, before condensation occurs at the nozzle exit
# Calculation is done based on Ideal Rocket Theory, and saturation curve of propellant

import numpy as np
import matplotlib.pyplot as plt
import basic.IRT as IRT
from thermo.prop import FluidProperties

def plot_condensation_curve(AR_min, AR_max, p_c, fp, T_ref):
    gamma = fp.get_specific_heat_ratio(T=T_ref,p=p_c) # [-] Specific heat ratio

    # Determine maximum Mach number at exit
    M_exit_max = IRT.Mach_from_area_ratio(AR=AR_max, gamma=gamma) # [-] Max. Mach number at exit
    M_exit_min = IRT.Mach_from_area_ratio(AR=AR_min, gamma=gamma) # [-] Min. Mach number at exit
    M_exit = np.linspace(start=M_exit_min, stop=M_exit_max, num=50) # [-] Range of exit Mach numbers to evalulate

    # Determine pressure at exit
    PR_exit = IRT.pressure_ratio(M=M_exit, gamma=gamma) # [-] Pressure ratio at exit
    p_exit = p_c / PR_exit # [-] Exit pressure

    # Determine saturation temperature at exit
    T_sat_exit = fp.get_saturation_temperature(p=p_exit) # [K] Saturation temperature at exit
    # print(p_exit)
    # print(T_sat_exit)
    TR_exit = IRT.temperature_ratio(M=M_exit, gamma=gamma) # [-] Temperature ratio at exit
    T_chamber = T_sat_exit*TR_exit # [K] Chamber temperature that would result precisely in saturation temperature at nozzle exit

    AR = IRT.area_ratio(M=M_exit, gamma=gamma) # [-] Area ratios corresponding to all specified mach numbers

    plt.plot(AR, T_chamber, label ="{:2.0f} bar , $T_c={:3.0f}$ K, $\\gamma={:1.2f}$".format(p_c*1e-5,T_ref,gamma))

def plot_pressure_curve(AR_min, AR_max, p_c, fp, T_ref):
    gamma = fp.get_specific_heat_ratio(T=T_ref,p=p_c) # [-] Specific heat ratio

    # Determine maximum Mach number at exit
    M_exit_max = IRT.Mach_from_area_ratio(AR=AR_max, gamma=gamma) # [-] Max. Mach number at exit
    M_exit_min = IRT.Mach_from_area_ratio(AR=AR_min, gamma=gamma) # [-] Min. Mach number at exit
    M_exit = np.linspace(start=M_exit_min, stop=M_exit_max, num=50) # [-] Range of exit Mach numbers to evalulate

    # Determine pressure at exit
    PR_exit = IRT.pressure_ratio(M=M_exit, gamma=gamma) # [-] Pressure ratio at exit
    p_exit = p_c / PR_exit # [-] Exit pressure

    AR = IRT.area_ratio(M=M_exit, gamma=gamma) # [-] Area ratios corresponding to all specified mach numbers

    plt.plot(AR, p_exit*1e-5, label ="{:2.0f} bar , $\\gamma={:1.2f}$".format(p_c*1e-5,gamma))

def plot_TR_curve(AR_min, AR_max, p_c, fp, T_ref):
    gamma = fp.get_specific_heat_ratio(T=T_ref,p=p_c) # [-] Specific heat ratio

    # Determine maximum Mach number at exit
    M_exit_max = IRT.Mach_from_area_ratio(AR=AR_max, gamma=gamma) # [-] Max. Mach number at exit
    M_exit_min = IRT.Mach_from_area_ratio(AR=AR_min, gamma=gamma) # [-] Min. Mach number at exit
    M_exit = np.linspace(start=M_exit_min, stop=M_exit_max, num=50) # [-] Range of exit Mach numbers to evalulate

    # Determine pressure at exit
    #PR_exit = IRT.pressure_ratio(M=M_exit, gamma=gamma) # [-] Pressure ratio at exit
    #p_exit = p_c / PR_exit # [-] Exit pressure

    # Determine saturation temperature at exit
    #T_sat_exit = fp.get_saturation_temperature(p=p_exit) # [K] Saturation temperature at exit
    # print(p_exit)
    # print(T_sat_exit)
    TR_exit = IRT.temperature_ratio(M=M_exit, gamma=gamma) # [-] Temperature ratio at exit
    #T_chamber = T_sat_exit*TR_exit # [K] Chamber temperature that would result precisely in saturation temperature at nozzle exit

    AR = IRT.area_ratio(M=M_exit, gamma=gamma) # [-] Area ratios corresponding to all specified mach numbers

    plt.plot(AR, TR_exit, label ="{:2.0f} bar , $\\gamma={:1.2f}$".format(p_c*1e-5,gamma))


def run():
    AR_max = 100 # [-] Maximum exit area ratio to consider
    AR_min = 5 # [-] Minimum exit area ratio to consider
    fp = FluidProperties("HEOS::Water") # Object to retrieve fluid properties from
    p_c = 5e5 # [Pa] Chamber pressure is important because it determines the saturation pressure at the exit
    T_ref = 600 # [K] Arbitrary choice, must be varied to check effect on results
    plot_condensation_curve(AR_min=AR_min, AR_max=AR_max, p_c=p_c, fp=fp, T_ref=T_ref)
    T_ref = 800 # [K] Arbitrary choice, must be varied to check effect on results
    plot_condensation_curve(AR_min=AR_min, AR_max=AR_max, p_c=p_c, fp=fp, T_ref=T_ref)
    T_ref = 1000
    plot_condensation_curve(AR_min=AR_min, AR_max=AR_max, p_c=p_c, fp=fp, T_ref=T_ref)
    p_c = 10e5 # [Pa] Chamber pressure is important because it determines the saturation pressure at the exit
    T_ref = 600 # [K] Arbitrary choice, must be varied to check effect on results
    plot_condensation_curve(AR_min=AR_min, AR_max=AR_max, p_c=p_c, fp=fp, T_ref=T_ref)
    T_ref = 800
    plot_condensation_curve(AR_min=AR_min, AR_max=AR_max, p_c=p_c, fp=fp, T_ref=T_ref)
    T_ref = 1000
    plot_condensation_curve(AR_min=AR_min, AR_max=AR_max, p_c=p_c, fp=fp, T_ref=T_ref)
    plt.grid()
    plt.title("Minimum chamber temperature to prevent condensation")
    plt.xlabel("Exit area ratio $\\frac{{A_e}}{{A_t}}$ [-]")
    plt.ylabel("Chamber temperature $T_c$ [K]")
    plt.legend()
    

    fig = plt.figure()
    p_c = 5e5
    T_ref = 600
    plot_pressure_curve(AR_min=AR_min, AR_max=AR_max, p_c=p_c, fp=fp, T_ref=T_ref)
    T_ref = 800
    plot_pressure_curve(AR_min=AR_min, AR_max=AR_max, p_c=p_c, fp=fp, T_ref=T_ref)
    T_ref = 1000
    plot_pressure_curve(AR_min=AR_min, AR_max=AR_max, p_c=p_c, fp=fp, T_ref=T_ref)
    plt.grid()
    plt.title("Exit pressure with varying area ratio")
    plt.xlabel("Exit area ratio $\\frac{{A_e}}{{A_t}}$ [-]")
    plt.ylabel("Exit pressure $p_e$ [bar]")
    plt.legend()

    fig = plt.figure()
    p_c = 5e5
    T_ref = 600
    plot_TR_curve(AR_min=AR_min, AR_max=AR_max, p_c=p_c, fp=fp, T_ref=T_ref)
    T_ref = 800
    plot_TR_curve(AR_min=AR_min, AR_max=AR_max, p_c=p_c, fp=fp, T_ref=T_ref)
    T_ref = 1000
    plot_TR_curve(AR_min=AR_min, AR_max=AR_max, p_c=p_c, fp=fp, T_ref=T_ref)
    plt.grid()
    plt.title("Temperature ratio with varying area ratio")
    plt.xlabel("Exit area ratio $\\frac{{A_e}}{{A_t}}$ [-]")
    plt.ylabel("Temperature ratio [-]")
    plt.legend()
    
    plt.show()


if __name__ == "__main__":
    run()