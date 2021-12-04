# File to automatically run optimization for various temperatures
import main_1D_rectangular_multichannel_optimization as optim
import sys
import numpy as np



# Default values in case it is not run from command line
if len(sys.argv)>1:
    F_input = float(sys.argv[1])*1e-3  # [N] Thrust is passed as micronewtons
    T_low = int(sys.argv[2]) # [K] Low temp
    #T_high = int(sys.argv[3]) # [K] High temp

F = F_input
T_l= T_low
T_h = 1100
s = round((T_h-T_l)/25)+1

#str_f = "optimization_results_SA_AR5/optimization_results-{:1.0f}mN/".format(F*1e3)
#str_f = "optimization_results_SA_AR15/optimization_results-{:1.0f}mN/".format(F*1e3)
str_f = "optimization_results_SA_DIV30/optimization_results-{:1.0f}mN/".format(F*1e3)
#str_f = "optimization_results_SA_w_channel_10/optimization_results-{:1.0f}mN/".format(F*1e3)
#str_f = "optimization_results_SA_w_channel_spacing_25/optimization_results-{:1.0f}mN/".format(F*1e3)
#str_f = "optimization_results-{:1.0f}mN/".format(F*1e3)

def run(F_desired=F, T_low=T_l, T_high=T_h, steps=s, str_folder=str_f):
    print(s)
    T_range = np.linspace(start=T_low, stop=T_high, num=steps)
    print(T_range)
    T_iter = np.nditer(T_range)
    for T in T_iter:
        optim.run(F_desired=F_desired, T_chamber=round(float(T)), str_folder=str_folder)

if __name__ == "__main__":
    run()