# File to automatically run optimization for various temperatures
import main_1D_rectangular_multichannel_optimization as optim
import math
import numpy as np

# Default values in case it is not run from command line
F = 5e-3
T_l= 500
T_h = 575
s = round((T_h-T_l)/25)+1

str_f = "optimization_results-5mN/"

def run(F_desired=F, T_low=T_l, T_high=T_h, steps=s, str_folder=str_f):
    print(s)
    T_range = np.linspace(start=T_low, stop=T_high, num=steps)
    print(T_range)
    T_iter = np.nditer(T_range)
    for T in T_iter:
        optim.run(F_desired=F_desired, T_chamber=round(float(T)), str_folder=str_folder)

if __name__ == "__main__":
    run()