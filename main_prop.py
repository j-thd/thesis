# Script to quickly get some fluid properties at specified temperature and pressure
import numpy as np
import matplotlib.pyplot as plt
from thermo.prop import FluidProperties

def run():
    printBasicStuff()
    # plotSaturationCurve()
    # plotGammaCurve()
    # plt.show()

def plotSaturationCurve():
    fp = FluidProperties("water")
    p_sat = np.linspace(start=0.01e5, stop=10.0e5, num=100) # [P] Pressure
    T_sat = fp.get_saturation_temperature(p_sat) # [K]

    plt.figure()
    plt.plot(T_sat, p_sat*1e-5)
    plt.xlabel("Saturation temperature - $T_{{sat}}$ [K]")
    plt.ylabel("Saturation pressure - $p_{{sat}}$ [bar]")
    plt.grid()
    plt.title("Saturation curve for water")
    plt.tight_layout()

def plotGammaCurve():
    fp = FluidProperties("water")
    temperature = np.linspace(start=500,stop=1000, num=200)
    pressure = np.linspace(start=1.0e5, stop=5.0e5, num=3)

    gamma = np.zeros((pressure.size,temperature.size))
    
    p_iter = np.nditer(pressure, flags=['c_index'])
    T_iter = np.nditer(temperature, flags= ['c_index'])

    for p in p_iter:
        i = p_iter.index
        T_iter.reset()
        for T in T_iter:
            j = T_iter.index
            gamma[i,j] = fp.get_specific_heat_ratio(T=float(T), p=float(p))
        

    plt.figure()
    p_iter.reset()
    for p in p_iter:
        i = p_iter.index
        plt.plot(temperature, gamma[i,:], label="{:1.0f} bar".format(p*1e-5))
    
    plt.xlabel("Temperature - $T$ [K]")
    plt.ylabel("Specific heat ratio - $\\gamma$ [-] ")
    plt.grid()
    plt.legend(title="Pressure")
    plt.title("Specific heat ratio of water")
    plt.tight_layout()
    


def printBasicStuff():
    fp = FluidProperties("water")
    print("R: {:3.2f} J/(kg*K)".format(fp.get_specific_gas_constant()))

    

if __name__ == "__main__":
    run()