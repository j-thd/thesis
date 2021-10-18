import os
import math
import numpy as np
import matplotlib.pyplot as plt


def run():
    str_4mN_folder = "optimization_results-4mN/"
    str_5mN_folder = "optimization_results-5mN/"
    npz_files = discover_npz_files(str_5mN_folder)
    npz_data = read_and_order_npz_data(npz_files)
    data = process_data(npz_data)

    plotIspVsPower(data)
    plotOptimalDesign(data)
    plotThroatPressureResults(data)
    plotHighLevelStuff(data)

    plt.show()

def discover_npz_files(str_folder):
    # Discover and return list of all npz files
    npz_files = [] # Store all found npz files in here
    dirs = os.listdir(str_folder)
    for f in dirs:
        if f.endswith(".npz"):
            npz_files.append(str_folder + f)
    
    return npz_files

def read_and_order_npz_data(npz_files):
    npz_data = [] # Store all data sets in here
    T_chamber = [] # Read and store T_chamber in here, so npz_data can be sorted based on that

    for f in npz_files:
        npz_file = open(f, "rb")
        nd = np.load(npz_file)
        npz_data.append(nd)
        T_chamber.append(float(nd['T_chamber']))

    npz_data = [ x for _,x in sorted(zip(T_chamber,npz_data))]
    return npz_data

def process_data(npz_data):
    # Count data points and construct numpy arrays to put results in
    data_points = len(npz_data)
    

    
    # High-level performance
    T_chamber = np.zeros(data_points)
    Isp = np.zeros_like(T_chamber)
    F_desired = np.zeros_like(T_chamber)
    p_inlet = np.zeros_like(T_chamber)
    P_ideal = np.zeros_like(T_chamber)
    P_total = np.zeros_like(T_chamber)
    m_dot = np.zeros_like(T_chamber)

    # Optimal design
    channel_amount = np.zeros_like(T_chamber)
    w_channel = np.zeros_like(T_chamber)
    w_channel_spacing = np.zeros_like(T_chamber)
    T_wall = np.zeros_like(T_chamber)

    # Throat/pressure drop characteristics
    pressure_drop = np.zeros_like(T_chamber)
    w_throat = np.zeros_like(T_chamber)
    l_channel = np.zeros_like(T_chamber)
    l_total = np.zeros_like(T_chamber)

    T_iter = np.nditer(T_chamber, flags=['c_index'])

    for T in T_iter:
        i = T_iter.index # Index
        d = npz_data[i] # Data set to process this iteration
        # Find id of result with minimal power consumption in data point
        id = np.argmin(d['P_total'])
        # Store high-level values
        T_chamber[i] = d['T_chamber']
        F_desired[i] = d['F_desired']
        Isp[i] = d['Isp'][id]
        P_ideal[i] = d['P_ideal'][id]
        P_total[i] = d['P_total'][id]
        m_dot[i] = d['m_dot'][id]
        # Store optimal design input
        channel_amount[i] = d['channel_amount_range'][id]
        w_channel[i] = d['w_channel'][id]
        w_channel_spacing[i] = d['w_channel_spacing'][id]
        T_wall[i] = d['T_wall'][id]
        # Throat/pressure drop characteristcs
        pressure_drop[i] = d['pressure_drop'][id]
        w_throat[i] = d['w_throat_new'][id]
        l_channel[i] = d['l_channel'][id]
        l_total[i] = d['l_total'][id]

    return {
        'F_desired': F_desired,
        'p_inlet': p_inlet,
        'T_chamber': T_chamber,
        'Isp': Isp,
        'P_ideal': P_ideal,
        'P_total': P_total,
        'm_dot': m_dot,
        'channel_amount': channel_amount,
        'T_wall': T_wall,
        'w_channel': w_channel,
        'w_channel_spacing': w_channel_spacing,
        'pressure_drop': pressure_drop,
        'w_throat': w_throat,
        'l_channel': l_channel,
        'l_total': l_total,
    }

def plotIspVsPower(data):
    plt.figure()

    plt.plot(data['Isp'], data['P_total'], label="Total -$P_{{total}}$")
    plt.plot(data['Isp'], data['P_ideal'], label="Ideal - $P_{{\Delta h}}$")
    plt.ylabel("Power consumption - P [W]")
    plt.xlabel("Specific impulse - $I_{{sp}}$ [s]")
    plt.title("Optimal power consumption for {:2.1f} mN".format(data['F_desired'][0]*1e3))
    plt.legend()
    plt.grid()

def plotHighLevelStuff(data):
    fig, axs = plt.subplots(2,2)
    # Power consumption vs. chamber temperature
    axs[0][0].plot(data['T_chamber'], data['P_total'])
    axs[0][0].set_ylabel("Total Power Consumption - $P_{{total}}$ [W]")
    #axs[0][0].set_xlabel("Chamber temperature - $T_c$ [K]")
    axs[0][0].grid()
    axs[0][1].plot(data['T_chamber'], data['m_dot']*1e6)
    axs[0][1].set_ylabel("Mass flow - $\\dot{{m}}$ [mg$\\cdot$s$^{{-1}}$]")
    axs[0][1].grid()
    axs[1][0].plot(data['T_chamber'], data['Isp'])
    axs[1][0].set_ylabel("Specific Impulse - $I_{{sp}}$ [s]")
    axs[1][0].set_xlabel("Chamber temperature - $T_c$ [K]")
    axs[1][0].grid()
    axs[1][1].plot(data['T_chamber'], (data['P_ideal']/data['P_total']))
    axs[1][1].set_ylabel("Heating efficiency - $\\mu$ [-]")
    axs[1][1].set_xlabel("Chamber temperature - $T_c$ [K]")
    axs[1][1].grid()
    fig.suptitle("High-level performance ({:2.1f}) mN".format(data['F_desired'][0]*1e3))
    #axs[0][0].legend()
    plt.tight_layout(pad=0.5)

def plotOptimalDesign(data):
    fig, axs = plt.subplots(2,2)
    # Total power consumption
    axs[0][0].plot(data['Isp'], data['channel_amount'])
    axs[0][0].set_ylabel("Number of channels - $N_c$ [-]")
    axs[0][0].grid()
    axs[0][1].plot(data['Isp'], data['T_wall']-data['T_chamber'])
    axs[0][1].set_ylabel("Top wall superheat - $(T_w-T_c)$ [K]")
    axs[0][1].grid()
    axs[1][0].plot(data['Isp'], data['w_channel']*1e6)
    axs[1][0].set_ylabel("Channel width - $w_c$ [$\\mu$m]")
    axs[1][0].set_xlabel("Specific Impulse - $I_{{sp}}$ [s]")
    axs[1][0].grid()
    axs[1][1].plot(data['Isp'], data['w_channel_spacing']*1e6)
    axs[1][1].set_ylabel("Channel spacing - $s_c$ [$\\mu$m]")
    axs[1][1].set_xlabel("Specific Impulse - $I_{{sp}}$ [s]")
    axs[1][1].grid()
    fig.suptitle("Optimal design for {:2.1f} mN for given $I_{{sp}}$".format(data['F_desired'][0]*1e3))
    #axs[0][0].legend()
    plt.tight_layout(pad=0.5)

def plotThroatPressureResults(data):
    fig, axs = plt.subplots(2,2)
    # Total power consumption
    axs[0][0].plot(data['Isp'], data['pressure_drop']*1e-5)
    axs[0][0].set_ylabel("Pressure drop - $\Delta p$ [bar]")
    axs[0][0].grid()
    axs[0][1].plot(data['Isp'], data['w_throat']*1e6)
    axs[0][1].set_ylabel("Throat width - $w_t$ [$\\mu$m]")
    axs[0][1].grid()
    axs[1][0].plot(data['Isp'], data['l_channel']*1e3)
    axs[1][0].set_ylabel("Channel length - $l_c$ [mm]")
    axs[1][0].set_xlabel("Specific Impulse - $I_{{sp}}$ [s]")
    axs[1][0].grid()
    axs[1][1].plot(data['Isp'], data['l_total']*1e3)
    axs[1][1].set_ylabel("Total chip length - $l_c$ [mm]")
    axs[1][1].set_xlabel("Specific Impulse - $I_{{sp}}$ [s]")
    axs[1][1].grid()
    fig.suptitle("???????? for {:2.1f} mN for given $I_{{sp}}$".format(data['F_desired'][0]*1e3))
    #axs[0][0].legend()
    plt.tight_layout(pad=0.5)


if __name__ == "__main__":
    run()