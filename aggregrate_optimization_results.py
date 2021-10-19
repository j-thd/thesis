import os
import numpy as np
import matplotlib.pyplot as plt


def run():
    # Strings to point towards results for different thrust levels
    str_1mN_folder = "optimization_results-1mN/"
    str_2mN_folder = "optimization_results-2mN/"
    str_3mN_folder = "optimization_results-3mN/"
    str_4mN_folder = "optimization_results-4mN/"
    str_5mN_folder = "optimization_results-5mN/"
    # Load data for different thrust levels
    d1 = load_data_in_folder(str_1mN_folder)
    d2 = load_data_in_folder(str_2mN_folder)
    d3 = load_data_in_folder(str_3mN_folder)
    d4 = load_data_in_folder(str_4mN_folder)
    d5 = load_data_in_folder(str_5mN_folder)
    # Make into a data list, to iterate over when plotting results
    dl = [d1, d2, d3, d4, d5]

    # Some plots to highlight lose results for one thrust level
    dh=d4
    #plotIspVsPower(dh)
    #plotOptimalDesign(dh)
    #plotThroatPressureResults(dh)
    #plotHighLevelStuff(dh)
    #plotMassFlowAndIsp(dh)
    # Some plots showing all thrust levels together
    #plotChannelNumbersVsIsp(dl)
    #plotTopWallSuperheatVsIsp(dl)
    plotThroatWidthVsIsp(dl)
    plotPressureDropVsIsp(dl)
    #plotTotalPowerVsIsp(dl=dl)
    #plotPowerLossVsIsp(dl=dl)
    # plotHeatingEfficiencyVsIsp(dl=dl)
    # plotChipAreaVsIsp(dl=dl)

    plt.show()

def load_data_in_folder(str_folder):
    npz_files = discover_npz_files(str_folder) # Finds all npz files in the folder
    npz_data = read_and_order_npz_data(npz_files) # Files are ordered by chamber temperature
    data = process_data(npz_data) # From each file only the optimum design must be extracted for each temperature

    return data # Return the data of each optimum associated with a thrust and chamber temperature

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
        #print(f)
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
    P_loss = np.zeros_like(T_chamber)
    m_dot = np.zeros_like(T_chamber)

    # Optimal design
    channel_amount = np.zeros_like(T_chamber)
    w_channel = np.zeros_like(T_chamber)
    w_channel_spacing = np.zeros_like(T_chamber)
    T_wall = np.zeros_like(T_chamber)

    # Throat/pressure drop characteristics
    pressure_drop = np.zeros_like(T_chamber)
    w_throat = np.zeros_like(T_chamber)
    w_nozzle = np.zeros_like(T_chamber)
    w_inlet_manifold = np.zeros_like(T_chamber)
    l_channel = np.zeros_like(T_chamber)
    l_total = np.zeros_like(T_chamber)
    A_chip = np.zeros_like(T_chamber)

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
        P_loss[i] = d['P_loss'][id]
        m_dot[i] = d['m_dot'][id]
        # Store optimal design input
        channel_amount[i] = d['channel_amount_range'][id]
        w_channel[i] = d['w_channel'][id]
        w_channel_spacing[i] = d['w_channel_spacing'][id]
        T_wall[i] = d['T_wall'][id]
        # Throat/pressure drop characteristcs
        pressure_drop[i] = d['pressure_drop'][id]
        w_throat[i] = d['w_throat_new'][id]
        w_nozzle[i] = d['w_nozzle'][id]
        w_inlet_manifold[i] = d['w_inlet'][id]
        l_channel[i] = d['l_channel'][id]
        l_total[i] = d['l_total'][id]
        A_chip[i] = d['A_chip'][id]

    return {
        'F_desired': F_desired,
        'p_inlet': p_inlet,
        'T_chamber': T_chamber,
        'Isp': Isp,
        'P_ideal': P_ideal,
        'P_total': P_total,
        'P_loss': P_loss,
        'm_dot': m_dot,
        'channel_amount': channel_amount,
        'T_wall': T_wall,
        'w_channel': w_channel,
        'w_channel_spacing': w_channel_spacing,
        'pressure_drop': pressure_drop,
        'w_throat': w_throat,
        'w_nozzle': w_nozzle,
        'w_inlet_manifold': w_inlet_manifold,
        'l_channel': l_channel,
        'l_total': l_total,
        'A_chip': A_chip,
    }

def plotTotalPowerVsIsp(dl):
    plt.figure()
    for data in dl:
        plt.plot(data['Isp'], data['P_total'], label="{:1.0f} mN".format(data['F_desired'][0]*1e3))
    plt.xlabel("Specific impulse - $I_{{sp}}$ [s]")
    plt.ylabel("Total power consumption - $P_t$ [W]")
    plt.grid()
    plt.legend(title="Thrust")
    plt.title("Total power consumption for given thrust and specific impulse")
    plt.tight_layout()

def plotThroatWidthVsIsp(dl):
    plt.figure()
    for data in dl:
        plt.plot(data['Isp'], data['w_throat']*1e6, label="{:1.0f} mN".format(data['F_desired'][0]*1e3))
        #plt.plot(data['Isp'], data['w_inlet_manifold']*1e6,label="{:1.0f} mN".format(data['F_desired'][0]*1e3) , linestyle='dashed')
    plt.xlabel("Specific impulse - $I_{{sp}}$ [s]")
    plt.ylabel("Throat width - $w_t$ [$\\mu$m]")
    plt.grid()
    plt.legend(title="Thrust")
    plt.title("Throat width for given thrust and specific impulse")
    plt.tight_layout()

def plotPressureDropVsIsp(dl):
    plt.figure()
    for data in dl:
        plt.plot(data['Isp'], data['pressure_drop']*1e-5, label="{:1.0f} mN".format(data['F_desired'][0]*1e3))
    plt.xlabel("Specific impulse - $I_{{sp}}$ [s]")
    plt.ylabel("Pressure drop - $\Delta p$ [bar]")
    plt.grid()
    plt.legend(title="Thrust")
    plt.title("Pressure drop for given thrust and specific impulse")
    plt.tight_layout()

def plotPowerLossVsIsp(dl):
    plt.figure()
    for data in dl:
        plt.plot(data['Isp'], data['P_loss'], label="{:1.0f} mN".format(data['F_desired'][0]*1e3))
    plt.xlabel("Specific impulse - $I_{{sp}}$ [s]")
    plt.ylabel("Power Loss - $P_t$ [W]")
    plt.grid()
    plt.legend(title="Thrust")
    plt.title("Power loss for given thrust and specific impulse")
    plt.tight_layout()

def plotChannelNumbersVsIsp(dl):
    plt.figure()
    for data in dl:
        plt.plot(data['Isp'], data['channel_amount'], label="{:1.0f} mN".format(data['F_desired'][0]*1e3))
    plt.xlabel("Specific impulse - $I_{{sp}}$ [s]")
    plt.ylabel("Number of channels - $N_c$ [-]")
    plt.grid()
    plt.legend(title="Thrust")
    plt.title("Optimal number of channels for given thrust and specific impulse")
    plt.tight_layout()

def plotTopWallSuperheatVsIsp(dl):
    plt.figure()
    for data in dl:
        plt.plot(data['Isp'], data['T_wall']-data['T_chamber'], label="{:1.0f} mN".format(data['F_desired'][0]*1e3))
    plt.xlabel("Specific impulse - $I_{{sp}}$ [s]")
    plt.ylabel("Top wall superheat ($T_{{w}} -T_c$) [K]")
    plt.grid()
    plt.legend(title="Thrust")
    plt.title("Optimal top wall superheat for given thrust and specific impulse")
    plt.tight_layout()

def plotHeatingEfficiencyVsIsp(dl):
    plt.figure()
    for data in dl:
        plt.plot(data['Isp'], data['P_ideal']/data['P_total'], label="{:1.0f} mN".format(data['F_desired'][0]*1e3))
    plt.xlabel("Specific impulse - $I_{{sp}}$ [s]")
    plt.ylabel("Heating efficiency - $\\mu$ [-]")
    plt.grid()
    plt.legend(title="Thrust")
    plt.title("Heating efficiency for given thrust and specific impulse")
    plt.tight_layout()

def plotChipAreaVsIsp(dl):
    plt.figure()
    for data in dl:
        plt.plot(data['Isp'], data['A_chip']*1e6, label="{:1.0f} mN".format(data['F_desired'][0]*1e3))
    plt.xlabel("Specific impulse - $I_{{sp}}$ [s]")
    plt.ylabel("Chip area- $A_{{chip}}$ [mm$^2$]")
    plt.grid()
    plt.legend(title="Thrust")
    plt.title("Chip area for given thrust and specific impulse")
    plt.tight_layout()

def plotIspVsPower(data):
    plt.figure()

    plt.plot(data['T_chamber'], data['P_total'], label="Total -$P_{{total}}$")
    plt.plot(data['T_chamber'], data['P_ideal'], label="Ideal - $P_{{\Delta h}}$")
    plt.ylabel("Power consumption - P [W]")
    plt.xlabel("Chamber temperature - $T_c$ [K]")
    plt.title("Optimal power consumption for {:2.0f} mN".format(data['F_desired'][0]*1e3))
    plt.legend()
    plt.grid()
    plt.tight_layout()

def plotMassFlowAndIsp(data):
    fig, ax1 = plt.subplots()
    ax2 = plt.twinx()
    ax1.plot(data['T_chamber'], data['m_dot']*1e6, c='tab:red', linestyle='dashed', label="Mass flow")
    ax2.plot(data['T_chamber'], data['Isp'], c='tab:blue', label="$I_{{sp}}$")
    ax1.set_ylabel("Mass flow - $\\dot{{m}}$ [mg$\\cdot$s$^{{-1}}$]")
    ax1.set_xlabel("Chamber temperature - $T_c$ [K]")
    ax2.set_ylabel("Specific impulse - $I_{{sp}}$ [s]")
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines+lines2, labels+labels2)
    fig.suptitle("Mass flow and chamber temperature for {:2.0f} mN".format(data['F_desired'][0]*1e3))
    ax1.grid()
    fig.tight_layout()

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