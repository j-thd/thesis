# File to configure optimization and run it

import numpy as np
from matplotlib import pyplot as plt
import timeit # To time speed of optmization

import optimization
import optimization.settings




if __name__ == "__main__":
    F_desired = 4e-3   # [N] Desired thrust
    T_chamber = 750     # [K] Chamber temperature

    # Load all settings
    channel_amount_range = np.arange(1,50)
    P_total = np.zeros(len(channel_amount_range))
    w_channel = np.zeros_like(P_total)
    w_channel_spacing = np.zeros_like(P_total)
    w_channel = np.zeros_like(P_total)
    l_channel = np.zeros_like(P_total)
    Re_channel_l = np.zeros_like(P_total)
    Re_channel_tp = np.zeros_like(P_total)
    Re_channel_g = np.zeros_like(P_total)
    M_channel_exit_after_dP = np.zeros_like(P_total)
    Re_channel_exit_after_dP = np.zeros_like(P_total)
    pressure_drop = np.zeros_like(P_total)

    T_wall = np.zeros_like(P_total)
    T_wall_bottom = np.zeros_like(P_total)
    start = timeit.default_timer()

    # To speed up the optimization, start at the previous result, to save time.
    x_previous = None # The first iteration does a nomral search

    i_channel = np.nditer(channel_amount_range, flags=['c_index'])
    for i in i_channel:
        
        settings = optimization.settings.settings_1D_rectangular_multichannel
        
        results = optimization.run(
            F_desired=F_desired,
            T_chamber=T_chamber,
            channel_amount=i,
            settings=settings,
            x_guess = x_previous)

        # Store previous guess
        x_previous = results['minimize_results'].x
        # Report results
        print("\n--- Optimized for channel amount: {}".format(i))
        print("- Power consumption      {:2.3f} W".format(results['P_total']))
        print("- Channel width:         {:3.3f} micron".format(results['w_channel']*1e6))
        print("- Channel spacing:       {:3.3f} micron".format(results['w_channel_spacing']*1e6))
        print("- Top wall temperature:  {:3.2f} K".format(results['T_wall_top']))

        P_total[i_channel.index] = results['P_total']
        w_channel[i_channel.index] = results['w_channel']
        l_channel[i_channel.index] = results['l_channel']
        pressure_drop[i_channel.index] = results['pressure_drop']
        Re_channel_l[i_channel.index] = results['Re_channel_l']
        Re_channel_tp[i_channel.index] = results['Re_channel_tp']
        Re_channel_g[i_channel.index] = results['Re_channel_g']
        Re_channel_exit_after_dP[i_channel.index] = results['Re_channel_exit_after_dP']
        M_channel_exit_after_dP[i_channel.index] = results['M_channel_exit_after_dP']
        w_channel_spacing[i_channel.index] = results['w_channel_spacing']
        T_wall[i_channel.index] = results['T_wall_top']
        T_wall_bottom[i_channel.index] = results['T_wall_bottom']
        
    stop = timeit.default_timer()
    
    print("Time elapsed: {} seconds".format(stop-start))

    plt.figure()
    plt.plot(channel_amount_range, P_total)
    plt.title("Total power consumption")
    plt.grid()

    plt.figure()
    plt.plot(channel_amount_range, pressure_drop*1e-5)
    plt.title("pressure_drop")
    plt.grid()

    plt.figure()
    plt.plot(channel_amount_range, T_wall_bottom)
    plt.title("Bottom wall temperature")
    plt.grid()

    plt.figure()
    plt.plot(channel_amount_range, T_wall)
    plt.title("Top wall temperature")
    plt.grid()

    plt.figure()
    plt.plot(channel_amount_range, w_channel*1e6)
    plt.title("Channel width")
    plt.grid()

    plt.figure()
    plt.plot(channel_amount_range, l_channel*1e6)
    plt.title("Channel length")
    plt.grid()

    plt.figure()
    plt.plot(channel_amount_range, w_channel_spacing*1e6)
    plt.title("Channel spacing")
    plt.grid()
    plt.show()

    plt.figure()
    plt.plot(channel_amount_range, M_channel_exit_after_dP)
    plt.title("Mach")
    plt.grid()
    plt.show()

    plt.figure()
    plt.plot(channel_amount_range, Re_channel_exit_after_dP)
    plt.title("Reynolds")
    plt.grid()
    plt.show()