# File to configure optimization and run it
import sys
import numpy as np
from matplotlib import pyplot as plt
import timeit # To time speed of optmization

import optimization
import optimization.settings




if __name__ == "__main__":
    F_desired = 4e-3   # [N] Desired thrust
    T_chamber = 1000     # [K] Chamber temperature
    if len(sys.argv)>1:
        F_desired = float(sys.argv[1])*1e-3  # [N] Thrust is passed as micronewtons
        T_chamber = int(sys.argv[2]) # [K] Chamber temperature

    save_file = open("optimization_results/optimization_results-F{:1.0f}mN-{:3.0f}K.npz".format(F_desired*1e3, T_chamber), "wb")

    # Load all settings
    channel_amount_range = np.arange(1,30)
    P_total = np.zeros(len(channel_amount_range))
    P_loss = np.zeros_like(P_total)
    w_channel = np.zeros_like(P_total)
    w_channel_spacing = np.zeros_like(P_total)
    w_channel = np.zeros_like(P_total)
    l_channel = np.zeros_like(P_total)
    w_total = np.zeros_like(P_total)
    Re_channel_l = np.zeros_like(P_total)
    Re_channel_tp = np.zeros_like(P_total)
    Re_channel_g = np.zeros_like(P_total)
    M_channel_exit_after_dP = np.zeros_like(P_total)
    Re_channel_exit_after_dP = np.zeros_like(P_total)
    pressure_drop = np.zeros_like(P_total)
    m_dot = np.zeros_like(P_total)
    Isp = np.zeros_like(P_total)
    p_inlet = np.zeros_like(P_total)
    T_wall = np.zeros_like(P_total)
    T_wall_bottom = np.zeros_like(P_total)
    w_throat_new = np.zeros_like(P_total)
    Re_throat_new = np.zeros_like(P_total)
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
        P_loss[i_channel.index] = results['P_loss']
        w_channel[i_channel.index] = results['w_channel']
        l_channel[i_channel.index] = results['l_channel']
        w_total[i_channel.index] = results['w_total']
        pressure_drop[i_channel.index] = results['pressure_drop']
        Re_channel_l[i_channel.index] = results['Re_channel_l']
        Re_channel_tp[i_channel.index] = results['Re_channel_tp']
        Re_channel_g[i_channel.index] = results['Re_channel_g']
        Re_channel_exit_after_dP[i_channel.index] = results['Re_channel_exit_after_dP']
        M_channel_exit_after_dP[i_channel.index] = results['M_channel_exit_after_dP']
        w_channel_spacing[i_channel.index] = results['w_channel_spacing']
        T_wall[i_channel.index] = results['T_wall_top']
        T_wall_bottom[i_channel.index] = results['T_wall_bottom']
        m_dot[i_channel.index] = results['m_dot']
        Isp[i_channel.index] = results['Isp']
        p_inlet[i_channel.index] = results['p_inlet']
        w_throat_new[i_channel.index] = results['w_throat_new']
        Re_throat_new[i_channel.index] = results['Re_throat_new']
        
    stop = timeit.default_timer()
    print("Time elapsed: {} seconds".format(stop-start))

    # Save results
    np.savez(save_file,
        F_desired=F_desired,
        T_chamber=T_chamber,
        channel_amount_range=channel_amount_range,
        P_total=P_total,
        P_loss=P_loss,
        w_channel=w_channel,
        w_channel_spacing=w_channel_spacing,
        T_wall=T_wall,
        T_wall_bottom=T_wall_bottom,
        l_channel=l_channel,
        w_total=w_total,
        pressure_drop=pressure_drop,
        Re_channel_l=Re_channel_l,
        Re_channel_tp=Re_channel_tp,
        Re_channel_g=Re_channel_g,
        Re_channel_exit_after_dP=Re_channel_exit_after_dP,
        M_channel_exit_after_dP=M_channel_exit_after_dP,
        m_dot=m_dot,
        Isp=Isp,
        p_inlet=p_inlet,
        w_throat_new=w_throat_new,
        Re_throat_new=Re_throat_new
        )

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
    plt.plot(channel_amount_range, l_channel*1e3)
    plt.title("Channel length")
    plt.grid()

    plt.figure()
    plt.plot(channel_amount_range, w_channel_spacing*1e6)
    plt.title("Channel spacing")
    plt.grid()
    

    plt.figure()
    plt.plot(channel_amount_range, M_channel_exit_after_dP)
    plt.title("Mach")
    plt.grid()
 

    plt.figure()
    plt.plot(channel_amount_range, Re_channel_exit_after_dP)
    plt.title("Reynolds")
    plt.grid()

    plt.show()