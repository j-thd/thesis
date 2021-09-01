''' Early experiment to generate the first results when combining the 1D- heat transfer model with some heat loss models
'''

import math

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.style as style  

import models.one_D as oneD
import thrusters.thruster_data
import thermo.convection
import thermo.two_phase as tp
import basic.chamber as chamber
from thermo.prop import FluidProperties
from basic.IRT import engine_performance_from_F_and_T
import physical_constants

# Set plotting style to colourblind-friendly style
style.use('tableau-colorblind10')

def calc_and_plot_channel_width(w_channel, ax_P_loss, fig_P_loss, ax_w_total, fig_w_total, ax_l_total, fig_l_total, ax_pressure_drop, fig_pressure_drop, ax_is_choked, fig_is_choked,\
    ax_Re_channel_exit, fig_Re_channel_exit):
    # Taking some thruster data to make the optmization a bit similar to existing cases
    td = None#thrusters.thruster_data.Cen2010_6

    # Desired chamber outlet/nozzle inlet conditions
    # Switch between an existing thruster or a desired one
    if td is not None: 
        m_dot = td['m_dot'] # [kg/s] Mass flow
        T_chamber = td['T_chamber_guess'] # [K] Chamber temperature
        T_inlet = 300 # [K] Inlet temperature
        p_inlet = td['p_inlet'] # [Pa] Inlet pressure, assumed to be constant through-out the channel
        F_desired = 7.71e-3 # td['F'] # [N] Thrust
        p_back = td['p_back']
        fp = FluidProperties(td['propellant']) # Object to access fluid properties
        AR_exit = td['AR_exit'] # [-] Exit area of nozzle 
        w_throat = td['w_throat'] # [m] Throat width
        #w_channel = 1.7*td['w_channel'] #  [m] Channel width
        inlet_manifold_length_factor = 2 # [m] Multiplication factor with inlet manifold width to determine manifold length
        inlet_manifold_width_factor = 5.5 # [-] Multiplication factor (with channel width to determine margin in chamber)
        l_exit_manifold = 1e-3 # [m] Length between the end of multiple channels and start of convergent nozzle
        convergent_half_angle  = td['convergent_half_angle'] # [rad] Half-angle of the convergent part of the nozzle
        divergent_half_angle = td['divergent_half_angle'] # [rad] Half-angle of divergent part of the nozzle
        w_channel_spacing = td['w_channel_spacing'] # [m] Spacing between channels (wall-to-wall)
        w_outer_margin = 2e-3 # [m] Margin around the outer channels for structural integrity
        channel_amount = td['channel_amount'] # [-] Number of channels
        h_channel = td['h_channel'] # [m] Channel depth
        emissivity_chip_top = 1 # [-] Assumed emissivity of chip at top-side


    else:
        F_desired = 10e-3 # [N] Desired thrust
        p_inlet = 5.0e5 # [Pa] Inlet pressure
        T_inlet = 300 # [K] Inlet temperature
        T_chamber = 500 # [K] Chamber temperature
        p_back = 0 # [Pa] Back pressure
        AR_exit = 10 # [-] Exit area ratio
        fp = FluidProperties("water") # (obj) Object to access water properties
        h_channel = 100e-6 # [m] Channel depth
        #w_channel = 25e-6 #  [m] Channel width

        # Find the mass flow and throat area IRT requires
        ep_ideal = engine_performance_from_F_and_T(F_desired=F_desired, p_chamber=p_inlet, T_chamber=T_chamber, AR_exit=AR_exit, p_back=p_back, fp=fp)
        m_dot = ep_ideal['m_dot'] # [kg/s] Mass flow
        w_throat = ep_ideal['A_throat']/h_channel # [m] Throat width
        Isp = ep_ideal['thrust']/ep_ideal['m_dot']/physical_constants.g0
        # Print ideal values
        print("\n --- IDEAL PARAMETERS ---")
        print(" Thrust: {:1.1f} mN".format(ep_ideal['thrust']*1e3))
        print(" Isp: {:3.1f} s".format(Isp))
        print(" Mass flow: {:3.3f} mg/s".format(m_dot*1e6))
        print(" Throat width: {:4.1f} micron".format(w_throat*1e6))

        
        
        # Remaining paramters
        inlet_manifold_length_factor = 2 # [m] Multiplication factor with inlet manifold width to determine manifold length
        inlet_manifold_width_factor = 5.5 # [-] Multiplication factor (with channel width to determine margin in chamber)
        l_exit_manifold = 0 # [m] Length between the end of multiple channels and start of convergent nozzle
        convergent_half_angle  = math.radians(45) # [rad] Half-angle of the convergent part of the nozzle
        divergent_half_angle = math.radians(22.5) # [rad] Half-angle of divergent part of the nozzle
        w_channel_spacing = 1e-6 # [m] Spacing between channels (wall-to-wall)
        w_outer_margin = 2e-3 # [m] Margin around the outer channels for structural integrity
        channel_amount = 24 # [-] Number of channels
        emissivity_chip_top = 0.5 # [-] Assumed emissivity of chip at top-side

    # Desired geometric values
    
    
    
    
    # Used heat transfer relations
    Nusselt_relations = {
        'Nu_func_gas': thermo.convection.Nu_laminar_developed_constant_wall_temp_square, # [-] Function to calculate Nusselt number (gas phase)
        'Nu_func_liquid': thermo.convection.Nu_laminar_developed_constant_wall_temp_square,  # [-] Function to caculate Nusselt number (liquid phase)
        'Nu_func_two_phase': tp.Nu_Kandlikar_NBD_dryout, # [-] Function to calculate Nusselt number (two-phase)
        'Nu_func_le': thermo.convection.Nu_laminar_developed_constant_wall_temp_square, # [-] Function to calculate Nusselt in two-phase, AS IF the flow was entirely liquid (two-phase, le)
        'Nu_func_dryout': thermo.convection.Nu_laminar_developed_constant_wall_temp_square, #thermo.two_phase.Nu_DB_two_phase # [-] Function to calculate the Nusselt number after dry-out. It is up to Nu_func_two_phase to decide if/how to apply it'
    }

    pressure_drop_relations = {
        'l': oneD.calc_single_phase_frictional_pressure_drop_low_Reynolds,
        'tp': oneD.calc_two_phase_frictional_pressure_drop_low_Reynolds,
        'g': oneD.calc_single_phase_frictional_pressure_drop_low_Reynolds,
        'contraction': oneD.calc_single_phase_contraction_pressure_drop_Kawahara2015
    }

    # Fidelity of simulation
    steps_per_section = 50 # [-] Amount of subdivision in each section 
    steps_l = steps_per_section
    steps_tp = steps_per_section
    steps_g = steps_per_section

    prepared_values = oneD.full_homogenous_preparation(
        T_inlet=T_inlet,
        T_outlet=T_chamber, # Outlet of the channel is equal to chamber temperature in IRT
        m_dot=m_dot,
        p_ref=p_inlet, # The pressure in the channel is assumed to be constant
        steps_l=steps_l,
        steps_tp=steps_tp,
        steps_g=steps_g,
        fp=fp
    )


    # Power loss function to optimzie
    func_P_loss = lambda T_superheat : optim_P_total( 
        channel_amount=channel_amount,
        w_channel=w_channel,
        h_channel=h_channel,
        inlet_manifold_length_factor=inlet_manifold_length_factor,
        inlet_manifold_width_factor=inlet_manifold_width_factor,
        l_exit_manifold=l_exit_manifold,
        w_channel_spacing=w_channel_spacing,
        w_outer_margin=w_outer_margin,
        T_wall = T_chamber + T_superheat, # [K] Wall temperature must of course be higher than chamber temperature, so a variable T_superheater that must be larger than 0 is useful
        p_ref=p_inlet,
        m_dot=m_dot,
        prepared_values=prepared_values,
        Nusselt_relations=Nusselt_relations,
        pressure_drop_relations=pressure_drop_relations,
        convergent_half_angle=convergent_half_angle,
        divergent_half_angle=divergent_half_angle,
        F_desired=F_desired,
        p_back=p_back,
        AR_exit=AR_exit,
        w_throat=w_throat,
        emissivity_top=emissivity_chip_top,
        fp=fp
    )

    # Fist make a plot to see what's going on
    T_range = np.linspace(start=10, stop=500, num=50) # [K]
    
    # Store calculated heat loss here
    P_loss = np.zeros_like(T_range) * np.nan # [W] Heat loss
    l_channel = np.zeros_like(T_range) * np.nan # [m] Heating channel length (not the chip length)
    l_total = np.zeros_like(T_range) * np.nan # [m] Total chip length
    w_total = np.zeros_like(T_range) * np.nan # [m] Total chip width.
    pressure_drop = np.zeros_like(T_range) *np.nan # [Pa] Pressure drop over channel (only if it drops stagnation pressure)
    is_channel_choked = np.zeros_like(T_range) * np.nan # [-] Going to be 0 if not choked, and 1 if chocked
    Re_channel_exit = np.zeros_like(T_range) * np.nan # [-] Reynolds number at throat
    dP_contraction = np.zeros_like(T_range) * np.nan # [Pa] Pressure drop over sudden contraction at entrance
    #F = np.zeros_like(T_range) # [N] To store thrust that is iterated to, for debugging purposes
    #CF = np.zeros_like(T_range) # [-] To store thrust coefficient

    iter_T = np.nditer(T_range, flags=['c_index'])

    for T in iter_T:
        res = func_P_loss(T)
        

        P_loss[iter_T.index] = res['P_loss']
        l_channel[iter_T.index] = res['l_channel']
        pressure_drop[iter_T.index] = res['pressure_drop']
        dP_contraction[iter_T.index] = res['dP_contraction']
        l_total[iter_T.index] = res['l_total']
        w_total[iter_T.index] = res['w_total']


        # Check if channels are choked? 
        is_channel_choked[iter_T.index] = res['is_channel_choked'] # []
        # See if Reynolds number still makes sense (laminar or turbulent?)
        Re_channel_exit[iter_T.index] = res['Re_channel_exit']



    print("Minumum heat loss: {:2.3f} W".format(np.nanmin(P_loss)))

    #ax2 = ax1.twinx()
    ax_P_loss.plot(T_range+T_chamber, P_loss, label="{:3.0f}".format(w_channel*1e6))
    
    fig_P_loss.suptitle("Heat loss ($\\dot{{m}}={:3.2f}$ mg$\\cdot$s$^{{-1}}$, $p={:1.2f}$ bar, $T={:3.0f}$ K)".format(m_dot*1e6, p_inlet*1e-5, T_chamber))
    ax_P_loss.set_xlabel("Wall temp [K]")
    ax_P_loss.set_ylabel("Heat loss [W]")

    ax_l_total.plot(T_range+T_chamber, l_channel*1e3, label="{:3.0f}".format(w_channel*1e6))
    ax_l_total.set_ylabel("Channel length [mm]")
    ax_l_total.set_xlabel("Wall temp [K]")
    fig_l_total.suptitle("Channel length ($\\dot{{m}}={:3.2f}$ mg$\\cdot$s$^{{-1}}$, $p={:1.2f}$ bar, $T={:3.0f}$ K)".format(m_dot*1e6, p_inlet*1e-5, T_chamber))

    ax_pressure_drop.plot(T_range+T_chamber, dP_contraction/pressure_drop, label="{:3.0f}".format(w_channel*1e6))
    ax_pressure_drop.set_xlabel("Wall temp [K]")
    ax_pressure_drop.set_ylabel("Pressure drop contraction (fraction) [-]")
    fig_pressure_drop.suptitle("Pressure drop contraction percentage($\\dot{{m}}={:3.2f}$ mg$\\cdot$s$^{{-1}}$, $p={:1.2f}$ bar, $T={:3.0f}$ K)".format(m_dot*1e6, p_inlet*1e-5, T_chamber))



    ax_is_choked.plot(T_range+T_chamber, is_channel_choked, marker='o', linestyle=" ", label="{:3.0f}".format(w_channel*1e6))
    ax_is_choked.set_xlabel("Wall temp [K]")
    ax_is_choked.set_ylabel("Choked ?")
    fig_is_choked.suptitle("Choked? ($\\dot{{m}}={:3.2f}$ mg$\\cdot$s$^{{-1}}$, $p={:1.2f}$ bar, $T={:3.0f}$ K)".format(m_dot*1e6, p_inlet*1e-5, T_chamber))
 


    ax_w_total.plot(T_range+T_chamber, w_total*1e3, label="{:3.0f}".format(w_channel*1e6))
    ax_w_total.set_xlabel("Wall temp [K]")
    ax_w_total.set_ylabel("Chip width [mm] ")
    fig_w_total.suptitle("Chip width ($\\dot{{m}}={:3.2f}$ mg$\\cdot$s$^{{-1}}$, $p={:1.2f}$ bar, $T={:3.0f}$ K)".format(m_dot*1e6, p_inlet*1e-5, T_chamber))

    ax_Re_channel_exit.plot(T_range+T_chamber, Re_channel_exit, label="{:3.0f}".format(w_channel*1e6))
    ax_Re_channel_exit.set_xlabel("Wall temp [K]")
    ax_Re_channel_exit.set_ylabel("Reynolds Channel [mm] ")
    fig_Re_channel_exit.suptitle("Reynolds channel($\\dot{{m}}={:3.2f}$ mg$\\cdot$s$^{{-1}}$, $p={:1.2f}$ bar, $T={:3.0f}$ K)".format(m_dot*1e6, p_inlet*1e-5, T_chamber))
  

    # plt.show()


def optim_P_total(channel_amount, w_channel, h_channel, inlet_manifold_length_factor, inlet_manifold_width_factor, l_exit_manifold, w_channel_spacing, w_outer_margin, T_wall, p_ref, m_dot, \
    prepared_values, Nusselt_relations, pressure_drop_relations, convergent_half_angle, divergent_half_angle, F_desired, p_back, AR_exit, w_throat, emissivity_top, fp: FluidProperties):
    # First calculate area ratio at entrance so it can be provided to contraction pressure drop function
    w_inlet_manifold = chamber.inlet_manifold_width(\
            channel_amount=channel_amount,
            w_channel=w_channel,
            w_channel_spacing=w_channel_spacing,
            inlet_manifold_width_factor=inlet_manifold_width_factor) # [m] Total width of inlet manifold
    area_ratio_contraction = chamber.area_ratio_contraction(\
        w_inlet_manifold=w_inlet_manifold,
        w_channel=w_channel,
        channel_amount=channel_amount) # [-] Area ratio of sudden contraction
    
    # Calculate heat transfer to determine channel length
    res = oneD.rectangular_multi_channel_homogenous_calculation(\
        channel_amount=channel_amount,
        prepared_values=prepared_values,
        Nusselt_relations=Nusselt_relations,
        pressure_drop_relations=pressure_drop_relations,
        w_channel=w_channel,
        h_channel=h_channel,
        m_dot=m_dot,
        T_wall=T_wall,
        p_inlet=p_ref,
        fp=fp,
        area_ratio_contraction=area_ratio_contraction)
    
    l_channel = res['L_total'] # [m] Channel length
    p_chamber = res['p_chamber'] # [Pa] Pressure out channel outlet/inlet of nozzle according to IRT
    T_chamber = res['T_chamber']

    # Nozzle performance must be recalculated, if pressure drop occured
    assert(pressure_drop_relations is not None) # Just check is one was calculated

    if p_chamber <= 0: # If pressure drop was so large p_chamber was negative, the solution is invalid, and NaNs must be returned
        # Return the same dictionary, but all NaN.
        return {
            'P_loss': np.nan, # [W] Current approximation of heat loss
            'l_channel': np.nan, # [m] Channel length
            'pressure_drop': np.nan, # [Pa] Total pressure drop
            'dP_contraction': np.nan, # [Pa] Pressure drop de to contraction
            'w_total': np.nan, # [m] Total chip width
            'l_total': np.nan, # [m] Total chip length
            'is_channel_choked': np.nan,
            'Re_channel_exit': np.nan, # [-] Reynolds number at channel exit
    }
    # This only runs if pressure drop was not too large.
    # In that case the throat area must be recalculated still
    else:
            

        ep = engine_performance_from_F_and_T(
                F_desired=F_desired,
                p_chamber=p_chamber,
                T_chamber=T_chamber,
                AR_exit=AR_exit,
                p_back=p_back,
                fp=fp
            )

        # New throat widtl_inlet_manifoldmined from desired performance parameter
        w_throat_new = ep['A_throat']/h_channel # [m]
        # Check if new pressure, does not cause choked heating channels
        is_channel_choked = (w_channel*h_channel*channel_amount < ep['A_throat']) # [bool] If combined channel area is smaller than throat, they are choked.
        # Check new Reynolds number at channel exit for laminar/turbulent flow assumptions
        hydraulic_diameter = chamber.hydraulic_diameter_rectangular(w_channel=w_channel,h_channel=h_channel)
        Re_channel_exit = fp.get_Reynolds_from_mass_flow(T=T_chamber, p=p_chamber,L_ref=hydraulic_diameter, m_dot=m_dot, A=h_channel*w_channel)

        l_outlet = chamber.outlet_length(\
            w_channel=w_channel,
            w_channel_spacing=w_channel_spacing,
            channel_amount=channel_amount,
            convergent_half_angle=convergent_half_angle,
            w_throat=w_throat_new,
            divergent_half_angle=divergent_half_angle,
            AR_exit=AR_exit,
            l_exit_manifold=l_exit_manifold
            )
        
        
        l_inlet_manifold = chamber.inlet_manifold_length(w_inlet_manifold=w_inlet_manifold, inlet_manifold_length_factor=inlet_manifold_length_factor)
        l_total = chamber.total_chip_length(l_inlet_manifold=l_inlet_manifold, l_channel=l_channel, l_outlet=l_outlet)
        w_total = chamber.total_chip_width(\
            w_inlet_manifold=w_inlet_manifold,
            w_outer_margin=w_outer_margin,
            w_throat=w_throat_new,
            AR_exit=AR_exit)

        A_chip = l_total*w_total # [m^2]

        P_rad_loss_top = 2*chamber.radiation_loss(T=T_wall, A=A_chip, emmisivity=emissivity_top) # [W]  Radiation loss through top of chip

        return {'P_loss': 2*P_rad_loss_top, # [W] Current approximation of heat loss
                'l_channel': l_channel, # [m] Channel length
                'pressure_drop': res['dP_total'], # [Pa] Total pressure drop
                'dP_contraction': res['dP_contraction'], # [Pa] Pressure drop due to contraction
                'w_total': w_total, # [m] Total chip width
                'l_total': l_total, # [m] Total chip length
                'is_channel_choked': is_channel_choked, # [bool] 0 and 1 in this case, to check if channels are not too small
                'Re_channel_exit': Re_channel_exit, # [-] Check Reynolds number to see if laminar flow assumptions still hold 
        }

# Run script
if __name__ == "__main__":

    # Go through several channel widths to plot

    w_channel = ( 5e-6, 10e-6, 25e-6, 50e-6, 75e-6, 100e-6, 200e-6, 300e-6)

    # Create and figures and axes to plot on
    fig_P_loss = plt.figure()
    ax_P_loss = fig_P_loss.gca()
    ax_P_loss.set_ylim(0, 12)

    fig_w_total = plt.figure()
    ax_w_total = fig_w_total.gca()
    ax_w_total.set_ylim(5,10)

    fig_pressure_drop = plt.figure()
    ax_pressure_drop = fig_pressure_drop.gca()
    ax_pressure_drop.set_ylim(0,5)

    fig_is_choked = plt.figure()
    ax_is_choked = fig_is_choked.gca()
    ax_is_choked.set_ylim(-0.1,1.1)

    fig_l_total = plt.figure()
    ax_l_total = fig_l_total.gca()
    ax_l_total.set_ylim(0,50)

    # Reynolds channel exit would be interesting
    fig_Re_channel_exit = plt.figure()
    ax_Re_channel_exit = fig_Re_channel_exit.gca()
    #ax_Re_channel_exit.set_ylim(10,50)
    
    for w in w_channel:
        calc_and_plot_channel_width(w_channel=w, ax_P_loss=ax_P_loss, fig_P_loss=fig_P_loss,\
            ax_w_total=ax_w_total, fig_w_total=fig_w_total, ax_l_total=ax_l_total,fig_l_total=fig_l_total,\
                ax_pressure_drop=ax_pressure_drop, fig_pressure_drop=fig_pressure_drop,\
                    ax_is_choked=ax_is_choked, fig_is_choked=fig_is_choked,
                    ax_Re_channel_exit=ax_Re_channel_exit, fig_Re_channel_exit=fig_Re_channel_exit)
    
    ax_P_loss.legend(title="$w_{{channel}} [\\mu$m$]$")
    ax_P_loss.grid()
    fig_P_loss.tight_layout()

    ax_l_total.legend(title="$w_{{channel}} [\\mu$m$]$")
    ax_l_total.grid()
    fig_l_total.tight_layout()

    ax_w_total.legend(title="$w_{{channel}} [\\mu$m$]$")
    ax_w_total.grid()
    fig_w_total.tight_layout()

    ax_pressure_drop.legend(title="$w_{{channel}} [\\mu$m$]$")
    ax_pressure_drop.grid()
    fig_pressure_drop.tight_layout()

    ax_is_choked.legend(title="$w_{{channel}} [\\mu$m$]$")
    ax_is_choked.grid()
    fig_is_choked.tight_layout()

    ax_Re_channel_exit.legend(title="$w_{{channel}} [\\mu$m$]$")
    ax_Re_channel_exit.grid()
    fig_Re_channel_exit.tight_layout()

    plt.show()