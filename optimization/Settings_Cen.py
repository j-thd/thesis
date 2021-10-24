"""This file contains the settings for various optimization tasks
    """

import math
from thermo.prop import FluidProperties
import thermo.convection
import thermo.two_phase as tp
from thermo.prop import FluidProperties
import models.one_D as oneD

#####################################################################################
###                        1D RECTANGULAR MULTICHANNEL 4mN                         ##
#####################################################################################
settings_1D_rectangular_multichannel = {}
## Fixed design parameters
settings_1D_rectangular_multichannel['FDP'] = {
    'fp': FluidProperties('water'),                 # (object) Interface with Coolprop that returns all necessary properties of water
    'p_inlet': 5e5,                                 # [Pa] Pressure at inlet
    'p_back': 0,                                    # [Pa] NOTE: back pressure MUST be vacuum/zero for nozzle adjustment to work correctly
    'T_inlet': 300,                                 # [K] Heating chamber inlet temperature
    'h_channel': 100e-6,                            # [m] Channel depth/height
    'AR_exit': 10,                                  # [-] Exit area ratio of nozzles
    'inlet_manifold_length_factor' : 2,             # [m] Multiplication factor with inlet manifold width to determine manifold length
    'inlet_manifold_width_factor' : 5.5,            # [-] Multiplication factor (with channel width to determine margin in chamber)
    'l_exit_manifold' : 0,                          # [m] Length between the end of multiple channels and start of convergent nozzle
    'convergent_half_angle'  : math.radians(45),    # [rad] Half-angle of the convergent part of the nozzle
    'divergent_half_angle' : math.radians(22.5),    # [rad] Half-angle of divergent part of the nozzle
    'w_outer_margin' : 0.5e-3,                        # [m] Margin around the outer channels for structural integrity
    'emissivity_chip_top' : 0.7,                    # [-] Emissivity of chip at top-side
    'emissivity_chip_bottom' : 0.7,                 # [-] Emissivity of chip at bottom-side
    'kappa_wall': 50,                              # [W/(m*K)] Thermal conductivty of chip/wall in general

}
## BOUNDARIES ON DESIGN SPACE
settings_1D_rectangular_multichannel['bounds'] = {
    'w_channel_spacing': (50e-6,200e-6),
    'channel_amount': (1,25),
    'w_channel': (20e-6, 1e-3),
    'T_wall_superheat': (1, 500),           # [K] Wall temperature range (T_wall = T_chamber + T_wall_superheat)
}
## HEAT TRANSFER RELATIONS IN DIFFERENT SECTIONS OF THE CHANNEL
settings_1D_rectangular_multichannel['Nusselt_relations'] = {
        'Nu_func_gas': thermo.convection.Nu_laminar_developed_constant_wall_temp_square, # [-] Function to calculate Nusselt number (gas phase)
        'Nu_func_liquid': thermo.convection.Nu_laminar_developed_constant_wall_temp_square,  # [-] Function to caculate Nusselt number (liquid phase)
        'Nu_func_two_phase': tp.Nu_Kandlikar_NBD_dryout, # [-] Function to calculate Nusselt number (two-phase)
        'Nu_func_le': thermo.convection.Nu_laminar_developed_constant_wall_temp_square, # [-] Function to calculate Nusselt in two-phase, AS IF the flow was entirely liquid (two-phase, le)
        'Nu_func_dryout': thermo.convection.Nu_laminar_developed_constant_wall_temp_square, #thermo.two_phase.Nu_DB_two_phase # [-] Function to calculate the Nusselt number after dry-out. It is up to Nu_func_two_phase to decide if/how to apply it'
    }
## PRESSURE DROP RELATIONS IN DIFFERENT SECTIONS OF THE CHANNEL
settings_1D_rectangular_multichannel['pressure_drop_relations'] = {
        'l': oneD.calc_single_phase_frictional_pressure_drop_low_Reynolds,
        'tp': oneD.calc_two_phase_frictional_pressure_drop_low_Reynolds,
        'g': oneD.calc_single_phase_frictional_pressure_drop_low_Reynolds,
        'contraction': oneD.calc_single_phase_contraction_pressure_drop_Kawahara2015
    }
## SIMULATION SETTINGS
steps_per_section = 1000 # [-] Amount of subdivisions in the channel
settings_1D_rectangular_multichannel['steps'] = {
    'steps_l' : steps_per_section, # [-] Liquid section
    'steps_tp' : steps_per_section, # [-] Two-phase section
    'steps_g' : steps_per_section, # [-] Gas section
}

settings_1D_rectangular_multichannel['function_scaling'] = 1e10 # [-] Amount of scaling of optimization function 

#####################################################################################
###                        1D RECTANGULAR MULTICHANNEL 1mN                         ##
#####################################################################################
settings_1D_rectangular_multichannel_1mN = {}
## Fixed design parameters
settings_1D_rectangular_multichannel_1mN['FDP'] = {
    'fp': FluidProperties('water'),                 # (object) Interface with Coolprop that returns all necessary properties of water
    'p_inlet': 5e5,                                 # [Pa] Pressure at inlet
    'p_back': 0,                                    # [Pa] NOTE: back pressure MUST be vacuum/zero for nozzle adjustment to work correctly
    'T_inlet': 300,                                 # [K] Heating chamber inlet temperature
    'h_channel': 100e-6,                            # [m] Channel depth/height
    'AR_exit': 10,                                  # [-] Exit area ratio of nozzles
    'inlet_manifold_length_factor' : 2,             # [m] Multiplication factor with inlet manifold width to determine manifold length
    'inlet_manifold_width_factor' : 5.5,            # [-] Multiplication factor (with channel width to determine margin in chamber)
    'l_exit_manifold' : 0,                          # [m] Length between the end of multiple channels and start of convergent nozzle
    'convergent_half_angle'  : math.radians(45),    # [rad] Half-angle of the convergent part of the nozzle
    'divergent_half_angle' : math.radians(22.5),    # [rad] Half-angle of divergent part of the nozzle
    'w_outer_margin' : 0.5e-3,                        # [m] Margin around the outer channels for structural integrity
    'emissivity_chip_top' : 0.7,                    # [-] Emissivity of chip at top-side
    'emissivity_chip_bottom' : 0.7,                 # [-] Emissivity of chip at bottom-side
    'kappa_wall': 50,                              # [W/(m*K)] Thermal conductivty of chip/wall in general

}
## BOUNDARIES ON DESIGN SPACE
settings_1D_rectangular_multichannel_1mN['bounds'] = {
    'w_channel_spacing': (50e-6,500e-6),
    'channel_amount': (1,25),
    'w_channel': (20e-6, 1e-3),
    'T_wall_superheat': (1, 500),           # [K] Wall temperature range (T_wall = T_chamber + T_wall_superheat)
}
## HEAT TRANSFER RELATIONS IN DIFFERENT SECTIONS OF THE CHANNEL
settings_1D_rectangular_multichannel_1mN['Nusselt_relations'] = {
        'Nu_func_gas': thermo.convection.Nu_laminar_developed_constant_wall_temp_square, # [-] Function to calculate Nusselt number (gas phase)
        'Nu_func_liquid': thermo.convection.Nu_laminar_developed_constant_wall_temp_square,  # [-] Function to caculate Nusselt number (liquid phase)
        'Nu_func_two_phase': tp.Nu_Kandlikar_NBD_dryout, # [-] Function to calculate Nusselt number (two-phase)
        'Nu_func_le': thermo.convection.Nu_laminar_developed_constant_wall_temp_square, # [-] Function to calculate Nusselt in two-phase, AS IF the flow was entirely liquid (two-phase, le)
        'Nu_func_dryout': thermo.convection.Nu_laminar_developed_constant_wall_temp_square, #thermo.two_phase.Nu_DB_two_phase # [-] Function to calculate the Nusselt number after dry-out. It is up to Nu_func_two_phase to decide if/how to apply it'
    }
## PRESSURE DROP RELATIONS IN DIFFERENT SECTIONS OF THE CHANNEL
settings_1D_rectangular_multichannel_1mN['pressure_drop_relations'] = {
        'l': oneD.calc_single_phase_frictional_pressure_drop_low_Reynolds,
        'tp': oneD.calc_two_phase_frictional_pressure_drop_low_Reynolds,
        'g': oneD.calc_single_phase_frictional_pressure_drop_low_Reynolds,
        'contraction': oneD.calc_single_phase_contraction_pressure_drop_Kawahara2015
    }
## SIMULATION SETTINGS
steps_per_section = 1000 # [-] Amount of subdivisions in the channel
settings_1D_rectangular_multichannel_1mN['steps'] = {
    'steps_l' : steps_per_section, # [-] Liquid section
    'steps_tp' : steps_per_section, # [-] Two-phase section
    'steps_g' : steps_per_section, # [-] Gas section
}

settings_1D_rectangular_multichannel_1mN['function_scaling'] = 1e10 # [-] Amount of scaling of optimization function 