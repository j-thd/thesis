"""Thrusters data is very often re-used in scripts, so it would be nice to have some dictionaries with their values preset
"""

import math

## Thruster data from Huib's thruster
td_Huib_TTH_4_1 = {
    "name": "useful name here",
    "propellant": "nitrogen",
    "F": 11.839e-3, # [N]
    "p_chamber": 1.66e5, # [Pa]
    "p_back": 209, # [Pa]
    "T_chamber": 670.7, # [K]
    "m_dot": 2.083462e-8, # [kg/s]
    "h_throat": 500e-6, # [m]
    "w_throat": 130e-6, # [m]
    "throat_roc": 260e-6, # [m]
    "AR_exit": 8.25, # [-]
    "w_nozzle_exit": 1072.5, # [m]
    "divergent_half_angle": math.radians(20) # [rad]
    }

## Thruster number 5 from Silva2017
td_Silva_5 = {
    "name": "Silva2017 - #5 Large serpentine 5-channel",
    "propellant": "water",
    "F": 0.67e-3, # [N]
    "m_dot": 0.55e-6, # [kg/s]
    "T_chamber": 423.03 , # [K]
    "T_wall": 273.15+200, # [K] Wall temperature is hard to determine from paper (seems low)
    "channel_amount": 5, # [-]
    "p_inlet": 4.8e5, # [Pa]
    "T_inlet": 297.15, # [K] Room temperature of 24 degrees mentioned in report
    "p_back": 0, # [Pa]
    "h_channel": 100e-6, # [m]
    "w_channel": 0.212e-3, # [m]
    "L_channel": 14*2*math.pi*0.160e-3, # [m]  # Serpentine channel. The length is path length of center-line of the channel.
    "h_throat": 100e-6, # [m]
    "w_throat": 130e-6, # [m]
    "throat_roc": 260e-6, # [m]
    "AR_exit": 8.25, # [-]
    "w_nozzle_exit": 1072.5, # [m]
    "divergent_half_angle": math.radians(20) # [rad]
    }

## CEN2010 THRUSTERS ## STRAIGHT MULTI-CHANNEL THRUSTERS
Cen2010_1 = {
    "name": "Cen2010 2.33 mg/s",
    "propellant": "water",
    "F": 2.3e-3, # [N]
    "m_dot": 2.33e-6, # [kg/s]
    "T_chamber": None , # [K]
    "T_wall": 573.15, # [K] TBD
    "channel_amount": 9, # [-]
    "p_inlet": 1.28e5, # [Pa]
    "T_inlet": 300, # [K] TBD
    "p_back": 0, # [Pa]
    "h_channel": 120e-6, # [m]
    "w_channel": 80e-6, # [m]
    "L_channel": 6e-3, # [m]  
    "h_throat": 120e-6, # [m]
    "w_throat": 150e-6, # [m]
    "AR_exit": 11.72, # [-]
    "w_nozzle_exit": 1758e-6, # [m]
    "L_nozzle": 3000e-6, # [m] Nozzle length
    "divergent_half_angle": math.radians(15), # [rad]
    "T_chamber_guess": 571.90, # [K] Chamber temperature is not measured, this is based on laminar 1D model for validation
    }

Cen2010_2 = {
    "name": "Cen2010 3.00mg/s",
    "propellant": "water",
    "F": 3.0e-3, # [N]
    "m_dot": 3.00e-6, # [kg/s]
    "T_chamber": None , # [K]
    "T_wall": 573.15, # [K] TBD
    "channel_amount": 9, # [-]
    "p_inlet": 1.52e5, # [Pa]
    "T_inlet": 300, # [K] TBD
    "p_back": 0, # [Pa]
    "h_channel": 120e-6, # [m]
    "w_channel": 80e-6, # [m]
    "L_channel": 6e-3, # [m]  
    "h_throat": 120e-6, # [m]
    "w_throat": 150e-6, # [m]
    "AR_exit": 11.72, # [-]
    "w_nozzle_exit": 1758e-6, # [m]
    "divergent_half_angle": math.radians(15), # [rad]
    "T_chamber_guess": 570.80, # [K] Chamber temperature is not measured, this is based on laminar 1D model for validation
    }

Cen2010_3 = {
    "name": "Cen2010 3.66mg/s",
    "propellant": "water",
    "F": 3.8e-3, # [N]
    "m_dot": 3.66e-6, # [kg/s]
    "T_chamber": None , # [K]
    "T_wall": 573.15, # [K] TBD
    "channel_amount": 9, # [-]
    "p_inlet": 1.80e5, # [Pa]
    "T_inlet": 300, # [K] TBD
    "p_back": 0, # [Pa]
    "h_channel": 120e-6, # [m]
    "w_channel": 80e-6, # [m]
    "L_channel": 6e-3, # [m]  
    "h_throat": 120e-6, # [m]
    "w_throat": 150e-6, # [m]
    "AR_exit": 11.72, # [-]
    "w_nozzle_exit": 1758e-6, # [m]
    "divergent_half_angle": math.radians(15), # [rad]
    "T_chamber_guess": 568.86, # [K] Chamber temperature is not measured, this is based on laminar 1D model for validation
    }

Cen2010_4 = {
    "name": "Cen2010 4.33mg/s",
    "propellant": "water",
    "F": 4.5e-3, # [N]
    "m_dot": 4.33e-6, # [kg/s]
    "T_chamber": None , # [K]
    "T_wall": 573.15, # [K] TBD
    "channel_amount": 9, # [-]
    "p_inlet": 2.08e5, # [Pa]
    "T_inlet": 300, # [K] TBD
    "p_back": 0, # [Pa]
    "h_channel": 120e-6, # [m]
    "w_channel": 80e-6, # [m]
    "L_channel": 6e-3, # [m]  
    "h_throat": 120e-6, # [m]
    "w_throat": 150e-6, # [m]
    "AR_exit": 11.72, # [-]
    "w_nozzle_exit": 1758e-6, # [m]
    "divergent_half_angle": math.radians(15), # [rad]
    "T_chamber_guess": 565.95, # [K] Chamber temperature is not measured, this is based on laminar 1D model for validation
    }

Cen2010_5 = {
    "name": "Cen2010 5.00mg/s",
    "propellant": "water",
    "F": 5.2e-3, # [N]
    "m_dot": 5.00e-6, # [kg/s]
    "T_chamber": None , # [K]
    "T_wall": 573.15, # [K] TBD
    "channel_amount": 9, # [-]
    "p_inlet": 2.28e5, # [Pa]
    "T_inlet": 300, # [K] TBD
    "p_back": 0, # [Pa]
    "h_channel": 120e-6, # [m]
    "w_channel": 80e-6, # [m]
    "L_channel": 6e-3, # [m]  
    "h_throat": 120e-6, # [m]
    "w_throat": 150e-6, # [m]
    "AR_exit": 11.72, # [-]
    "w_nozzle_exit": 1758e-6, # [m]
    "divergent_half_angle": math.radians(15), # [rad]
    "T_chamber_guess": 562.22, # [K] Chamber temperature is not measured, this is based on laminar 1D model for validation
    }

Cen2010_6 = {
    "name": "Cen2010 5.66mg/s",
    "propellant": "water",
    "F": 5.7e-3, # [N]
    "m_dot": 5.66e-6, # [kg/s]
    "T_chamber": None , # [K]
    "T_wall": 573.15, # [K] TBD
    "channel_amount": 9, # [-]
    "p_inlet": 2.52e5, # [Pa]
    "T_inlet": 300, # [K] TBD
    "p_back": 0, # [Pa]
    "h_channel": 120e-6, # [m]
    "w_channel": 80e-6, # [m]
    "L_channel": 6e-3, # [m]  
    "h_throat": 120e-6, # [m]
    "w_throat": 150e-6, # [m]
    "AR_exit": 11.72, # [-]
    "w_nozzle_exit": 1758e-6, # [m]
    "divergent_half_angle": math.radians(15), # [rad]
    "T_chamber_guess": 557.80, # [K] Chamber temperature is not measured, this is based on laminar 1D model for validation
    }

td_verification_one = {
    "name":             "Fictional thruster for verification of two_phase_single_channel()",
    "propellant":   "water",
    "m_dot":            1e-6, # [kg/s]
    "p_inlet":          5e5, # [Pa]
    "T_inlet":          300, # [K]
    "T_chamber":        500,  # [K]
    "T_wall":           600, # [K]
    "h_channel":        100e-6, # [m]
    "w_channel":        100e-6, # [m]
}

td_high_kinetic_energy = {
    "name":             "Fictional thruster for worst-case pressure drop calculations",
    "propellant":       "water",
    "m_dot":            3.6e-6, # [kg/s]
    "channel_amount":   5, # [-]
    "p_inlet":          5e5, # [Pa]
    "T_inlet":          300, # [K]
    "T_chamber":        450, # [K]
    "h_channel":        100e-6, # [m]
    "w_channel":        100e-6, # [m]
    "L_channel":        20e-3, # [m]
}