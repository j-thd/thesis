import unittest

import models.zero_D as zD
import thrusters.thruster_data
from thermo.prop import FluidProperties
import thermo.convection

class TestTwoPhaseSingleChannel(unittest.TestCase):
    def testDummyWithWebBookData(self):
        td = thrusters.thruster_data.td_verification_one # Dummy dictionary with thruster data for verifcation

        Nu_func_gas = thermo.convection.Nu_DB
        Nu_func_liquid = thermo.convection.Nu_DB

        T_wall = td['T_wall']
        w_channel = td['w_channel']
        T_inlet = td['T_inlet']
        T_chamber = td['T_chamber']
        p_ref = td['p_inlet']
        m_dot = td['m_dot']
        h_channel = td['h_channel']
        fp = FluidProperties(td['propellant'])

        # Expected values for liquid and multi/phase portion of channel
        exp_Q_dot_liquid_multi = 2.635 # [W]
        exp_T_bulk_liquid_multi = 362.49 # [K]
        exp_u_bulk_liquid_multi = 0.1035271 # [m/s]
        exp_rho_bulk_liquid_multi = 965.93 # [kg/m^3]
        exp_Re_liquid_multi = 31.5556 # [-]
        exp_Pr_liquid_multi = 1.973033 # [-]
        exp_Nu_liquid_multi = 0.4775 # [-]
        exp_St_liquid_multi = 0.00767 # [-]
        exp_h_conv_liquid_multi = 3224.5 # [W/(m^2*K)]
        exp_A_heater_liquid_multi = 3.4407e-6 # [m^2]
        exp_L_channel_liquid_multi = 0.0086018 # [m]

        # Expected values for gas portion of channel
        exp_Q_dot_gas = 0.1646 # [W]
        exp_T_bulk_gas = 462.49 # [K]
        exp_u_bulk_gas = 41.423 # [m/2]
        exp_rho_bulk_gas = 2.4141 # [kg/m^3]
        exp_Re_bulk_gas = 640.70 # [-]
        exp_Pr_bulk_gas = 0.99236 # [-]
        exp_Nu_gas = 4.0338 # [-]
        exp_St_gas = 0.0063444 # [-]
        exp_h_conv_gas = 1377.4 # [W/(m^2*K)]
        exp_A_heater_gas = 8.6903e-7 # [m^2]
        exp_L_channel_gas = 2.1726e-3 # [m]

        # Expected value for total
        exp_L_channel = 10.774e-3 # [m]


        res = zD.two_phase_single_channel(
            T_wall=T_wall,
            w_channel=w_channel,
            Nu_func_gas=Nu_func_gas,
            Nu_func_liquid=Nu_func_liquid,
            T_inlet=T_inlet,
            T_chamber=T_chamber,
            p_ref=p_ref,
            m_dot=m_dot,
            h_channel=h_channel,
            fp=fp,
            print_info=False)
        
        
        # Checking liquid/multi-phase results
        self.assertAlmostEqual(exp_Q_dot_liquid_multi, res['Q_dot_liquid_multi'], delta=0.001*exp_Q_dot_liquid_multi)
        self.assertAlmostEqual(exp_T_bulk_liquid_multi, res['T_bulk_liquid_multi'], delta=0.001*exp_T_bulk_liquid_multi)
        self.assertAlmostEqual(exp_rho_bulk_liquid_multi, res['rho_bulk_liquid_multi'], delta=0.001*exp_rho_bulk_liquid_multi)
        self.assertAlmostEqual(exp_u_bulk_liquid_multi, res['u_bulk_liquid_multi'], delta=0.001*exp_u_bulk_liquid_multi)
        self.assertAlmostEqual(exp_Re_liquid_multi, res['Re_bulk_liquid_multi'], delta=0.001*exp_Re_liquid_multi)
        self.assertAlmostEqual(exp_Pr_liquid_multi, res['Pr_bulk_liquid_multi'], delta=0.01*exp_Pr_liquid_multi)
        self.assertAlmostEqual(exp_Nu_liquid_multi, res['Nu_liquid_multi'], delta=0.01*exp_Nu_liquid_multi)
        self.assertAlmostEqual(exp_St_liquid_multi, res['St_liquid_multi'], delta=0.01*exp_St_liquid_multi)
        self.assertAlmostEqual(exp_h_conv_liquid_multi, res['h_conv_liquid_multi'], delta=0.01*exp_h_conv_liquid_multi)
        self.assertAlmostEqual(exp_A_heater_liquid_multi, res['A_heater_liquid_multi'], delta=0.01*exp_A_heater_liquid_multi)
        self.assertAlmostEqual(exp_L_channel_liquid_multi, res['L_channel_liquid_multi'], delta=0.01*exp_L_channel_liquid_multi)

        # Again, but for gas phase
        self.assertAlmostEqual(exp_Q_dot_gas, res['Q_dot_gas'], delta=0.001*exp_Q_dot_gas)
        self.assertAlmostEqual(exp_T_bulk_gas, res['T_bulk_gas'], delta=0.001*exp_T_bulk_gas)
        self.assertAlmostEqual(exp_u_bulk_gas, res['u_bulk_gas'], delta=0.001*exp_u_bulk_gas)
        self.assertAlmostEqual(exp_rho_bulk_gas, res['rho_bulk_gas'], delta=0.001*exp_rho_bulk_gas)
        self.assertAlmostEqual(exp_Re_bulk_gas, res['Re_bulk_gas'], delta=0.001*exp_Re_bulk_gas)
        self.assertAlmostEqual(exp_Pr_bulk_gas, res['Pr_bulk_gas'], delta=0.013*exp_Pr_bulk_gas)
        self.assertAlmostEqual(exp_Nu_gas, res['Nu_gas'], delta=0.01*exp_Nu_gas)
        self.assertAlmostEqual(exp_St_gas, res['St_gas'], delta=0.01*exp_St_gas)
        self.assertAlmostEqual(exp_h_conv_gas, res['h_conv_gas'], delta=0.01*exp_h_conv_gas)
        self.assertAlmostEqual(exp_A_heater_gas, res['A_heater_gas'], delta=0.01*exp_A_heater_gas)
        self.assertAlmostEqual(exp_L_channel_gas, res['L_channel_gas'], delta=0.01*exp_L_channel_gas)

        # Checking the total
        self.assertAlmostEqual(exp_L_channel, res['L_channel'], delta = 0.01*exp_L_channel)

    def testIntegrationWithKandlikarRelation(self):
        # Same test but only focussing on Nusselt number for Kandlikar relation, as that would prove thar arguments are correctly passed to function
        td = thrusters.thruster_data.td_verification_one # Dummy dictionary with thruster data for verifcation

        Nu_func_gas = thermo.convection.Nu_DB
        Nu_func_liquid = thermo.convection.Nu_Kandlikar_NDB_Re_low_sat_gas_constant_wall_temp_square_water

        T_wall = td['T_wall']
        w_channel = td['w_channel']
        T_inlet = td['T_inlet']
        T_chamber = td['T_chamber']
        p_ref = td['p_inlet']
        m_dot = td['m_dot']
        h_channel = td['h_channel']
        fp = FluidProperties(td['propellant'])

        # Check just the liquid/multi-phase Nusselt number of the Kandlikar relation, to ensure integration is correct
        exp_Nu_liquid_multi = 46.74856022286813
        

        res = zD.two_phase_single_channel(
            T_wall=T_wall,
            w_channel=w_channel,
            Nu_func_gas=Nu_func_gas,
            Nu_func_liquid=Nu_func_liquid,
            T_inlet=T_inlet,
            T_chamber=T_chamber,
            p_ref=p_ref,
            m_dot=m_dot,
            h_channel=h_channel,
            fp=fp,
            print_info=False)

        self.assertAlmostEqual(exp_Nu_liquid_multi, res['Nu_liquid_multi'], delta=0.01*exp_Nu_liquid_multi)