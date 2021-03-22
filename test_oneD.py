# File to test functions of oneD.py

import unittest
import numpy as np

from thermo.prop import FluidProperties
import thermo.convection
import thermo.two_phase as tp
import models.one_D as oneD


class TestPrepareSinglePhaseLiquid(unittest.TestCase):
    def setUp(self):
        # Inputs
        fp = FluidProperties('water')
        T_inlet = 300 # [K]
        p_inlet = 2e5  # [Pa] Inlet pressure
        steps = 10 # [-] Number of data point
        m_dot = 1e-3 # [kg/s]

        self.prep = oneD.prepare_single_phase_liquid(T_inlet=T_inlet,steps=steps, p_ref = p_inlet, m_dot=m_dot, fp=fp)
        return super().setUp()
    
    def testSaturationValues(self):
        # Reference: https://webbook.nist.gov/cgi/fluid.cgi?Action=Load&ID=C7732185&Type=SatT&Digits=5&PLow=2&PHigh=2&PInc=1&RefState=DEF&TUnit=K&PUnit=bar&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm
        # Saturation values which must be true for the LAST element in the array
        exp_T = 393.36 # [K] Saturation temperature
        exp_rho = 942.94 # [kg/m^3] Density
        exp_h = 504.70e3 # [J/kg] Enthalpy
        exp_mu = 0.00023162 # [Pa*s] Viscosity
        exp_kappa = 0.68321 # [W/(m*K)] Thermal conductivity
        exp_Pr_gas = 1.438755460253802 # [-] Prandtl

        self.assertAlmostEqual(exp_T, self.prep['T'][-1], delta=exp_T*1e-5)
        self.assertAlmostEqual(exp_rho, self.prep['rho'][-1], delta=exp_rho*1e-5)
        self.assertAlmostEqual(exp_h, self.prep['h'][-1], delta=exp_h*1e-5)
        self.assertAlmostEqual(exp_mu, self.prep['mu'][-1], delta=exp_mu*1e-4)
        self.assertAlmostEqual(exp_kappa, self.prep['kappa'][-1], delta=exp_kappa*2.5e-3)
        self.assertAlmostEqual(exp_Pr_gas, self.prep['Pr'][-1], delta=exp_Pr_gas*2.5e-3)

    def testTotalQDot(self):
        # The total power put in must of course equal the specific enthalpy change times mass flow
        # Reference: https://webbook.nist.gov/cgi/fluid.cgi?Action=Load&ID=C7732185&Type=IsoBar&Digits=5&P=2&THigh=500&TLow=300&TInc=100&RefState=DEF&TUnit=K&PUnit=bar&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm
        exp_total_Q_dot = 391.95 # [W] 
        res_total_Q_dot = np.sum(self.prep['Q_dot'])

        self.assertAlmostEqual(exp_total_Q_dot, res_total_Q_dot, delta=exp_total_Q_dot*1e-4)

    def testThirdValue(self):
        # Just to check if values are proper at intermediate sections
        # There are ten values in the array, so 9 steps from T_inlet to T_sat
        #  so the third value should 2/9 along the between saturation and inlet values
        exp_T = (393.36-300)*(2/9) + 300 # [K] 320.746666 

        self.assertAlmostEqual(exp_T, self.prep['T'][2], delta=exp_T*1e-5)
        # This expected value can be used as reference to check other values at this point
        # Ref: https://webbook.nist.gov/cgi/fluid.cgi?Action=Load&ID=C7732185&Type=IsoBar&Digits=5&P=2&THigh=320.746666&TLow=320.746666&TInc=1&RefState=DEF&TUnit=K&PUnit=bar&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm
        # NOTE: Slight difference in reference value as webbook rounds temperature to 320.75 K
        exp_rho = 989.15 # [kg/m^3]
        exp_h = 199.46e3 # [J/kg]
        exp_mu = 0.00056965 # [Pa*s]
        exp_kappa = 0.64072 # [W/(m*K)]
        exp_Pr = 3.7167902125

        self.assertAlmostEqual(exp_rho, self.prep['rho'][2], delta=exp_rho*1e-5)
        self.assertAlmostEqual(exp_h, self.prep['h'][2], delta=exp_h*1e-4)
        self.assertAlmostEqual(exp_mu, self.prep['mu'][2], delta=exp_mu*1e-3)
        self.assertAlmostEqual(exp_kappa, self.prep['kappa'][2], delta=exp_kappa*5-3)
        self.assertAlmostEqual(exp_Pr, self.prep['Pr'][2], delta=exp_Pr*5e-3)
        
    def testDeltaT(self):
        # Check if delta T is evenly spaced as expected
        exp_dT = (393.36-300)/9 # [K] temperature step
        self.assertAlmostEqual(exp_dT, self.prep['dT'], delta=1e-6*exp_dT)

    def testReferenceTemperature(self):
        # Check if if T_ref_i is indeed half-way between T_i-1 and T_i ()
        exp_T_ref = 300 + (393.36-300)*(1/9)/2
        self.assertAlmostEqual(exp_T_ref, self.prep['T_ref'][1],delta=exp_T_ref*1e-5)


class TestCalcSinglePhase(unittest.TestCase):
    # Test single-phase calculations with liquid example as input
    # Results based on manual calculation in excelsheet Verification_one_d.xlsx
    def setUp(self):
        # Inputs
        fp = FluidProperties('water')
        T_inlet = 300 # [K]
        p_inlet = 2e5  # [Pa] Inlet pressure
        steps = 3 # [-] 3 steps for convenience
        m_dot = 0.1e-3 # [kg/s]
        # Get the prepared thermodynamic values
        self.prep = oneD.prepare_single_phase_liquid(T_inlet=T_inlet,steps=steps, p_ref = p_inlet, m_dot=m_dot, fp=fp)
        # set the Nusselt function
        Nu_func_liquid = thermo.convection.Nu_DB # [-] Nusselt number for liquid phase
        # Set remaining inputs
        T_wall = 500 # [K] Wall temperature\

        #Geometry
        w_h = 100e-6 # [m] Width and height
        A_channel = w_h**2 # [m^2]
        wetted_perimeter = 4*w_h # [m]
        D_hydr = 4*A_channel/wetted_perimeter # [m]

        self.res = oneD.calc_channel_single_phase(\
            T = self.prep['T'],
            T_ref = self.prep['T_ref'],
            Q_dot= self.prep['Q_dot'],
            rho = self.prep['rho'],
            Pr = self.prep['Pr'],
            kappa = self.prep['kappa'],
            mu = self.prep['mu'],
            p_ref=p_inlet,
            m_dot=m_dot,
            T_wall=T_wall,
            D_hydr=D_hydr,
            wetted_perimeter=wetted_perimeter,
            A_channel=A_channel,
            Nu_func=Nu_func_liquid,
            fp=fp
            )
        return super().setUp()

    def testLength(self):
        # Results based on manual calculation in excelsheet Verification_one_d.xlsx
        exp_L0 = 0 # [m]
        exp_L1 = 2.709706792e-3 # [m]
        exp_L2 = 5.840860059e-3 # [m]

        self.assertEqual(exp_L0, self.res['L'][0])
        self.assertAlmostEqual(exp_L1, self.res['L'][1], delta=2e-3*exp_L1)
        self.assertAlmostEqual(exp_L2, self.res['L'][2], delta=2e-3*exp_L2)

    def testReynolds(self):
        # Results based on manual calculation in excelsheet Verification_one_d.xlsx
        exp_Re1 = 2.60E+03 # [-]
        exp_Re2 = 4.32E+03 # [-]

        self.assertAlmostEqual(exp_Re1, self.res['Re'][1], delta=1e-3*exp_Re1)
        self.assertAlmostEqual(exp_Re2, self.res['Re'][2], delta=1e-3*exp_Re2)

    def testNusselt(self):
        # Results based on manual calculation in excelsheet Verification_one_d.xlsx
        exp_Nu1 = 1.7665E+01 # [-]
        exp_Nu2 = 2.1533E+01 # [-]

        self.assertAlmostEqual(exp_Nu1, self.res['Nu'][1], delta=5e-3*exp_Nu1)
        self.assertAlmostEqual(exp_Nu2, self.res['Nu'][2], delta=5e-3*exp_Nu2)

    def testVelocity(self):
        # Results based on manual calculation in excelsheet Verification_one_d.xlsx
        exp_u1 = 1.024842E+01 # [-]
        exp_u2 = 1.060513E+01 # [-]

        self.assertAlmostEqual(exp_u1, self.res['u'][1], delta=1e-5*exp_u1)
        self.assertAlmostEqual(exp_u2, self.res['u'][2], delta=1e-5*exp_u2)

class TestPrepareSinglePhaseGas(unittest.TestCase):
    def setUp(self):
        # Inputs
        fp = FluidProperties('water')
        T_outlet = 450 # [K]
        p_inlet = 2e5  # [Pa] Inlet pressure
        steps = 10 # [-] Number of data point
        m_dot = 1e-3 # [kg/s]

        self.prep = oneD.prepare_single_phase_gas(T_outlet=T_outlet,steps=steps, p_ref = p_inlet, m_dot=m_dot, fp=fp)
        return super().setUp()
    
    def testSaturationValues(self):
        # Reference: https://webbook.nist.gov/cgi/fluid.cgi?Action=Load&ID=C7732185&Type=SatT&Digits=5&PLow=2&PHigh=2&PInc=1&RefState=DEF&TUnit=K&PUnit=bar&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm
        # Saturation values which must be true for the FIRST element in the array
        exp_T = 393.36 # [K] Saturation temperature
        exp_rho = 1.1291 # [kg/m^3] Density
        exp_h = 2706.2e3 # [J/kg] Enthalpy
        exp_mu = 1.2963e-05 # [Pa*s] Viscosity
        exp_kappa = 0.027493 # [W/(m*K)] Thermal conductivity
        exp_Pr = 1.0270253009857053 # [-] Prandtl

        self.assertAlmostEqual(exp_T, self.prep['T'][0], delta=exp_T*1e-5)
        self.assertAlmostEqual(exp_rho, self.prep['rho'][0], delta=exp_rho*1e-4)
        self.assertAlmostEqual(exp_h, self.prep['h'][0], delta=exp_h*1e-4)
        self.assertAlmostEqual(exp_mu, self.prep['mu'][0], delta=exp_mu*2.5e-3)
        self.assertAlmostEqual(exp_kappa, self.prep['kappa'][0], delta=exp_kappa*3e-2)
        self.assertAlmostEqual(exp_Pr, self.prep['Pr'][0], delta=exp_Pr*3e-2)

    def testTotalQDot(self):
        # The total power put in must of course equal the specific enthalpy change times mass flow
        # Reference: https://webbook.nist.gov/cgi/fluid.cgi?Action=Load&ID=C7732185&Type=IsoBar&Digits=5&P=2&THigh=500&TLow=300&TInc=100&RefState=DEF&TUnit=K&PUnit=bar&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm
        exp_total_Q_dot = 117.8 # [W] 
        res_total_Q_dot = np.sum(self.prep['Q_dot'])

        self.assertAlmostEqual(exp_total_Q_dot, res_total_Q_dot, delta=exp_total_Q_dot*5e-4)

    def testThirdValue(self):
        # Just to check if values are proper at intermediate sections
        # There are ten values in the array, so 9 steps from T_inlet to T_sat
        #  so the third value should 2/9 along the between saturation and inlet values
        exp_T = (450-393.36)*(2/9) + 393.36 # [K] 405.95

        self.assertAlmostEqual(exp_T, self.prep['T'][2], delta=exp_T*1e-5)
        # This expected value can be used as reference to check other values at this point
        # Ref: https://webbook.nist.gov/cgi/fluid.cgi?Action=Load&ID=C7732185&Type=IsoBar&Digits=5&P=2&THigh=450&TLow=393.36&TInc=6.293333333333332&RefState=DEF&TUnit=K&PUnit=bar&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm
        # NOTE: Slight difference in reference value as webbook rounds temperature to 320.75 K
        exp_rho = 1.0901 # [kg/m^3]
        exp_h = 2733.2e3 # [J/kg]
        exp_mu = 1.3454e-05 # [Pa*s]
        exp_kappa = 0.028312 # [W/(m*K)]
        exp_Pr = 1.004250430912687 # [-]

        self.assertAlmostEqual(exp_rho, self.prep['rho'][2], delta=exp_rho*1e-4)
        self.assertAlmostEqual(exp_h, self.prep['h'][2], delta=exp_h*1e-4)
        self.assertAlmostEqual(exp_mu, self.prep['mu'][2], delta=exp_mu*2e-3)
        self.assertAlmostEqual(exp_kappa, self.prep['kappa'][2], delta=exp_kappa*2.5e-2)
        self.assertAlmostEqual(exp_Pr, self.prep['Pr'][2], delta=exp_Pr*2.5e-2)
        
    def testDeltaT(self):
        # Check if delta T is evenly spaced as expected
        exp_dT = (450-393.36)/9 # [K] temperature step
        self.assertAlmostEqual(exp_dT, self.prep['dT'], delta=1e-5*exp_dT)

    def testReferenceTemperature(self):
        # Check if if T_ref_i is indeed half-way between T_i-1 and T_i ()
        exp_T_ref = 393.36 + (450-393.36)*(1/9)/2
        self.assertAlmostEqual(exp_T_ref, self.prep['T_ref'][1],delta=exp_T_ref*1e-5)

class TestPrepareHomogeneousTransition(unittest.TestCase):
    def setUp(self):
        # Inputs
        fp = FluidProperties('water')
        p_inlet = 2e5  # [Pa] Inlet pressure
        steps = 10 # [-] Number of data point
        m_dot = 1e-4 # [kg/s]

        self.prep = oneD.prepare_homogenous_transition(steps=steps, p = p_inlet, m_dot=m_dot, fp=fp)
        return super().setUp()
    
    def testSaturationValues(self):
        # Reference: https://webbook.nist.gov/cgi/fluid.cgi?Action=Load&ID=C7732185&Type=SatT&Digits=5&PLow=2&PHigh=2&PInc=1&RefState=DEF&TUnit=K&PUnit=bar&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm
        # Saturation values which must be true for the FIRST element in the array
        exp_T_sat = 393.36 # [K] Saturation temperature
        exp_rho_gas = 1.1291 # [kg/m^3] Density
        exp_h_gas = 2706.2e3 # [J/kg] Enthalpy
        exp_mu_gas = 1.2963e-05 # [Pa*s] Viscosity
        exp_kappa_gas = 0.027493 # [W/(m*K)] Thermal conductivity
        exp_Pr_gas = 1.0270253009857053 # [-] Prandtl

        self.assertAlmostEqual(exp_T_sat, self.prep['T_sat'], delta=exp_T_sat*1e-5)
        self.assertAlmostEqual(exp_rho_gas, self.prep['rho'][-1], delta=exp_rho_gas*1e-4)
        self.assertAlmostEqual(exp_h_gas, self.prep['h'][-1], delta=exp_h_gas*1e-4)
        self.assertAlmostEqual(exp_mu_gas, self.prep['mu'][-1], delta=exp_mu_gas*2.5e-3)
        self.assertAlmostEqual(exp_kappa_gas, self.prep['kappa'][-1], delta=exp_kappa_gas*3e-2)
        self.assertAlmostEqual(exp_Pr_gas, self.prep['Pr'][-1], delta=exp_Pr_gas*3e-2)

        # Reference: https://webbook.nist.gov/cgi/fluid.cgi?Action=Load&ID=C7732185&Type=SatT&Digits=5&PLow=2&PHigh=2&PInc=1&RefState=DEF&TUnit=K&PUnit=bar&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm
        # Saturation values which must be true for the FIRST element in the array
        exp_rho_liquid = 942.94 # [kg/m^3] Density
        exp_h_liquid = 504.70e3 # [J/kg] Enthalpy
        exp_mu_liquid = 0.00023162 # [Pa*s] Viscosity
        exp_kappa_liquid = 0.68321 # [W/(m*K)] Thermal conductivity
        exp_Pr_liquid = 1.438755460253802 # [-] Prandtl

        self.assertAlmostEqual(exp_rho_liquid, self.prep['rho'][0], delta=exp_rho_liquid*1e-5)
        self.assertAlmostEqual(exp_h_liquid, self.prep['h'][0], delta=exp_h_liquid*1e-5)
        self.assertAlmostEqual(exp_mu_liquid, self.prep['mu'][0], delta=exp_mu_liquid*1e-3)
        self.assertAlmostEqual(exp_kappa_liquid, self.prep['kappa'][0], delta=exp_kappa_liquid*2e-3)
        self.assertAlmostEqual(exp_Pr_liquid, self.prep['Pr'][0], delta=exp_Pr_liquid*2e-3)

    def testTotalQDot(self):
        # The total power put in must of course equal the specific enthalpy change times mass flow
        # Reference: https://webbook.nist.gov/cgi/fluid.cgi?Action=Load&ID=C7732185&Type=IsoBar&Digits=5&P=2&THigh=500&TLow=300&TInc=100&RefState=DEF&TUnit=K&PUnit=bar&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm
        exp_total_Q_dot = 220.15 # [W] 
        res_total_Q_dot = np.sum(self.prep['Q_dot'])

        self.assertAlmostEqual(exp_total_Q_dot, res_total_Q_dot, delta=exp_total_Q_dot*5e-4)

    def testThirdValue(self):
        # Just to check if values are proper at intermediate sections
        # There are ten values in the array, so 9 steps from x=0 to x=1
        #  The third value depends both on void fraction alpha and vapour quality
        # Results based on manual calculation in excelsheet Verification_one_d.xlsx
        exp_x = 2/9 # [-] Vapour quality
        exp_alpha = 0.995826503 # [-] Void fraction

        self.assertAlmostEqual(exp_x, self.prep['x'][2], delta=exp_x*1e-5)
        self.assertAlmostEqual(exp_alpha, self.prep['alpha'][2], delta=exp_alpha*1e-7)
        # This expected value can be used as reference to check other values at this point
        # Ref: https://webbook.nist.gov/cgi/fluid.cgi?Action=Load&ID=C7732185&Type=SatT&Digits=5&PLow=2&PHigh=2&PInc=1&RefState=DEF&TUnit=K&PUnit=bar&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm
        # Calculating in steps of 2/9 from the data points in the reference
        
        
        exp_rho = 5.059744672 # [kg/m^3]
        exp_h = 993.9222222e3 # [J/kg]
        exp_mu = 1.3876E-05 # [Pa*s]
        exp_kappa =0.030229633 # [W/(m*K)]
        exp_Pr = 1.028743655# [-]

        self.assertAlmostEqual(exp_rho, self.prep['rho'][2], delta=exp_rho*1e-4)
        self.assertAlmostEqual(exp_h, self.prep['h'][2], delta=exp_h*1e-4)
        self.assertAlmostEqual(exp_mu, self.prep['mu'][2], delta=exp_mu*2.5e-3)
        self.assertAlmostEqual(exp_kappa, self.prep['kappa'][2], delta=exp_kappa*3e-2)
        self.assertAlmostEqual(exp_Pr, self.prep['Pr'][2], delta=exp_Pr*3e-2)

class TestCalcHomogeneousTransition(unittest.TestCase):
    def setUp(self):
        # Inputs
        fp = FluidProperties('water')
        p_inlet = 2e5  # [Pa] Inlet pressure
        steps = 3 # [-] Number of data point
        m_dot = 0.1e-3 # [kg/s]
        T_wall = 500 # [K]

        self.prep = oneD.prepare_homogenous_transition(steps=steps, p = p_inlet, m_dot=m_dot, fp=fp)
        
        # Nusselt functions
        Nu_func_two_phase = tp.Nu_Kandlikar_NBD_dryout # [-] Function to calculate Nusselt number (two-phase)
        Nu_func_le = thermo.convection.Nu_DB # [-] Function to calculate Nusselt in two-phase, AS IF the flow was entirely liquid (two-phase, le)
        # NOTE: Can be set to one, so only the last value should change from 0 to 1. As it allows for better testing
        Nu_func_dryout = thermo.two_phase.Nu_DB_two_phase #thermo.two_phase.Nu_DB_two_phase # [-] Function to calculate the Nusselt number after dry-out. It is up to Nu_func_two_phase to decide if/how to apply it

        #Geometry
        w_h = 100e-6 # [m] Width and height
        A_channel = w_h**2 # [m^2]
        wetted_perimeter = 4*w_h # [m]
        D_hydr = 4*A_channel/wetted_perimeter # [m]

        x_tp = self.prep['x']
        alpha_tp = self.prep['alpha']
        T_sat = self.prep['T_sat']
        rho_tp_l = self.prep['rho_l']
        rho_tp_g = self.prep['rho_g']
        rho_tp = self.prep['rho']
        mu_tp_l = self.prep['mu_l']
        mu_tp = self.prep['mu']
        Pr_tp_l = self.prep['Pr_l']
        Pr_tp = self.prep['Pr']
        kappa_tp_l = self.prep['kappa_l']
        kappa_tp = self.prep['kappa']
        Q_dot_tp = self.prep['Q_dot']

        self.res = oneD.calc_homogenous_transition(
            p_sat=p_inlet,
            x=x_tp,
            alpha=alpha_tp,
            T_sat=T_sat,
            rho_l=rho_tp_l,
            rho_g=rho_tp_g,
            rho=rho_tp,
            m_dot=m_dot,
            mu_l=mu_tp_l,
            mu=mu_tp,
            Pr_l=Pr_tp_l,
            Pr=Pr_tp,
            kappa_l=kappa_tp_l,
            kappa=kappa_tp,
            Q_dot=Q_dot_tp,
            T_wall=T_wall,
            D_hydr=D_hydr,
            wetted_perimeter=wetted_perimeter,
            A_channel=A_channel,
            Nu_func_tp=Nu_func_two_phase,
            Nu_func_le=Nu_func_le,
            Nu_func_dryout=Nu_func_dryout,
            fp=fp
        )

        return super().setUp()

    def testVelocity(self):
        # Results based on manual calculation in excelsheet Verification_one_d.xlsx
        exp_u0 = 10.6051 # [m/s]
        exp_u1 = 4433.61 # [m/s]
        exp_u2 = 8856.61 # [m/s]

        self.assertAlmostEqual(exp_u0, self.res['u'][0], delta=1e-4*exp_u0)
        self.assertAlmostEqual(exp_u1, self.res['u'][1], delta=1e-4*exp_u1)
        self.assertAlmostEqual(exp_u2, self.res['u'][2], delta=1e-4*exp_u2)

    def testReynolds(self):
        # Results based on manual calculation in excelsheet Verification_one_d.xlsx
        exp_Re0 = 4.3174E+03 # [m/s]
        exp_Re1 = 7.5617E+04 # [m/s]
        exp_Re2 = 7.7143E+04 # [m/s]

        self.assertAlmostEqual(exp_Re0, self.res['Re'][0], delta=1e-3*exp_Re0)
        self.assertAlmostEqual(exp_Re1, self.res['Re'][1], delta=2.5e-3*exp_Re1)
        self.assertAlmostEqual(exp_Re2, self.res['Re'][2], delta=2.5e-3*exp_Re2)

    def testNusselt(self):
        # Results based on manual calculation in excelsheet Verification_one_d.xlsx
        exp_Nu0 = 2.6043E+02 # [m/s]
        exp_Nu1 = 1.6577E+02 # [m/s]
        # Dry-out function is set to 1 to only change the last Nusselt number, to allow for better testing of Kandlikar function
        exp_Nu2 = 7.6009E+00 # [m/s]

        self.assertAlmostEqual(exp_Nu0, self.res['Nu'][0], delta=1e-3*exp_Nu0)
        self.assertAlmostEqual(exp_Nu1, self.res['Nu'][1], delta=2.5e-3*exp_Nu1)
        self.assertAlmostEqual(exp_Nu2, self.res['Nu'][2], delta=2e-2*exp_Nu2)

    def testNusseltLiquidEntire(self):
        # Results based on manual calculation in excelsheet Verification_one_d.xlsx
        # The Nusselt number for the case as if the flow was entirely liquid
        exp_Nu_le = 2.1533E+01 # [-] In this case Dittus Boelter, as seen in setup.

        self.assertAlmostEqual(exp_Nu_le, self.res['Nu_le'], delta=1e-3*exp_Nu_le)

    def testLength(self):
         # Results based on manual calculation in excelsheet Verification_one_d.xlsx
        exp_L0 = 0 # [m]
        exp_L1 = 2.2785E-03 # [m]
        exp_L2 = 5.1971E-02 # [m]

        self.assertEqual(exp_L0, self.res['L'][0])
        self.assertAlmostEqual(exp_L1, self.res['L'][1], delta=1e-3*exp_L1)
        self.assertAlmostEqual(exp_L2, self.res['L'][2], delta=2e-2*exp_L2)
        # So, 2 percent difference in length calculations