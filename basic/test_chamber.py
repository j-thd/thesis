import unittest
import basic.chamber as chamber
import thermo.prop 
from thermo.prop import FluidProperties

class TestIdealEnthalpyChange(unittest.TestCase):
    def test_from_table(self):
        # Set water to be the fluid
        fp = thermo.prop.FluidProperties("HEOS::Water")
        # Test a few cases by calculating enthalpy differences from table
        # T1, p1, T2, p2, delta_h, places
        in_out = (
        (300, 0.5e6, 1000, 0.5e6, 3875.448e3, -1),
        (300, 0.5e6, 400, 0.5e6, 420.105e3, 0),
        (330, 0.25e6, 410, 0.25e6, 2498.946e3,-1)
        )

        for T_1, p_1, T_2, p_2, delta_h, places in in_out:
            res = chamber.ideal_enthalpy_change(T_inlet=T_1, p_inlet=p_1, T_outlet=T_2, p_outlet=p_2, fp=fp)
            self.assertAlmostEqual(res, delta_h, places)

class TestIdealPowerConsumption(unittest.TestCase):
    def test_from_table(self):

        # Set water to be the fluid
        fp = thermo.prop.FluidProperties("HEOS::Water")

        # Same as previous test, but with mass flow simply added
        in_out = (
        (0, 300, 0.5e6, 1000, 0.5e6, 3875.448e3*0, -1),
        (1, 300, 0.5e6, 400, 0.5e6, 420.105e3, 0),
        (2, 330, 0.25e6, 410, 0.25e6, 2498.946e3*2,-2)
        )

        for m, T_1, p_1, T_2, p_2, delta_h, places in in_out:
            res = chamber.ideal_power_consumption(mass_flow=m, T_inlet=T_1, p_inlet=p_1, T_outlet=T_2, p_outlet=p_2, fp=fp)
            self.assertAlmostEqual(res, delta_h, places)

class TestIdealHeaterTemperature(unittest.TestCase):
    def test_one_input(self):
        T_expected = 2
        res = chamber.ideal_heater_temperature(P_mh=1, T_inlet=1,T_outlet=1,A=1,k=1,d=1)
        self.assertEqual(T_expected,res)
        T_expected = 106
        res = chamber.ideal_heater_temperature(P_mh=1, T_inlet=50,T_outlet=150,A=2,k=3,d=36)
        self.assertEqual(T_expected,res)
        T_expected = 1
        res = chamber.ideal_heater_temperature(P_mh=1, T_inlet=1,T_outlet=1,A=1,k=1,d=0)
        self.assertEqual(T_expected,res)

class TestConvectiveHeatFlow(unittest.TestCase):
    def test_simple_input(self):
        # Simple inputs to check mathematical behaviour is correct
        # It is expected that when heat flows from station 1 to station 2 that this is NEGATIVE heat flow
        A = 1 # [m^2] Area
        T1 = 20 # [K] Heat at station 1 
        T2 = 10 # [K] Heat at station 2
        h_conv = 2 # [W/(m^2*K)] Convective heat transfer coefficient

        exp_Q = -20 # [W] Expected heat flow
        res_Q = chamber.convective_heat_flow(heat_transfer_coefficient=h_conv,T_wall=T1,T_ref=T2,A_wall=A)
        self.assertEqual(exp_Q,res_Q)

        A = 0.5
        exp_Q = -10 # [W] Expected heat flow
        res_Q = chamber.convective_heat_flow(heat_transfer_coefficient=h_conv,T_wall=T1,T_ref=T2,A_wall=A)
        self.assertEqual(exp_Q,res_Q)

        T2 = 40
        exp_Q = 20 # [W] Expected heat flow
        res_Q = chamber.convective_heat_flow(heat_transfer_coefficient=h_conv,T_wall=T1,T_ref=T2,A_wall=A)
        self.assertEqual(exp_Q,res_Q)

        h_conv = 30
        exp_Q = 300 # [W] Expected heat flow
        res_Q = chamber.convective_heat_flow(heat_transfer_coefficient=h_conv,T_wall=T1,T_ref=T2,A_wall=A)
        self.assertEqual(exp_Q,res_Q)

        h_conv = -30
        with self.assertRaises(AssertionError):
            res_Q = chamber.convective_heat_flow(heat_transfer_coefficient=h_conv,T_wall=T1,T_ref=T2,A_wall=A)

class TestHConvFromStanton(unittest.TestCase):
    def test_water_input(self):
        # Source: https://webbook.nist.gov/cgi/fluid.cgi?Action=Load&ID=C7732185&Type=IsoBar&Digits=5&P=0.1&THigh=1000&TLow=300&TInc=100&RefState=DEF&TUnit=K&PUnit=MPa&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm
        fp = FluidProperties("HEOS::Water")

        T = 500 # [K]
        p = 5e5 # [Pa]
        u = 1e-3 # [m/s]
        St = 1

        exp_h_conv = 4.6448084 # [W/(m^2*K)] Convective heat transfer
        res_h_conv = chamber.h_conv_from_Stanton(Stanton=St,T_ref=T,p_ref=p,u=u,fp=fp)
        self.assertAlmostEqual(exp_h_conv,res_h_conv,delta=exp_h_conv/res_h_conv*0.0001)



