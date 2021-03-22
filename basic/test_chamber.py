import unittest
import basic.chamber as chamber
import thermo.prop 
from thermo.prop import FluidProperties
import numpy as np

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

class TestWettedPerimeterRectangular(unittest.TestCase):
    def test_simple_input(self):
        h = 1
        w = 1

        exp = 4
        res = chamber.wetter_perimeter_rectangular(w_channel=w, h_channel=h)
        self.assertEqual(res,exp)

        h=4

        exp = 10
        res = chamber.wetter_perimeter_rectangular(w_channel=w, h_channel=h)
        self.assertEqual(res,exp)

        w = 2
        exp= 12
        res = chamber.wetter_perimeter_rectangular(w_channel=w, h_channel=h)
        self.assertEqual(res,exp)

class TestHydraulicDiameterRectangular(unittest.TestCase):
    def test_simple_input(self):
        h=1
        w=1

        # Should revert to side's length for a square channel
        exp = 1
        res = chamber.hydraulic_diameter_rectangular(w_channel=w,h_channel=h)
        self.assertEqual(exp,res)

        h=2
        exp= 4*2/6
        res = chamber.hydraulic_diameter_rectangular(w_channel=w,h_channel=h)
        self.assertEqual(exp,res)

        w = 4
        exp = 4*8/12
        res = chamber.hydraulic_diameter_rectangular(w_channel=w,h_channel=h)
        self.assertEqual(exp,res)



class TestVelocityFromMassFlow(unittest.TestCase):
    def test_simple_input(self):
        A=1
        m_dot=1
        rho=1

        exp = 1
        res = chamber.velocity_from_mass_flow(A=A, m_dot=m_dot, rho=rho) 
        self.assertEqual(exp,res)

        A=2
        exp=0.5
        res = chamber.velocity_from_mass_flow(A=A, m_dot=m_dot, rho=rho) 
        self.assertEqual(exp,res)

        rho= 5
        exp=0.1
        res = chamber.velocity_from_mass_flow(A=A, m_dot=m_dot, rho=rho) 
        self.assertEqual(exp,res)

        m_dot=100
        exp=10
        res = chamber.velocity_from_mass_flow(A=A, m_dot=m_dot, rho=rho) 
        self.assertEqual(exp,res)

class TestMassFlow(unittest.TestCase):
    def test_simple_input(self):
        rho=3
        A=5
        u=7

        exp = 105
        res = chamber.mass_flow(A=A, u=u, rho=rho)
        self.assertEqual(exp,res)

        u =1
        exp=15
        res = chamber.mass_flow(A=A, u=u, rho=rho)
        self.assertEqual(exp,res)

        A = 10
        exp = 30
        res = chamber.mass_flow(A=A, u=u, rho=rho)
        self.assertEqual(exp,res)
        
        rho = 1
        exp = 10
        res = chamber.mass_flow(A=A, u=u, rho=rho)
        self.assertEqual(exp,res)


class TestRadiationLoss(unittest.TestCase):
    def test_simple_input(self):
        # inputs of 1 should return stefan-boltzmann constant
        emmisivity = 1
        A = 1
        T = 1

        exp = 5.670374419e-8 # [W/(m^2*K^4)]  Source: https://en.wikipedia.org/wiki/Stefan%E2%80%93Boltzmann_constant
        res = chamber.radiation_loss(T=T,A=A, emmisivity=emmisivity)
        self.assertEqual(exp,res)

        emmisivity = 0.5
        exp = 0.5*exp
        res = chamber.radiation_loss(T=T,A=A, emmisivity=emmisivity)
        self.assertEqual(exp,res)

        A = 4
        exp = 4*exp
        res = chamber.radiation_loss(T=T,A=A, emmisivity=emmisivity)
        self.assertEqual(exp,res)

        T = 2
        exp= 16*exp
        res = chamber.radiation_loss(T=T,A=A, emmisivity=emmisivity)
        self.assertEqual(exp,res)

class TestRequiredPower(unittest.TestCase):
    def test_simple_input(self):
        m_dot = 1
        delta_h = 1

        exp = 1
        res = chamber.required_power(m_dot=m_dot, delta_h=delta_h)
        self.assertEqual(exp,res)

        m_dot = 2
        exp = 2
        res = chamber.required_power(m_dot=m_dot, delta_h=delta_h)
        self.assertEqual(exp,res)

        delta_h = 0.5
        exp=1
        res = chamber.required_power(m_dot=m_dot, delta_h=delta_h)
        self.assertEqual(exp,res)

class TestRequiredHeaterArea(unittest.TestCase):
    def test_simple_input(self):
        h=1
        Q=1 
        T_wall = 2
        T_ref = 1

        exp= 1
        res = chamber.required_heater_area(Q_dot=Q,h_conv=h,T_wall=T_wall,T_ref=T_ref)
        self.assertEqual(exp,res)

        T_wall = 3
        exp = 0.5
        res = chamber.required_heater_area(Q_dot=Q,h_conv=h,T_wall=T_wall,T_ref=T_ref)
        self.assertEqual(exp,res)
        
        h=0.25
        exp = 2
        res = chamber.required_heater_area(Q_dot=Q,h_conv=h,T_wall=T_wall,T_ref=T_ref)
        self.assertEqual(exp,res)

        Q=10
        exp = 20
        res = chamber.required_heater_area(Q_dot=Q,h_conv=h,T_wall=T_wall,T_ref=T_ref)
        self.assertEqual(exp,res)
        
        T_wall = 5
        exp = 10
        res = chamber.required_heater_area(Q_dot=Q,h_conv=h,T_wall=T_wall,T_ref=T_ref)
        self.assertEqual(exp,res)

        T_ref = 4
        exp = 40
        res = chamber.required_heater_area(Q_dot=Q,h_conv=h,T_wall=T_wall,T_ref=T_ref)
        self.assertEqual(exp,res)

class TestEnthalpyDifferencePerSection(unittest.TestCase):
    def test_simple_input(self):
        a= np.array((0,1,2,3,4))
        exp = np.array((0,1,1,1,1))
        res = chamber.delta_enthalpy_per_section(a)
        np.testing.assert_equal(exp,res)

        a= np.array((1,1,2,3,5))
        exp = np.array((0,0,1,1,2))
        res = chamber.delta_enthalpy_per_section(a)
        np.testing.assert_equal(exp,res)

class TestReynoldsFromMassFlow(unittest.TestCase):
    def test_simple_input(self):
        m_dot = 1 # [kg/s]
        A = 1 # [kg/s] 
        mu = 1e-6 # [kg/s]
        D_h = 1e-3 # [kg/s]

        exp = 1e3
        res = chamber.Reynolds_from_mass_flow(m_dot=m_dot,A_channel=A, L_ref=D_h,mu=mu)
        self.assertAlmostEqual(exp,res,delta=1e-13*exp)

        m_dot = 2
        exp = 2e3
        res = chamber.Reynolds_from_mass_flow(m_dot=m_dot,A_channel=A, L_ref=D_h,mu=mu)
        self.assertAlmostEqual(exp,res,delta=1e-13*exp)

        A = 0.5
        exp = 4e3
        res = chamber.Reynolds_from_mass_flow(m_dot=m_dot,A_channel=A, L_ref=D_h,mu=mu)
        self.assertAlmostEqual(exp,res,delta=1e-13*exp)

        D_h = 1
        exp = 4e6
        res = chamber.Reynolds_from_mass_flow(m_dot=m_dot,A_channel=A, L_ref=D_h,mu=mu)
        self.assertAlmostEqual(exp,res,delta=1e-13*exp)
        
        mu = 1e-9
        exp = 4e9
        res = chamber.Reynolds_from_mass_flow(m_dot=m_dot,A_channel=A, L_ref=D_h,mu=mu)
        self.assertAlmostEqual(exp,res,delta=1e-13*exp)




