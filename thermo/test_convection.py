import unittest
from unittest.case import TestCase
from thermo.convection import Nu_DB
from thermo.prop import FluidProperties
import thermo.convection

class TestNusseltDittusBoelter(unittest.TestCase):
    def test_one(self):
        # Calculate a known case manually, and check for similarity
        # Source: https://webbook.nist.gov/cgi/fluid.cgi?Action=Load&ID=C7732185&Type=IsoBar&Digits=5&P=0.5&THigh=1000&TLow=300&TInc=100&RefState=DEF&TUnit=K&PUnit=MPa&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm

        fp = FluidProperties("HEOS::Water")
        T_inlet = 700 # [K]
        T_outlet = 900 # [K]
        # This makes the bulk temperature 800 K
        p = 5e5 # [Pa]

        D_h = 1e-3
        L = 1000e-3
        
        m_dot = 1e-6 # [kg/s]
        A= 0.0007359976448075365 # [m^2]

        exp_Nu = 0.0018785745665208552
        res_Nu = thermo.convection.Nusselt_Dittus_Boelter(T_inlet=T_inlet,T_outlet=T_outlet,p=p, D_hydraulic=D_h,L_channel=L,m_dot=m_dot,A=A,fp=fp, supressExceptions=True)
        self.assertAlmostEqual(exp_Nu,res_Nu,delta=0.01*exp_Nu)

        # Same but with cooling assumption

        exp_Nu = 0.001896461381677763
        res_Nu = thermo.convection.Nusselt_Dittus_Boelter(T_inlet=T_inlet,T_outlet=T_outlet,p=p, D_hydraulic=D_h,L_channel=L,m_dot=m_dot, A=A,fp=fp, heating=False, supressExceptions=True)
        self.assertAlmostEqual(exp_Nu,res_Nu,delta=0.01*exp_Nu)

class TestNuDB_func(unittest.TestCase):
    def test_one(self):
        # Simply re-using case above, only heating
        fp = FluidProperties("HEOS::Water")
        T_inlet = 700 # [K]
        T_outlet = 900 # [K]
        # This makes the bulk temperature 800 K
        p = 5e5 # [Pa]

        D_h = 1e-3
        u=1e-3

        T_bulk = (T_outlet+T_inlet)/2
        Re = fp.get_Reynolds_from_velocity(T=T_bulk,p=p,L_ref=D_h,u=u)
        Pr = fp.get_Prandtl(T=T_bulk,p=p)
        arguments = {   'fp':fp,
                        'Re': Re,
                        'Pr': Pr}
        exp_Nu = 0.0018785745665208552
        res_Nu = thermo.convection.Nu_DB(args=arguments)
        self.assertAlmostEqual(exp_Nu,res_Nu,delta=0.00001*exp_Nu/res_Nu)

class TestStantonFromNuFunc(unittest.TestCase):
    def test_one(self):
        #Simply Re-using case above again
        fp = FluidProperties("water")
        T_inlet = 700 # [K]
        T_outlet = 900 # [K]
        # This makes the bulk temperature 800 K
        p = 5e5 # [Pa]

        D_h = 1e-3
        m_dot = 1e-6 # [kg/s]
        A= 0.0007359976448075365 # [m^2]

        Nu_func = Nu_DB
        T_bulk = (T_inlet + T_outlet)/2

        Pr= 0.9095872167254094
        Re = 0.04578292954139569

        exp_St = 0.0018785745665208552/(Re*Pr)
        res_St = thermo.convection.Stanton_from_Nusselt_func_and_velocity(Nu_func=Nu_func, m_dot=m_dot, A=A, T_ref=T_bulk, p_ref=p, L_ref=D_h, fp=fp)
        self.assertAlmostEqual(exp_St,res_St,delta=exp_St*0.01)

        def dummy_func(args):
            return 5

        # Another test with a dummy func
        exp_St = 5/(Re*Pr)
        res_St = thermo.convection.Stanton_from_Nusselt_func_and_velocity(Nu_func=dummy_func, m_dot=m_dot,A=A, T_ref=T_bulk, p_ref=p, L_ref=D_h, fp=fp)
        self.assertAlmostEqual(exp_St,res_St,delta=exp_St*0.01)

class TestNuKandlikar_NDB_Re_low_sat_gas_constant_wall_temp_square_water(unittest.TestCase):
    def test_simple_inputs(self):
        args = {'Bo': 1}
        exp_Nu = 3823.612

        res_Nu = thermo.convection.Nu_Kandlikar_NDB_Re_low_sat_gas_constant_wall_temp_square_water(args=args)
        self.assertAlmostEqual(exp_Nu, res_Nu, delta=1e-10*exp_Nu)
        
        args = {'Bo': 0.5}
        exp_Nu = 2353.709276299291
        
        res_Nu = thermo.convection.Nu_Kandlikar_NDB_Re_low_sat_gas_constant_wall_temp_square_water(args=args)
        self.assertAlmostEqual(exp_Nu, res_Nu, delta=1e-10*exp_Nu)
