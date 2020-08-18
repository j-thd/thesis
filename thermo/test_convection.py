import unittest
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
        u=1e-3

        exp_Nu = 0.0018785745665208552
        res_Nu = thermo.convection.Nusselt_Dittus_Boelter(T_inlet=T_inlet,T_outlet=T_outlet,p=p, D_hydraulic=D_h,L_channel=L,u=u,fp=fp, supressExceptions=True)
        self.assertAlmostEqual(exp_Nu,res_Nu,delta=0.00001*exp_Nu/res_Nu)

        # Same but with cooling assumption

        exp_Nu = 0.001896461381677763
        res_Nu = thermo.convection.Nusselt_Dittus_Boelter(T_inlet=T_inlet,T_outlet=T_outlet,p=p, D_hydraulic=D_h,L_channel=L,u=u,fp=fp, heating=False, supressExceptions=True)
        self.assertAlmostEqual(exp_Nu,res_Nu,delta=0.00001*exp_Nu/res_Nu)