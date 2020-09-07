import unittest

import thermo.prop

class TestEnthalpy(unittest.TestCase):

    def testTable(self):
        # Test with table from IAPWS 95 report (Wagner2002)
        
        # Set water to be the fluid
        fp = thermo.prop.FluidProperties("HEOS::Water")

        pressure = 1e5 # 
        # T - h
        in_out = (
            (300, 112.654e3, 0),
            (330, 238.067e3, 0),
            (370, 405.891e3, 0),
            (375, 2679.60e3, -1),
            (385, 2700.11e3, -1),
            (600, 3128.76e3, -1),
            (725, 3386.73e3, -1)
        )

        for T, h, places in in_out:
            res_h = fp.get_enthalpy(T=T,p=pressure)
            self.assertAlmostEqual(res_h, h, places=places)

        # again for a different pressure
        pressure = 1e6
        # T- h
        in_out = (
            (275, 8.768e3, 0),
            (290, 71.678e3, 0),
            (325, 217.926e3, 0),
            (450, 749.197e3, 0),
            (460, 2795.49e3, -1),
            (570, 3044.88e3, -1),
            (675, 3268.41e3, -1)
        )

        for T, h, places in in_out:
            res_h = fp.get_enthalpy(T=T,p=pressure)
            self.assertAlmostEqual(res_h, h, places=places)

class TestDensity(unittest.TestCase):
    def test_table(self):
        # Wagner2002 page 496 and onwards
        fp = thermo.prop.FluidProperties("HEOS::Water")
        # Set pressure to 1 bar
        p = 1e5 # [Pa]

        # T [K] in, rho [kg/m^3] out, places
        in_out = (  (300, 996.556, 3),
                    (335, 982.233, 3),
                    (370, 960.591, 3),
                    (375, 0.58653, 3),
                    (580, 0.37444, 3),
                    (750, 0.28915, 3),
                    (1000, 0.21673, 3),
                    (1250, 0.17335, 3))

        for T, rho, places in in_out:
            res_rho = fp.get_density(T=T,p=p) # [kg/m^3]
            self.assertAlmostEqual(rho, res_rho, places=places)

        # Same, but at 10 bar
        p = 10e5 # [Pa]

        # T [K] in, rho [kg/m^3] out, places
        in_out = (  (300, 996.960, 3),
                    (335, 982.627, 3),
                    (370, 961.010, 3),
                    (375, 957.434, 3),
                    (450, 890.386, 3),
                    (460, 5.0376, 3),
                    (1000, 2.1723, 3),
                    (1250, 1.7347, 3))

        for T, rho, places in in_out:
            res_rho = fp.get_density(T=T,p=p) # [kg/m^3]
            self.assertAlmostEqual(rho, res_rho, places=places)

class TestPressure(unittest.TestCase):
    def test_table(self):
        # Same as TestDensity but reversed
        # Wagner2002 page 496 and onwards
        fp = thermo.prop.FluidProperties("HEOS::Water")
        # Set pressure to 1 bar
        p = 1e5 # [Pa]

        # T [K] in, rho [kg/m^3] out, places
        in_out = (  (300, 996.556, -4),
                    (335, 982.233, -4),
                    (370, 960.591, -4),
                    (375, 0.58653, -4),
                    (580, 0.37444, -4),
                    (750, 0.28915, -4),
                    (1000, 0.21673, -4),
                    (1250, 0.17335, -4))

        for T, rho, places in in_out:
            res_p = fp.get_pressure(T=T,rho=rho) # [kg/m^3]
            self.assertAlmostEqual(p, res_p, places=places)

        # Same, but at 10 bar
        p = 10e5 # [Pa]

        # T [K] in, rho [kg/m^3] out, places
        in_out = (  (300, 996.960, -4),
                    (335, 982.627, -4),
                    (370, 961.010, -4),
                    (375, 957.434, -4),
                    (450, 890.386, -4),
                    (460, 5.0376, -4),
                    (1000, 2.1723, -4),
                    (1250, 1.7347, -4))

        
        for T, rho, places in in_out:
            res_p = fp.get_pressure(T=T,rho=rho) # [kg/m^3]
            self.assertAlmostEqual(p, res_p, places=places)

class TestCp(unittest.TestCase):
    def test_table(self):
        # Same as TestDensity but for cp
        # Wagner2002 page 496 and onwards
        fp = thermo.prop.FluidProperties("HEOS::Water")
        # Set pressure to 1 bar
        p = 1e5 # [Pa]

        # T [K] in, cp [kg/m^3] out, places
        in_out = (  (300, 4.1806e3, 1),
                    (335, 4.1858e3, 1),
                    (370, 4.2121e3, 1),
                    (375, 2.0686e3, 1),
                    (580, 2.0160e3, 1),
                    (750, 2.1191e3, 1),
                    (1000, 2.2921e3, 1),
                    (1250, 2.4632e3, 1))

        for T, cp, places in in_out:
            res_cp = fp.get_cp(T=T,p=p) # [kg/m^3]
            self.assertAlmostEqual(cp, res_cp, places=places)

        # Same, but at 10 bar
        p = 10e5 # [Pa]

        # T [K] in, rho [kg/m^3] out, places
        in_out = (  (300, 4.1781e3, 1),
                    (335, 4.1838e3, 1),
                    (370, 4.2101e3, 1),
                    (375, 4.2158e3, 1),
                    (450, 4.3924e3, 1),
                    (460, 2.5726e3, 1),
                    (1000, 2.3048e3, 1),
                    (1250, 2.4692e3, 1))

        
        for T, cp, places in in_out:
            res_cp = fp.get_cp(T=T,p=p) # [kg/m^3]
            self.assertAlmostEqual(cp, res_cp, places=places)

class TestSpecificHeatRatio(unittest.TestCase):
    def test_table(self):
        # Test reported values of gamma by taking cp/cv from Wagner2002 (page 496 and below)
        # Test if the table of the original paper can be reproduced Wagner2002
        # Page 496 of Wagner2002 
        fp = thermo.prop.FluidProperties("HEOS::Water")

        pressure = 1e5 # 1 bar of pressure

        # T in, gamma out
        # gamma not explicitly in table, but cp/cv is put in here instead
        in_out = (  (300, 1.01220279889593722, 2),
                    (330,  1.0479422889061443, 2),
                    (450, 1.3218229271230677, 3),
                    (600, 1.2992307692307692, 3),
                    (1000, 1.2527190249767721, 3))

        for T, gamma, places in in_out:
            res_gamma = fp.get_specific_heat_ratio(T=T,p=pressure) 
            self.assertAlmostEqual(res_gamma, gamma, places=places)

        # Again but with different pressure
        pressure = 1e6 # 10 bar of pressure

        # T in, gamma out
        # gamma not explicitly in table, but cp/cv is put in here instead
        in_out = (  (300, 1.012332816437294, 2),
                    (330,  1.0479675204250414, 2),
                    (450, 1.2890010564620262, 3),
                    (600, 1.3226487762454964, 3),
                    (1000, 1.2562956502779898, 3))

        for T, gamma, places in in_out:
            res_gamma = fp.get_specific_heat_ratio(T=T,p=pressure) 
            self.assertAlmostEqual(res_gamma, gamma, places=places)

class TestViscosity(unittest.TestCase):
    def test_table(self):
        # Test reported values of gamma by taking cp/cv from Wagner2002 (page 496 and below)
        # Test if the table of the original paper can be reproduced Wagner2002
        # Page 496 of Wagner2002 
        fp = thermo.prop.FluidProperties("HEOS::Water")

        pressure = 1e5 # 1 bar of pressure

        # T in, mu out, places
        in_out = (  (300, 0.00085383, 6),
                    (400, 1.3285e-05, 6),
                    (600, 2.1407e-05, 6),
                    (1000, 3.7592e-05, 6))

        for T, mu, places in in_out:
            res_mu = fp.get_viscosity(T=T,p=pressure) 
            self.assertAlmostEqual(res_mu, mu, places=places)


        # Again but with different pressure
        pressure = 1e6 # 10 bar of pressure
        
         # T in, mu out, places
        in_out = (  (300, 0.00085367, 6),
                    (400, 0.00021880, 6),
                    (600, 2.1329e-05, 6),
                    (1000, 3.7630e-05, 6))

        for T, mu, places in in_out:
            res_mu = fp.get_viscosity(T=T,p=pressure) 
            self.assertAlmostEqual(res_mu, mu, places=places)


        # Add new test at 2.6 bar, as it seems to return weird viscosities
        pressure = 0.26e6 

        in_out = (  (300, 853.80e-6),
                    (350, 368.83e-6),
                    (400, 218.60e-6),
                    (410, 13.579e-6),
                    (450, 15.184e-6),
                    (500, 17.231e-6),
                    (640, 23.064e-6),
                    (820, 30.484e-6))

        for T, mu in in_out:
            res_mu = fp.get_viscosity(T=T, p=pressure)
            self.assertAlmostEqual(res_mu, mu, delta=0.01*mu)
        
class TestPrandtl(unittest.TestCase):
    def test_webbook_data(self):
        # Pick some random points of NIST webbook to calculate Prandlt numbers, instead of re-using tables of different functions
        # Source: https://webbook.nist.gov/cgi/fluid.cgi?Action=Load&ID=C7732185&Type=IsoBar&Digits=5&P=0.1&THigh=1000&TLow=300&TInc=100&RefState=DEF&TUnit=K&PUnit=MPa&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm

        fp = thermo.prop.FluidProperties("HEOS::Water")

        T = 300 # [K]
        p = 1e5 # [Pa]

        exp_Pr = 5.848606793157688
        res_Pr = fp.get_Prandtl(T=T, p=p)
        self.assertAlmostEqual(exp_Pr, res_Pr, places=1)

        T = 1000 # [K]
        p = 1e5 # [Pa]

        exp_Pr = 0.8875173631353968
        res_Pr = fp.get_Prandtl(T=T, p=p)
        self.assertAlmostEqual(exp_Pr, res_Pr, delta=0.02*exp_Pr)


        # Source for 1000K : https://webbook.nist.gov/cgi/fluid.cgi?Action=Load&ID=C7732185&Type=IsoBar&Digits=5&P=1&THigh=1000&TLow=300&TInc=100&RefState=DEF&TUnit=K&PUnit=MPa&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm
        T = 300 # [K]
        p = 10e5 # [Pa]

        exp_Pr = 5.840090755325594
        res_Pr = fp.get_Prandtl(T=T, p=p)
        self.assertAlmostEqual(exp_Pr, res_Pr, places=1)

        T = 1000 # [K]
        p = 10e5 # [Pa]

        exp_Pr = 0.88886909
        res_Pr = fp.get_Prandtl(T=T, p=p)
        self.assertAlmostEqual(exp_Pr, res_Pr, places=1)

class TestReynolds(unittest.TestCase):
    def test_webbook_data(self):
        # Pick some random points of NIST webbook to calculate Reynolds numbers, instead of re-using tables of different functions
        # Source: https://webbook.nist.gov/cgi/fluid.cgi?Action=Load&ID=C7732185&Type=IsoBar&Digits=5&P=0.1&THigh=1000&TLow=300&TInc=100&RefState=DEF&TUnit=K&PUnit=MPa&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm

        fp = thermo.prop.FluidProperties("water")

        T = 300 # [K]
        p = 1e5 # [Pa]
        u = 1 # [m/s]
        L = 1 # [m]

        exp_Re = 1167164.423831442
        res_Re = fp.get_Reynolds_from_velocity(T=T,p=p,u=u,L_ref=L)
        self.assertAlmostEqual(exp_Re,res_Re,places=-3)

        u = u*2 # [m/s]
        exp_Re = 2*1167164.423831442
        res_Re = fp.get_Reynolds_from_velocity(T=T,p=p,u=u,L_ref=L)
        self.assertAlmostEqual(exp_Re,res_Re,places=-3)

        L = 1/3 # [m]
        exp_Re = 2/3*1167164.42383144
        res_Re = fp.get_Reynolds_from_velocity(T=T,p=p,u=u,L_ref=L)
        self.assertAlmostEqual(exp_Re,res_Re,places=-3)

        # Source: https://webbook.nist.gov/cgi/fluid.cgi?Action=Load&ID=C7732185&Type=IsoBar&Digits=5&P=1&THigh=1000&TLow=300&TInc=100&RefState=DEF&TUnit=K&PUnit=MPa&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm

        T = 1000 # [K]
        p = 10e5 # [Pa]
        u = 1e-3 # [m/s]
        L = 10e-3 # [m]

        exp_Re = 0.5772787669412702
        res_Re = fp.get_Reynolds_from_velocity(T=T,p=p,u=u,L_ref=L)
        self.assertAlmostEqual(exp_Re,res_Re,places=3)


class TestSpecificGasConstant(unittest.TestCase):
    def test_gas_constant(self):
        fp = thermo.prop.FluidProperties("HEOS::Water")
        
        specific_gas_constant = 461.52280831345604 # Calculated from wikipedia values
        res = fp.get_specific_gas_constant()
        self.assertAlmostEqual(specific_gas_constant,res,2)

