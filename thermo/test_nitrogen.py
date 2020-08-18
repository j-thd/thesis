import unittest

from thermo.prop import FluidProperties

class TestSpecificGasConstant(unittest.TestCase):

    def test_gas_constant(self):
        fp = FluidProperties("nitrogen")

        exp = 296.80305204485137
        res = fp.get_specific_gas_constant()
        self.assertAlmostEqual(exp,res, places=2)

class TestEnthalpy(unittest.TestCase):
    def test_table(self):
        # Page 1410 of Span2000 
        fp = FluidProperties("nitrogen")

        pressure = 1e5 # 1 bar of pressure
        molar_mass = 0.02801348 # [kg/mol] also from Span2000

        # T in, enthalpy out
        # Enthalpy taken in J/mol from table
        in_out = (  (300, 8717.7, 0),
                    (330, 9593.0, -1),
                    (450, 13105, -1),
                    (600, 17564, -1),
                    (1000, 30135, -2))

        for T, h, places in in_out:
            res_h = fp.get_enthalpy(T=T,p=pressure) 
            self.assertAlmostEqual(res_h, h/molar_mass, places=places) # Must be divided by molar mass as CoolProp gives enthalpy for J/kg

        # Testing again, but for another pressure

        pressure = 1e6
        in_out = (  (300, 8662.5, -1),
                    (320, 9253.6, -1),
                    (400, 11612, -1),
                    (550, 16060, -2),
                    (800, 23728, -2),
                    (1000, 30152, -2))

        for T, h, places in in_out:
            res_h = fp.get_enthalpy(T=T,p=pressure) 
            self.assertAlmostEqual(res_h, h/molar_mass, places=places) # Must be divided by molar mass as CoolProp gives enthalpy for J/kg

class TestCp(unittest.TestCase):
    def test_table(self):
        # Page 1410 of Span2000 
        fp = FluidProperties("nitrogen")

        pressure = 1e5 # 1 bar of pressureM
        # T in, cpout
        in_out = (  (300, 1041.2844102196514, 0),
                    (330, 1041.6413812207552, 0),
                    (450, 1049.4947432450376, 0),
                    (600, 1075.1966553245081, 0),
                    (1000, 1167.2951736092768, 0))

        for T, cp, places in in_out:
            res_cp = fp.get_cp(T=T,p=pressure) 
            self.assertAlmostEqual(res_cp, cp, places=places)

        # Testing again, but for another pressure

        pressure = 1e6
        in_out = (  (300, 1055.9202212649052, 0),
                    (320, 1054.1353662593865, 0),
                    (400, 1052.3505112538678, 0),
                    (550, 1068.4142063035367, 0),
                    (800, 1123.7447114746187, 0),
                    (1000, 1168.3660866125879, 0))

        for T, cp, places in in_out:
            res_cp = fp.get_cp(T=T,p=pressure) 
            self.assertAlmostEqual(res_cp, cp, places=places) # Must be divided by molar mass as CoolProp gives enthalpy for J/kg

class TestDensity(unittest.TestCase):
    def test_table(self):
        # Test if the table of the original paper can be reproduced Span2000
        # Page 1410 of Span2000 
        fp = FluidProperties("nitrogen")

        pressure = 1e5 # 1 bar of pressure

        # T in, density out
        # Density taken in kg/mol from table, divided by molar mass
        in_out = (  (300, 1.12328452104, 4),
                    (330, 1.0209512786000001, 4),
                    (450, 0.74846415864, 4),
                    (600, 0.5613060987599999, 4),
                    (1000, 0.33680607004, 4))

        for T, rho, places in in_out:
            res_rho = fp.get_density(T=T,p=pressure) 
            self.assertAlmostEqual(res_rho, rho, places=places)

        # Same again but for different pressure
        pressure = 1e6 # 10 bars of pressure

        # T in, density out
        # Density taken in kg/mol from table, divided by molar mass
        in_out = (  (300, 11.248812894, 4),
                    (330, 10.2058710336, 3),
                    (450, 7.459989724000001, 3),
                    (600, 5.591490608, 3),
                    (1000, 3.3576957128, 3))

        for T, rho, places in in_out:
            res_rho = fp.get_density(T=T,p=pressure) 
            self.assertAlmostEqual(res_rho, rho, places=places)

class TestPressure(unittest.TestCase):
    def test_table(self):
        # Page 1410 of Span2000 and further. Values must be converted to kg/m^3 from mol/dm^3 with molar mass = 0.02801348e3
        fp = FluidProperties("nitrogen")
        # Just use a single temperature
        T = 300 # [K]

        # density in, pressure out
        in_out = (  (1.1232845210400002, 1e5, -1),
                    (2.2469612308, 2e5, -1),
                    (5.6200643576000004, 5e5, -2),
                    (11.248812894, 10e5, -1),
                    (16.883724396, 15e5, -2),
                    (22.523398189599998, 20e5,-1))

        for rho, p, places in in_out:
            res_p = fp.get_pressure(T=T,rho=rho)
            self.assertAlmostEqual(res_p, p, places=places)

        # Again, but different temperature
        T=1000 # [K]

        in_out = (  (0.33680607004, 1e5, 0),
                    (0.67338803224, 2e5, -1),
                    (1.68170523136, 5e5, -1,),
                    (3.35775173976, 10e5, -2),
                    (5.027859390400001, 15e5, -2),
                    (6.6921402372, 20e5, -2))
        
        for rho, p, places in in_out:
            res_p = fp.get_pressure(T=T,rho=rho)
            self.assertAlmostEqual(res_p, p, places=places)

class TestSpecificHeatRatio(unittest.TestCase):
    def test_table(self):
        # Test reported values of gamma by taking cp/cv from Span2000 (page 1410 and below)
        # Test if the table of the original paper can be reproduced Span2000
        # Page 1410 of Span2000 
        fp = FluidProperties("nitrogen")

        pressure = 1e5 # 1 bar of pressure

        # T in, gamma out
        # gamma not explicitly in table, but cp/cv is put in here instead
        in_out = (  (300, 1.402784445511282, 2),
                    (330,  1.4021113243761996, 2),
                    (450, 1.3956356736242885, 3),
                    (600, 1.3821100917431193, 3),
                    (1000, 1.3411234112341124, 3))

        for T, gamma, places in in_out:
            res_gamma = fp.get_specific_heat_ratio(T=T,p=pressure) 
            self.assertAlmostEqual(res_gamma, gamma, places=places)

        # Again but with different pressure
        pressure = 1e6 # 10 bar of pressure

        # T in, gamma out
        # gamma not explicitly in table, but cp/cv is put in here instead
        in_out = (  (300, 1.41666666666666672, 2),
                    (330,  1.4126376256582096, 2),
                    (450, 1.400947867298578, 3),
                    (600, 1.3846859238881248, 3),
                    (1000, 1.3419434194341942, 3))

        for T, gamma, places in in_out:
            res_gamma = fp.get_specific_heat_ratio(T=T,p=pressure) 
            self.assertAlmostEqual(res_gamma, gamma, places=places)

class TestViscosity(unittest.TestCase):
    def test_table(self):
        # Table V from Lemmon2004 has density as inputs, which is not interesting for our purposes, so first the pressure must be obtained

        fp = FluidProperties("nitrogen")

        # T in, Density in mol/dm^3, viscosity out Pa*s
        in_out = (  (100, 25.0, 79.7418e-6,10),
                    (200, 10.0, 21.0810e-6,10),
                    (300, 5.0, 20.7430e-6,10))

        molar_mass = 0.02801348 # [kg/mol] also from Span2000

        for T , rho, mu, places in in_out:
            p = fp.get_pressure(T=T, rho=rho*molar_mass*1e3)
            res_mu = fp.get_viscosity(T=T, p=p)
            self.assertAlmostEqual(res_mu, mu, places=places)

class TestThermalConductivity(unittest.TestCase):
    def test_table(self):
        # Table V from Lemmon2004 has density as inputs, which is not interesting for our purposes, so first the pressure must be obtained

        fp = FluidProperties("nitrogen")

        # T in, Density in mol/dm^3, thermal conductivity in W/(m*K)
        in_out = (  (100, 25.0, 103.834e-3,6),
                    (200, 10.0, 36.0099e-3,7),
                    (300, 5.0, 32.7694e-3,7))

        molar_mass = 0.02801348 # [kg/mol] also from Span2000

        for T , rho, kappa, places in in_out:
            p = fp.get_pressure(T=T, rho=rho*molar_mass*1e3)
            res_kappa = fp.get_thermal_conductivity(T=T, p=p)
            self.assertAlmostEqual(res_kappa, kappa, places=places)

class TestPrandtl(unittest.TestCase):
    def test_webbook_data(self):
        # Pick some random points of NIST webbook to calculate Prandlt numbers, instead of re-using tables of different functions
        # Source: https://webbook.nist.gov/cgi/fluid.cgi?T=300&PLow=1&PHigh=10&PInc=1&Applet=on&Digits=5&ID=C7727379&Action=Load&Type=IsoTherm&TUnit=K&PUnit=bar&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm&RefState=DEF

        fp = FluidProperties("nitrogen")

        T = 300 # [K]
        p = 2e5 # [Pa]

        exp_Pr = 0.721667851210939
        res_Pr = fp.get_Prandtl(T=T, p=p)
        self.assertAlmostEqual(exp_Pr, res_Pr, places=2)

        p = 10e5
        exp_Pr= 0.7279411479231152
        res_Pr = fp.get_Prandtl(T=T, p=p)
        self.assertAlmostEqual(exp_Pr, res_Pr, places=2)

        # Source for 1000K : https://webbook.nist.gov/cgi/fluid.cgi?Action=Load&ID=C7727379&Type=IsoTherm&Digits=5&PLow=1&PHigh=10&PInc=1&T=1000&RefState=DEF&TUnit=K&PUnit=bar&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm
        T = 1000 # [K]
        p = 2e5 # [Pa]

        exp_Pr = 0.7359141666666666
        res_Pr = fp.get_Prandtl(T=T, p=p)
        self.assertAlmostEqual(exp_Pr, res_Pr, places=1)

        T = 1000 # [K]
        p = 10e5 # [Pa]

        exp_Pr = 0.7361587678944341
        res_Pr = fp.get_Prandtl(T=T, p=p)
        self.assertAlmostEqual(exp_Pr, res_Pr, places=1)

class TestReynolds(unittest.TestCase):
    def test_webbook_data(self):
        # Pick some random points of NIST webbook to calculate Reynolds numbers, instead of re-using tables of different functions
        # Source: https://webbook.nist.gov/cgi/fluid.cgi?T=300&PLow=1&PHigh=10&PInc=1&Applet=on&Digits=5&ID=C7727379&Action=Load&Type=IsoTherm&TUnit=K&PUnit=bar&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm&RefState=DEF

        fp = FluidProperties("nitrogen")

        T = 300 # [K]
        p = 1e5 # [Pa]
        u = 1 # [m/s]
        L = 1 # [m]

        exp_Re = 62764.70916913448
        res_Re = fp.get_Reynolds_from_velocity(T=T,p=p,u=u,L_ref=L)
        self.assertAlmostEqual(exp_Re,res_Re,places=-2)

        u = u*2 # [m/s]
        exp_Re = 2*62764.70916913448
        res_Re = fp.get_Reynolds_from_velocity(T=T,p=p,u=u,L_ref=L)
        self.assertAlmostEqual(exp_Re,res_Re,places=-2)

        L = 1/3 # [m]
        exp_Re = 2/3*62764.70916913448
        res_Re = fp.get_Reynolds_from_velocity(T=T,p=p,u=u,L_ref=L)
        self.assertAlmostEqual(exp_Re,res_Re,places=-2)

        T = 1000 # [K]
        p = 10e5 # [Pa]
        u = 1e-3 # [m/s]
        L = 10e-3 # [m]

        exp_Re = 0.806359422656644
        res_Re = fp.get_Reynolds_from_velocity(T=T,p=p,u=u,L_ref=L)
        self.assertAlmostEqual(exp_Re,res_Re,places=2)
