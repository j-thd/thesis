import unittest
import math
import basic.IRT as IRT

class TestVdk(unittest.TestCase):

    def test_zero_two(self):
        # Test if zero and two give expected results
        in_out = (
            (0, 0),
            (2, 0.769800358919501)
        )

        for input, output in in_out:
            res = IRT.vdk(input)
            self.assertEqual(res, output)

    def test_table_from_TRP(self):
        # Page 51 from the TRP reader has a table with in and outputs. It is used to verify the function here
        in_out = (
            (1.05, 0.6177),
            (1.10, 0.6284),
            (1.12, 0.6325),
            (1.19, 0.6466),
            (1.23, 0.6543),
            (1.31, 0.6691),
            (1.38, 0.6813),
            (1.65, 0.7238)
        )

        for input, output in in_out:
            res = IRT.vdk(input)
            self.assertEqual(round(res,4), output)

class TestMassFlow(unittest.TestCase):
    def test_zero_two(self):
        expected_result = 0.769800358919501
        res = IRT.mass_flow(p_chamber=1,A_throat=1,R=1,T_chamber=1,gamma=2)
        self.assertEqual(res,expected_result)

    def test_Anderson_example(self):
        # Example 10.5 from Anderson2016 (recalculated manually due to different formula, calculation approach)
        expected_result = 586.1 # [kg/s]
        res = IRT.mass_flow(p_chamber=3.03e6, A_throat=0.4, R=520,T_chamber=3500,gamma=1.22)
        self.assertEqual(round(res,1),expected_result)

class TestAreaRatio(unittest.TestCase):
    def test_one(self):
        # Mach number is 1 in the throat so area ratio must be 1.
        expected_result = 1
        res = IRT.area_ratio(M=1,gamma=1.4)
        self.assertEqual(res, expected_result)

    def test_Anderson_table(self):
        # Anderson2016 has a table in the back of his book for area ratio vs. Mach numbers
        in_out = (
            (0.2e-1, 0.2894e2, 2),
            (0.1, 0.5822e1, 3),
            (0.5, 0.134e1, 3),
            (0.78, 0.1047e1, 3),
            (1.22, 0.1037e1, 3),
            (2.45, 0.2517e1, 3),
            (5.2, 0.2928e2,2),
            (7.9, 0.1795e3,1),
            (50, 0.1455e7,-3)
        )
        for input, output, places in in_out:
                res = IRT.area_ratio(M=input, gamma=1.4)
                self.assertAlmostEqual(res, output, places)

class TestTemperatureAndPressureRatio(unittest.TestCase):

    def test_simple_input(self):
        M = 1
        gamma = 2
        expected_temperature_ratio = 1.5
        expected_pressure_ratio = 2.25
        res_temperature_ratio = IRT.temperature_ratio(M=M, gamma=gamma)
        res_pressure_ratio = IRT.pressure_ratio(M=M, gamma=gamma)
        self.assertEqual(expected_temperature_ratio, res_temperature_ratio)
        self.assertEqual(expected_pressure_ratio, res_pressure_ratio)
        
    def test_anderson_table(self):
        # Table in back of Anderson book gives values for gamma=1.4
        gamma = 1.4

        in_out = (
            (0.2e-1, 0.1e1, 3,  0.1e1, 3),
            (0.24e0, 0.1041e1, 3, 0.1012e1, 3),
            (0.56e0, 0.1237e1, 3, 0.1063e1, 3),
            (1, 0.1893e1, 3, 0.12e1, 3),
            (1.64e0, 0.4511e1, 3, 0.1538e1, 3),
            (2.45e0, 0.1581e2, 2, 0.22e1, 3),
            (50, 0.2815e10, -6, 0.5010e3, 1)
        )

        for M, expected_PR, PR_places, expected_TR, TR_places in in_out:
            res_PR = IRT.pressure_ratio(M=M, gamma=gamma)
            res_TR = IRT.temperature_ratio(M=M, gamma=gamma)
            self.assertAlmostEqual(res_PR, expected_PR, PR_places)
            self.assertAlmostEqual(res_TR, expected_TR, TR_places)

class TestIsThroatSonic(unittest.TestCase):
    def testAndersonTable(self): 
        # Pressure ratio slightly increased(by 0.001e6) from Anderson table to work around rounding errors
        res = IRT.is_throat_sonic(p_chamber=1.1291e6, p_back=1e6, exit_area_ratio=1.529, gamma=1.4)
        self.assertTrue(res)
        # Same test but now slightless less to check if it is indeed False
        res = IRT.is_throat_sonic(p_chamber=1.1290e6, p_back=1e6, exit_area_ratio=1.529, gamma=1.4)
        self.assertTrue(not res)

class TestNozzleStatus(unittest.TestCase):
    def testAndersonTable(self):
        # Using table from Anderon to test for different back pressures
        # Taking line for M=0.5 first for subsonic test
        p_chamber = 1.186e6
        exit_area_ratio = 1.340
        gamma = 1.4
        p_exit = 1e6
        res = IRT.nozzle_status(p_chamber=p_chamber,p_back=p_exit, AR_exit=exit_area_ratio,gamma=gamma)
        self.assertEqual(res, 'subsonic')
        p_chamber = 1.54e6 # Slight below value for which it will not have a shock in the nozzle
        res = IRT.nozzle_status(p_chamber=p_chamber,p_back=p_exit, AR_exit=exit_area_ratio,gamma=gamma)
        self.assertEqual(res, 'shock in nozzle')
        # Now test the same but with slightly higher pressure to see if it returns underexpanded
        p_chamber = 1.55e6 # Slight above value for which it will not have a shock in the nozzle
        res = IRT.nozzle_status(p_chamber=p_chamber,p_back=p_exit, AR_exit=exit_area_ratio,gamma=gamma)
        self.assertEqual(res, 'overexpanded')
        # Now isentropic expansion (tolerance is set 0.01)
        # Area ratio must be updated to a value present in Anderson table for M>1
        # M = 2.1 is selected
        exit_area_ratio = 1.837
        p_chamber = 9.145e6
        res = IRT.nozzle_status(p_chamber=p_chamber,p_back=p_exit, AR_exit=exit_area_ratio,gamma=gamma)
        self.assertEqual(res, 'isentropically expanded (margin: 0.01)')
        # Same but with exit pressure slightly below margin set, so it is overexpanded
        p_exit = 0.9899e6
        res = IRT.nozzle_status(p_chamber=p_chamber,p_back=p_exit, AR_exit=exit_area_ratio,gamma=gamma)
        self.assertEqual(res, 'underexpanded')
        # Same but slightly higher so it is underexpanded
        p_exit = 1.01e6
        res = IRT.nozzle_status(p_chamber=p_chamber,p_back=p_exit, AR_exit=exit_area_ratio,gamma=gamma)
        self.assertEqual(res, 'overexpanded')
        # Table values are rounded, so tests are fine if the exit pressure need tiny tweaks to pass.

class TestExitPressureRatioShockwave(unittest.TestCase):
    def testOne(self):
        # Should return pressure ratio for M=1, as shockwave is infinitely weak
        # Value taken from Anderson2006 Appendix A again
        expected = 0.1893e1
        res = IRT.exit_pressure_ratio_shockwave(exit_area_ratio=1, gamma=1.4)
        self.assertAlmostEqual(res, expected,3)
    
    def testAndersonTable(self): # Appendix A and B from Anderson2016
        # For this case the pressure ratios of two tables must be multiplied. Accuracy unknown after multiplication
        # Lets guess 1 digit loss is fine
        gamma = 1.4

        in_out = (
            ( 0.1760e1, 0.8458e1, 0.4736e1, 3), # M = 2.05
            ( 0.1395e2, 0.2247e3, 0.2140e2, 2) # M = 4.3
        )

        for AR, PR_exit, PR_shockwave, places in in_out:
            res = IRT.exit_pressure_ratio_shockwave(exit_area_ratio=AR, gamma=gamma)
            self.assertAlmostEqual(res, PR_exit/PR_shockwave, places)

class TestPressureRatioShockwave(unittest.TestCase):
    def testOne(self):
        # Mach is 1 input should return 1
        expected = 1
        res = IRT.pressure_ratio_shockwave(M=1,gamma=1.4)
        self.assertEqual(res, expected)

    def testAndersonTable(self):
        # Test values from Appendix B from Anderson
        gamma=1.4

        in_out = (
            (0.102e1, 0.1047e1, 3),
            (0.158e1, 0.2746e1, 3),
            (0.345e1, 0.1372e2, 2),
            (0.61e1, 0.4324e2, 2),
            (0.36e2, 0.1512e4, 0)
        )

        for M , PR, places in in_out:
            res = IRT.pressure_ratio_shockwave(M=M,gamma=gamma)
            self.assertAlmostEqual(res, PR, places)

class TestMachFromAreaRatio(unittest.TestCase):
    def testAreaRatioOne(self):
        expected_result = 1.0
        input = 1
        res = IRT.Mach_from_area_ratio(AR=input,gamma=1.4)
        self.assertEqual(res, expected_result)

    def testAndersonTable(self):
        # use Anderson2016 table Appendix for gamma-1.4 in order to verify function workings
        in_out = (
            (0.2894e2, 0.2e-1, False, 5),
            (0.2708e1, 0.22e0 , False, 4),
            (0.1213e1, 0.58e0, False, 4),
            #(0.1000e1, 0.98e0, False, 2), # Function is extremely insensitive to results close to M=1
            #(0.1000e1, 0.102e1, True, 3),
            (0.1126e1, 0.1420e1, True, 3),
            (0.2763e1, 0.2550e1, True, 3),
            (0.3272e3, 0.9e1, True, 3),
            (0.1538e5, 0.2e2, True, 2)

        )
        for AR, expected_M, supersonic, M_places in in_out:
            res_M = IRT.Mach_from_area_ratio(AR=AR,gamma=1.4,supersonic=supersonic)
            self.assertAlmostEqual(res_M, expected_M, places=M_places)


class TestExitVelocity(unittest.TestCase):
    def test_simple_input(self):
        expected = 1.0801234497346432 # Manually calculated
        # area ratio of 1 means the mach number is 1
        res = IRT.exit_velocity(AR_exit=1,T_chamber=1,R=1,gamma=1.4)
        self.assertEqual(res, expected)
    
    def test_TRP_example(self):
        # Slide 57-51 from Ideal Rocket Motor lecture from Thermal Rocket Propulsion course are used as example (Zandbergen)
        expected = 3021
        gamma = 1.3
        R = 390.4
        T_chamber = 3400
        AR_exit = 49
        res = IRT.exit_velocity(AR_exit=AR_exit,T_chamber=T_chamber,R=R, gamma=gamma)
        self.assertAlmostEqual(res, expected, -1) # 4 digits accuracy (-1, rounding difference)

class TestThrust(unittest.TestCase):
    def test_TRP_example(self):
        # Slide 57-51 from Ideal Rocket Motor lecture from Thermal Rocket Propulsion course are used as example (Zandbergen)
        expected = 401.8e3
        gamma = 1.3
        R = 390.4
        p_chamber = 5.6e6
        p_back = 6.03e3
        T_chamber = 3400
        AR_exit = 49
        d_exit = 1.6
        A_exit = 0.25*math.pi*d_exit**2
        A_throat = A_exit/AR_exit
        res = IRT.thrust(p_chamber=p_chamber,T_chamber=T_chamber,A_throat=A_throat, AR_exit=AR_exit,p_back=p_back,gamma=gamma,R=R)
        self.assertAlmostEqual(expected, res, -3) # Slight rounding error again and 1 digit accuracy less

    def test_simple_input(self):
        # Expected exit velocity and mass flow
        # Should be sonic throat too, due to vacuum at exit
        u_exit_expected = 1.0801234497346432
        m_dot_expected = 0.6847314563772704 # Should be equal to Gamma(1.4) (vandankerkhove function)
        p_exit = 0.5282817877171743 # Should be exit pressure for M=1
        AR_exit=1
        T_chamber=1
        R=1
        gamma=1.4
        p_chamber=1
        A_throat = 1
        AR_exit = 1
        p_back = 0

        # Expected thrust
        F_expected = m_dot_expected*u_exit_expected + p_exit*A_throat*AR_exit
        res = IRT.thrust(p_chamber=p_chamber,T_chamber=T_chamber,A_throat=A_throat, AR_exit=AR_exit,p_back=p_back,gamma=gamma,R=R)
        self.assertEqual(res, F_expected)

class TestGetEnginePerformance(unittest.TestCase):
    def test_simple_input(self):
        u_exit_expected = 1.0801234497346432
        m_dot_expected = 0.6847314563772704 # Should be equal to Gamma(1.4) (vandankerkhove function)
        p_exit = 0.5282817877171743 # Should be exit pressure for M=1
        AR_exit=1
        T_chamber=1
        R=1
        gamma=1.4
        p_chamber=1
        A_throat = 1
        AR_exit = 1
        p_back = 0

        # Expected thrust
        F_expected = m_dot_expected*u_exit_expected + p_exit*A_throat*AR_exit
        # Since the throat is M=1 and the back pressure is vacuum, it is clearly underexpanded
        NS = "underexpanded"
        res = IRT.get_engine_performance(p_chamber=p_chamber,T_chamber=T_chamber,A_throat=A_throat, AR_exit=AR_exit, p_back=p_back, gamma=gamma, R=R)
        self.assertEqual(res['thrust'], F_expected)
        self.assertEqual(res['m_dot'], m_dot_expected)
        self.assertEqual(res['u_exit'],u_exit_expected)
        self.assertEqual(res['nozzle_status'], NS)

    def test_TRP_example(self):
        # Slide 57-51 from Ideal Rocket Motor lecture from Thermal Rocket Propulsion course are used as example (Zandbergen)
        F_expected = 401.8e3
        u_exit_expected = 3021
        m_dot_expected = 132.9
        gamma = 1.3
        R = 390.4
        p_chamber = 5.6e6
        p_back = 6.03e3
        T_chamber = 3400
        AR_exit = 49
        d_exit = 1.6
        A_exit = 0.25*math.pi*d_exit**2
        A_throat = A_exit/AR_exit

        NS = "isentropically expanded (margin: 0.01)"
        res = IRT.get_engine_performance(p_chamber=p_chamber,T_chamber=T_chamber,A_throat=A_throat, AR_exit=AR_exit, p_back=p_back, gamma=gamma, R=R)
        self.assertAlmostEqual(res['thrust'], F_expected,-3)
        self.assertAlmostEqual(res['m_dot'], m_dot_expected,0)
        self.assertAlmostEqual(res['u_exit'],u_exit_expected,-1)
        self.assertEqual(res['nozzle_status'], NS)
        
if __name__ == '__main__':
    unittest.main()