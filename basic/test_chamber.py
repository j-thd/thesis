import unittest
import basic.chamber as chamber

class TestIdealEnthalpyChange(unittest.TestCase):
    def test_from_table(self):
        # Test a few cases by calculating enthalpy differences from table
        # T1, p1, T2, p2, delta_h, places
        in_out = (
        (300, 0.5e6, 1000, 0.5e6, 3875.448e3, -1),
        (300, 0.5e6, 400, 0.5e6, 420.105e3, 0),
        (330, 0.25e6, 410, 0.25e6, 2498.946e3,-1)
        )

        for T_1, p_1, T_2, p_2, delta_h, places in in_out:
            res = chamber.ideal_enthalpy_change(T_inlet=T_1, p_inlet=p_1, T_outlet=T_2, p_outlet=p_2)
            self.assertAlmostEqual(res, delta_h, places)

class TestIdealPowerConsumption(unittest.TestCase):
    def test_from_table(self):
        # Same as previous test, but with mass flow simply added
        in_out = (
        (0, 300, 0.5e6, 1000, 0.5e6, 3875.448e3*0, -1),
        (1, 300, 0.5e6, 400, 0.5e6, 420.105e3, 0),
        (2, 330, 0.25e6, 410, 0.25e6, 2498.946e3*2,-2)
        )

        for m, T_1, p_1, T_2, p_2, delta_h, places in in_out:
            res = chamber.ideal_power_consumption(mass_flow=m, T_inlet=T_1, p_inlet=p_1, T_outlet=T_2, p_outlet=p_2)
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

