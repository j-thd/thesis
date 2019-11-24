import unittest

import thermo.water as water

class TestEnthalpy(unittest.TestCase):

    def testTable(self):
        # Test with table from IAPWS 95 report (Wagner2002)

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
            res_h = water.get_enthalpy(T=T,p=pressure)
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
            res_h = water.get_enthalpy(T=T,p=pressure)
            self.assertAlmostEqual(res_h, h, places=places)

class TestSpecificGasConstant(unittest.TestCase):
    def test_gas_constant(self):
        specific_gas_constant = 461.52280831345604 # Calculated from wikipedia values
        res = water.get_specific_gas_constant()
        self.assertAlmostEqual(specific_gas_constant,res,2)

