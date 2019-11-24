# Unittest for shomate.py
import unittest
import thermo.shomate as shomate


class TestShomate(unittest.TestCase):

    def test_enthalpy_table(self):
        # Tuples with input temperature in K and output in kJ/mol
        # taken from table under shomate equation parameters
        in_out = (
            (298, -0.01),
            (300, 0.14),
            (400, 7.71),
            (500, 15.66))

        for input, output in in_out:
            # Results are given with 2 decimal places, so result must be
            #  rounded too
            res = shomate.getEnthalpy(shomate.water_L, input)
            self.assertEqual(round(res, 2), output)


if __name__ == '__main__':
    unittest.main()
