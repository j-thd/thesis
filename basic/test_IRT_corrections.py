import unittest
import math
import basic.IRT_corrections as IRT_c

class testRajeevModelIntegratin(unittest.TestCase):
    """Use Barry's sheet to see if a full calculation also yields the same results.

        The integration also uses the regular IRT calculation, as the corrections make no sense on their own

    """

    def test_one(self):
        # Try to replicate the first column of calculations (Column D)
        pass

class testDivergenceLossConical2D(unittest.TestCase):
    def test_zero(self):
        # Test if almost zero input gives expected result (singularity present), should converge to 1

        alpha = 1e-8
        res = IRT_c.divergence_loss_conical_2D(alpha)
        self.assertAlmostEqual(res, 1)

    def test_pi(self):
        # Should result in 1/pi
        alpha = 0.5*math.pi 

        res = IRT_c.divergence_loss_conical_2D(alpha)
        self.assertEqual(res, 1/(0.5*math.pi))
        # should give zero
        alpha = math.pi
        res = IRT_c.divergence_loss_conical_2D(alpha)
        self.assertAlmostEqual(res, 0)

    def test_sheet_rajeev(self):
        # Replicate some results from accompanying excel sheet of Barry, with values from Rajeev's model
        alpha = math.radians(20.5)
        exp = 0.9788
        res = IRT_c.divergence_loss_conical_2D(alpha)
        self.assertAlmostEqual(res, exp, places=4)


class testViscousLoss(unittest.TestCase):
    def test_simple(self):
        # Test some simple predicitbale inputs
        area_ratio = 0
        Reynolds_throat_wall = 1

        exp = 17.6
        res = IRT_c.viscous_loss(area_ratio=area_ratio,reynolds_throat_wall=Reynolds_throat_wall)
        self.assertEqual(res, exp)

        Reynolds_throat_wall = 4
        exp = 17.6/2
        res = IRT_c.viscous_loss(area_ratio=area_ratio,reynolds_throat_wall=Reynolds_throat_wall)
        self.assertEqual(res,exp)

        area_ratio = 1/0.0032
        exp = 17.6*math.e / 2
        res = IRT_c.viscous_loss(area_ratio=area_ratio, reynolds_throat_wall=Reynolds_throat_wall)
        self.assertAlmostEqual(res,exp,places=8)

    def test_sheet_barry(self):
        # Replicate some results from Barry's sheet ROW 62 (called thrust coefficient loss)
        # Barry takes into changes of throat reynolds at throat wall (just like Spisz1965)
        area_ratio = 30
        Re_throat = 830
        Re_throat_wall = IRT_c.Reynolds_throat_wall_cold(reynolds_throat=Re_throat)

        exp = 0.7647
        res = IRT_c.viscous_loss(area_ratio=area_ratio, reynolds_throat_wall=Re_throat_wall)
        self.assertAlmostEqual(exp,res, places=4)

        Re_throat = 1277
        Re_throat_wall = IRT_c.Reynolds_throat_wall_cold(reynolds_throat=Re_throat)

        exp = 0.6166
        res = IRT_c.viscous_loss(area_ratio=area_ratio, reynolds_throat_wall=Re_throat_wall)
        self.assertAlmostEqual(exp,res,places=3)
        
    def test_reynolds_throat_wall_hot(self):
        # No source was found for replicating this function, but it is a minor function and very simple.
        # Instead it is compared to a manual calculation after using the equation in the report

        Re_throat = 0
        res = IRT_c.Reynolds_throat_wall_hot(reynolds_throat=Re_throat)
        exp = 0
        self.assertEqual(exp,res)

        Re_throat = 1000
        exp = 1727.0932103916857
        res = IRT_c.Reynolds_throat_wall_hot(reynolds_throat=Re_throat)
        self.assertEqual(exp,res)
        
class testThroatBoundaryLoss(unittest.TestCase):
    # A simple test is lacking because the formula simply isn't very simple at all

    def test_sheet_barry(self):
        # Row 31 to 38 from Barry's sheet
        gamma = 1.403
        Re_throat = 830.193345707
        throat_roc = 1e-6
        throat_longitudinal_radius = 0.5*2.86463e-5

        exp = 0.95153  # This is after removing the extra ^2 from D39 "Help 1"
        res = IRT_c.throat_boundary_loss(gamma=gamma,reynolds_throat=Re_throat,throat_radius=throat_longitudinal_radius,throat_roc=throat_roc)
        self.assertAlmostEqual(exp,res, places=3)

class TestHydraulicDiameter(unittest.TestCase):

    def test_sheet_barry(self):
        # Row 7 to 11 from Barry's sheet
        A = 1.4094000E-09 # [m^2] area
        P = 1.96800E-04 # [m] perimeter

        exp_Dh = 2.86463414634E-05 # [m] Expected hydraulic diameter
        res_Dh = IRT_c.hydraulic_diameter(A=A,wetted_perimeter=P)
        self.assertAlmostEqual(exp_Dh,res_Dh,places=16)


