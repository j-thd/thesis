import unittest

import thermo.two_phase as tp

class TestDryoutQuality(unittest.TestCase):
    def test_simple_input(self):
        """Test with pressure of 98 bar an mass flux of G=1000 for simple result. Then vary on that
        Eq. 12.92 and 12.77 from Carey2008
        """
        p = 98e5 # [Pa]
        m_dot = 1000 # [kg/s]
        A_channel = 1 # [m^2]

        D_hydr = 8e-3 # [m] 

        exp_x = 0.60 # [-]
        res_x = tp.dryout_quality(p=p, m_dot=m_dot, A_channel=A_channel, D_hydr=D_hydr)
        self.assertAlmostEqual(exp_x, res_x, delta=1e-10*res_x)

        p = 49e5 # [Pa] Same, but halved pressure
        exp_x = 0.75 # [-]
        res_x = tp.dryout_quality(p=p, m_dot=m_dot, A_channel=A_channel, D_hydr=D_hydr)
        self.assertAlmostEqual(exp_x, res_x, delta=1e-10*res_x)

        m_dot = 500 # [kg/s] G = 500 kg/m^2
        exp_x = 1.0606601717798214
        res_x = tp.dryout_quality(p=p, m_dot=m_dot, A_channel=A_channel, D_hydr=D_hydr)
        self.assertAlmostEqual(exp_x, res_x, delta=1e-10*res_x)

        A_channel = 0.5
        exp_x = 0.75 # [-]
        res_x = tp.dryout_quality(p=p, m_dot=m_dot, A_channel=A_channel, D_hydr=D_hydr)
        self.assertAlmostEqual(exp_x, res_x, delta=1e-10*res_x)

        # try with 5 bar (import to know anyway)
        p=5e5
        exp_x = 0.46488206444593666
        res_x = tp.dryout_quality(p=p, m_dot=m_dot, A_channel=A_channel, D_hydr=D_hydr)
        self.assertAlmostEqual(exp_x, res_x, delta=1e-10*res_x)

        # Now decrease channel diameter to check if 12.77 adaption works properly
        D_hydr = 4e-3 # [m]
        exp_x = 0.5158189468210879
        res_x = tp.dryout_quality(p=p, m_dot=m_dot, A_channel=A_channel, D_hydr=D_hydr)
        self.assertAlmostEqual(exp_x, res_x, delta=1e-10*res_x)

        # Try 1 mm
        D_hydr = 1e-3 # [m]
        exp_x = 0.6350476146762407
        res_x = tp.dryout_quality(p=p, m_dot=m_dot, A_channel=A_channel, D_hydr=D_hydr)
        self.assertAlmostEqual(exp_x, res_x, delta=1e-10*res_x)
        
        # Try 100 micrometer (more realistic)
        D_hydr = 100e-6 # [m]
        exp_x = 0.897028598353314
        res_x = tp.dryout_quality(p=p, m_dot=m_dot, A_channel=A_channel, D_hydr=D_hydr)
        self.assertAlmostEqual(exp_x, res_x, delta=1e-10*res_x)

        
        m_dot = 100 # [kg/s] Now with more realistic mass flux, too ( G= 100kg/m^2 )
        A_channel = 1 # [m^2]
        exp_x = 2.8366534971048387
        res_x = tp.dryout_quality(p=p, m_dot=m_dot, A_channel=A_channel, D_hydr=D_hydr)
        self.assertAlmostEqual(exp_x, res_x, delta=1e-10*res_x)
        