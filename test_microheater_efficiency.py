import unittest

import numpy as np

import microheater_efficiency as me
import physical_constants

class TestCalcP_Delta_h(unittest.TestCase):
    def test_simple_input(self):
        T_in = 300 # [K]
        T_out = 400 # [K]
        m_dot = 1e-6 # [kg/s}]
        p_chamber = 1e5 # [Pa]
        # Expected output with enthalpy change from webbook at 1 bar
        exp = 2.61775 # [W]
        res = me.calc_P_delta_h(T_in=T_in,T_out=T_out,m_dot=m_dot,p_chamber=p_chamber) # [W]
        self.assertAlmostEqual(exp,res,places=4)

        # Check agian with different mass flow
        m_dot = 1e-3 # [kg/s]
        exp = 2.61775e3 # [W]
        res = me.calc_P_delta_h(T_in=T_in,T_out=T_out,m_dot=m_dot,p_chamber=p_chamber) # [W]
        self.assertAlmostEqual(exp, res, places = 1)

        # Different outlet temperature
        T_out = 500 # [K]
        res = me.calc_P_delta_h(T_in=T_in,T_out=T_out,m_dot=m_dot,p_chamber=p_chamber) # [W]
        exp = 2815.95 # [W]
        self.assertAlmostEqual(exp, res, places = 1)

        # Different pressure 
        p_chamber = 5e5 # [Pa]
        exp = 2799.68 # [W]
        res = me.calc_P_delta_h(T_in=T_in,T_out=T_out,m_dot=m_dot,p_chamber=p_chamber) # [W]
        self.assertAlmostEqual(exp, res, places=1)

class TestPRad(unittest.TestCase):
    def test_simple_input(self):
        # 400 K temperature
        emissivity = 1
        A_mh = 1
        T_mh_range = np.array([400], dtype=object)
        
        res = me.calc_P_rad(T_mh_range=T_mh_range, emissivity=emissivity,A_mh=A_mh)
        exp = 25600000000*physical_constants.stefan_boltzmann
        self.assertEqual(res, exp)

        # Same but temperature 1
        emissivity = 1
        A_mh = 1
        T_mh_range = np.array([1], dtype=object)
        
        res = me.calc_P_rad(T_mh_range=T_mh_range, emissivity=emissivity,A_mh=A_mh)
        exp = physical_constants.stefan_boltzmann
        self.assertEqual(res, exp)

        # same but with half emmisivity
        emissivity = 0.5
        res = me.calc_P_rad(T_mh_range=T_mh_range, emissivity=emissivity,A_mh=A_mh)
        exp = 0.5*physical_constants.stefan_boltzmann
        self.assertEqual(res, exp)
        # Same, but now with larger array, second temperature is T=2
        T_mh_range = np.array([1,2],dtype=object)
        exp = np.array([0.5, 8])*physical_constants.stefan_boltzmann
        res = me.calc_P_rad(T_mh_range=T_mh_range, emissivity=emissivity,A_mh=A_mh)
        # np.equal() works on arrays with different dtype, and does element-wise comparison, and then alltrue can check if all elements are the same
        self.assertTrue(np.alltrue(np.equal(res,exp)))

        # One more time but with different area
        A_mh = 100
        exp = np.array([50, 800])*physical_constants.stefan_boltzmann
        res = me.calc_P_rad(T_mh_range=T_mh_range, emissivity=emissivity,A_mh=A_mh)
        # np.equal() works on arrays with different dtype, and does element-wise comparison, and then alltrue can check if all elements are the same
        self.assertTrue(np.alltrue(np.equal(res,exp)))

        # Test to see if it can handle high values of T without overflow errors
        T_mh_range = np.array([1,10,100,1000,10000],dtype=object)
        exp = np.array([50, 50*1e4, 50*1e8, 50*1e12, 50*1e16])*physical_constants.stefan_boltzmann
        res = me.calc_P_rad(T_mh_range=T_mh_range, emissivity=emissivity,A_mh=A_mh)
        
        self.assertTrue(np.allclose(exp,res.astype('float64')))

        # Test if entering a non-object numpy array throws an error
        T_mh_range = np.array([1,10,100])
        with self.assertRaises(TypeError):
            res = me.calc_P_rad(T_mh_range=T_mh_range, emissivity=emissivity,A_mh=A_mh)

        # Now test with a regular python list
        T_mh_range = [1,10,100]
        with self.assertRaises(TypeError):
            res = me.calc_P_rad(T_mh_range=T_mh_range, emissivity=emissivity,A_mh=A_mh)

class TestCalcMicroheaterEfficiency(unittest.TestCase):
    # Only added functionality on top of helper functions is the temperature margin.
    # So tests can be repeated and, and then margin can be varied
    def test_single_input(self):
        # Previous tests combined into efficiency
        T_in = 300 # [K]
        m_dot = 1e-6 # [kg/s}]
        p_chamber = 1e5 # [Pa]
        emissivity = 1 #[-]
        A_mh = 1 # [m^2]
        T_mh_range = np.array([400], dtype=object)
        T_margin = 0 # [K]
        
        exp_P_delta_h = 2.61775
        exp_P_rad = 25600000000*physical_constants.stefan_boltzmann

        # Efficiency
        exp = exp_P_delta_h / (exp_P_delta_h+exp_P_rad) 
        res = me.calc_microheater_efficiency(T_mh_range=T_mh_range,emissivity=emissivity,A_mh=A_mh,T_in=T_in,T_margin=T_margin,m_dot=m_dot,p_chamber=p_chamber)
        self.assertAlmostEqual(res[0],exp,places=4)

        T_mh_range = np.array([400, 500], dtype=object)
        m_dot = 1e-3 
    
        res = me.calc_microheater_efficiency(T_mh_range=T_mh_range,emissivity=emissivity,A_mh=A_mh,T_in=T_in,T_margin=T_margin,m_dot=m_dot,p_chamber=p_chamber)
        exp_P_delta_h = np.array([2617.75, 2815.95])
        exp_P_rad = np.array([25600000000, 62500000000])*physical_constants.stefan_boltzmann
        exp = exp_P_delta_h / (exp_P_delta_h + exp_P_rad)
        self.assertTrue(np.allclose(res.astype('float64'),exp))
        # Same as previous function tests some mass_range, however, two inputs now

        # Test with T_margin non-zero
        T_mh_range = np.array([500, 600], dtype=object)
        T_margin = 100 #  [K]
        exp_P_delta_h = np.array([2617.75, 2815.95])
        exp_P_rad = np.array([62500000000,129600000000])*physical_constants.stefan_boltzmann
        exp = exp_P_delta_h / (exp_P_delta_h + exp_P_rad)
        res = me.calc_microheater_efficiency(T_mh_range=T_mh_range,emissivity=emissivity,A_mh=A_mh,T_in=T_in,T_margin=T_margin,m_dot=m_dot,p_chamber=p_chamber)
        self.assertTrue(np.allclose(res.astype('float64'),exp,rtol=1e-4))

class TestPCond(unittest.TestCase):
    def test_simple_input(self):
        # Test proportionality with simple inputs (should be sufficient)
        k = 1 # [W/mk] thermal conductivity
        T2 = 2 # [K] Temperature on one side
        T1 = 1 # [K] Temperature on other side
        L = 1 # [m] Length or thickness
        A = 1 # [m^2] Area through which heat is conducted

        exp = -1 # Should be minus one, due to sign convention of flow direction. (I.e.: losing heat)
        res = me.calc_P_cond(A_chip=A,T1=T1, T2=T2, L=L, thermal_conductivity=k)
        self.assertEqual(exp,res)

        L = 2 # Should halve heat flow
        exp = -0.5
        res = me.calc_P_cond(A_chip=A,T1=T1, T2=T2, L=L, thermal_conductivity=k)
        self.assertEqual(exp,res)
        
        A = 2 # Should double heat flow
        exp = -1
        res = me.calc_P_cond(A_chip=A,T1=T1, T2=T2, L=L, thermal_conductivity=k)
        self.assertEqual(exp,res)

        T2 = 3 # Should double flow again
        exp = -2
        res = me.calc_P_cond(A_chip=A,T1=T1, T2=T2, L=L, thermal_conductivity=k)
        self.assertEqual(exp,res)
        
        T1 = 2 # Should halve it again
        exp = -1
        res = me.calc_P_cond(A_chip=A,T1=T1, T2=T2, L=L, thermal_conductivity=k)
        self.assertEqual(exp,res)

        k = 2 # Should double it again
        exp = -2
        res = me.calc_P_cond(A_chip=A,T1=T1, T2=T2, L=L, thermal_conductivity=k)
        self.assertEqual(exp,res)

