import unittest

from dscribe.descriptors import TurboSOAP
from dscribe.descriptors.turbosoap import TurboSOAPSpecie
import numpy as np
from ase.build import molecule, bulk

class TestTurboSoapInterface(unittest.TestCase):
    def testWaterMoleculeSymbols(self):
        turbosoap = TurboSOAP(species = {'H': None,  'O': None}, lmax=8)
        water = molecule("H2O")
        turbosoap_water = turbosoap.create(water)
        self.assertEqual(turbosoap_water.shape[0], 3)
        self.assertEqual(turbosoap_water.shape[1], 1224)

    def testWaterMoleculeAtomicNumbers(self):
        turbosoap = TurboSOAP(species = {1: None,  8: None}, lmax=8)
        water = molecule("H2O")
        turbosoap_water = turbosoap.create(water)
        self.assertEqual(turbosoap_water.shape[0], 3)
        self.assertEqual(turbosoap_water.shape[1], 1224)

    def testInvalidSymbols(self):
        with self.assertRaises(ValueError):
            turbosoap = TurboSOAP(species = {'G': None,  'O': None}, lmax=8)

    def testSpeciesMissing(self):
        turbosoap = TurboSOAP(species = {'H': None}, lmax=8)
        water = molecule("H2O")
        with self.assertRaises(ValueError):
            turbosoap_water = turbosoap.create(water)

    def testMultipleDefs(self):
        with self.assertRaises(ValueError):
            turbosoap = TurboSOAP(species = {'H': None,  'O': None, 1: None}, lmax=8)
    
    def testMultipleSpeciesDefs(self):
        with self.assertRaises(ValueError):
            turbosoap = TurboSOAP(species = {'G': None,  'O': None}, lmax=8)
        with self.assertRaises(ValueError):
            turbosoap = TurboSOAP(species = {'G': None,  'O': None}, lmax=8)

    def testRcut(self):
        config = TurboSOAPSpecie(rcut=0.9)
        with self.assertRaises(ValueError):
            turbosoap = TurboSOAP(species = {'H': config,  'O': None}, lmax=8)




class TestTurboSoapDescriptor(unittest.TestCase):
    def testWaterMolecule(self):
        H_config = TurboSOAPSpecie(rcut=5.0, nmax=8, buffer=0.5, atom_sigma_r = 0.3, atom_sigma_t = 0.5,
            atom_sigma_r_scaling = 0.05, atom_sigma_t_scaling = 0.025)
        O_config = TurboSOAPSpecie(rcut=5.0, nmax=8, buffer=0.5, atom_sigma_r = 0.2, atom_sigma_t = 0.3,
            atom_sigma_r_scaling = 0.03, atom_sigma_t_scaling = 0.020)
        turbosoap = TurboSOAP(species = {'H': H_config,  'O': O_config}, lmax=10)
        water = molecule("H2O")
        turbosoap_water = turbosoap.create(water)
        self.assertAlmostEqual(np.dot(turbosoap_water[1], turbosoap_water[2]), 1.0, delta=1e-7)
        self.assertAlmostEqual(np.dot(turbosoap_water[0], turbosoap_water[1]), 0.0, delta=1e-3)
    
    def testCuFcc(self):
        cu_fcc = bulk('Cu', 'fcc', a=3.6, cubic=True)
        config = TurboSOAPSpecie(rcut=2.5)
        turbosoap = TurboSOAP(species = {'Cu': config}, lmax=8)
        turbosoap_fcc = turbosoap.create(cu_fcc, periodic=True)
        self.assertEqual(turbosoap_fcc.shape[0], len(cu_fcc))
        self.assertEqual(turbosoap_fcc.shape[1], 324)
        self.assertAlmostEqual(np.dot(turbosoap_fcc[0], turbosoap_fcc[1]), 1.0, delta=1e-7)

    #More periodic tests to make sure everything is working properly


if __name__ == "__main__":
    unittest.main()