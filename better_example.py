from dscribe.descriptors.turbosoap import TurboSOAP

turbosoap = TurboSOAP(species = {'H': None,  'O': None}, lmax=8)

from ase.build import molecule

#Water molecule
water = molecule("H2O")

#Create TurboSOAP output for water
turbosoap_water = turbosoap.create(water)

print(turbosoap_water.shape)
print(turbosoap_water)

from dscribe.descriptors.turbosoap import TurboSOAPSpecie

# More control on each species
H_config = TurboSOAPSpecie(rcut=5.0, nmax=8, buffer=0.5, atom_sigma_r = 0.3, atom_sigma_t = 0.5,
        atom_sigma_r_scaling = 0.05, atom_sigma_t_scaling = 0.025)
O_config = TurboSOAPSpecie(rcut=5.0, nmax=8, buffer=0.5, atom_sigma_r = 0.2, atom_sigma_t = 0.3,
        atom_sigma_r_scaling = 0.03, atom_sigma_t_scaling = 0.020)

turbosoap = TurboSOAP(species = {'H': H_config,  'O': O_config}, lmax=10)

from dscribe.descriptors.turbosoap import TurboSOAP
from ase.build import molecule
import numpy as np

turbosoap = TurboSOAP(species = {'H': None,  'O': None}, lmax=8)
water = molecule("H2O")
turbosoap_water = turbosoap.create(water)
#Compare similarity of hydrogen and oxygen environments
#Index 0 has oxygen atom and indices 1,2 have hydrogen atoms.
print(water.get_chemical_symbols())

#Compare hydrogen environments. Should be close to one 
print(np.dot(turbosoap_water[1], turbosoap_water[2]))

#Compare oxygen-hydrogen. 
#Should be quite small as the hydrogen and oxygen environments are quite different
print(np.dot(turbosoap_water[0], turbosoap_water[1]))
