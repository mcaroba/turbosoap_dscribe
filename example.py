#Simple example with default parameters
from dscribe.descriptors.turbosoap import TurboSOAP
from ase.build import molecule

turbosoap = TurboSOAP(species = {'H': None,  'O': None}, lmax=8)

#Water molecule
water = molecule("H2O")

#Create TurboSOAP output for water
turbosoap_water = turbosoap.create(water)

print(turbosoap_water.shape)
print(turbosoap_water)

#More involved example with explicit parameters
from turbosoap_dscribe import TurboSOAPSpecie

# More control on each species
H_config = TurboSOAPSpecie(rcut=5.0, nmax=8, buffer=0.5, atom_sigma_r = 0.3, atom_sigma_t = 0.5,
        atom_sigma_r_scaling = 0.05, atom_sigma_t_scaling = 0.025)
O_config = TurboSOAPSpecie(rcut=5.0, nmax=8, buffer=0.5, atom_sigma_r = 0.2, atom_sigma_t = 0.3,
        atom_sigma_r_scaling = 0.03, atom_sigma_t_scaling = 0.020)

turbosoap = TurboSOAP(species = {'H': H_config,  'O': O_config}, lmax=10)

turbosoap_water = turbosoap.create(water)

print(turbosoap_water.shape)
print(turbosoap_water)





