TurboSOAP-DScribe is a python package is an optional dependency of the to DScribe python package [https://singroup.github.io/dscribe/](https://singroup.github.io/dscribe/) which provides TurboSOAP machine learning descriptor for DScribe. Requires DScribe and fortran compiler.

# Example
For a quick example please see DScribe documentation of TurboSOAP descriptor [link-here]. 
Some examples below
```python
#Simple example with default parameters
from dscribe.descriptors import TurboSOAP
from ase.build import molecule

turbosoap = TurboSOAP(species = {'H': None,  'O': None}, lmax=8)

#Water molecule
water = molecule("H2O")

#Create TurboSOAP output for water
turbosoap_water = turbosoap.create(water)

print(turbosoap_water.shape)
print(turbosoap_water)

#More involved example with explicit parameters
from dscribe.descriptors import TurboSOAPSpecie

# More control on each species
H_config = TurboSOAPSpecie(rcut=5.0, nmax=8, buffer=0.5, atom_sigma_r = 0.3, atom_sigma_t = 0.5,
        atom_sigma_r_scaling = 0.05, atom_sigma_t_scaling = 0.025)
O_config = TurboSOAPSpecie(rcut=5.0, nmax=8, buffer=0.5, atom_sigma_r = 0.2, atom_sigma_t = 0.3,
        atom_sigma_r_scaling = 0.03, atom_sigma_t_scaling = 0.020)

turbosoap = TurboSOAP(species = {'H': H_config,  'O': O_config}, lmax=10)

turbosoap_water = turbosoap.create(water)

print(turbosoap_water.shape)
print(turbosoap_water)
```

# Installation
The latest stable release is available through pip: (add the --user flag if root access is not available)

```sh
pip install turbosoap-dscribe
```
