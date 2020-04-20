TurboSOAP-DScribe is a python package which provides TurboSOAP machine learning descriptor functionality to DScribe python package [https://singroup.github.io/dscribe/](https://singroup.github.io/dscribe/). Requires DScribe and fortran compiler.

# Example
For a quick example please see DScribe documentation of TurboSOAP descriptor [link-here]. 
Some examples below
```python
import numpy as np
from ase.build import molecule
from dscribe.descriptors.turbosoap import TurboSOAP
from turbosoap_dscribe import TurboSOAPSpecie
# Define atomic structures
water = molecule("H2O")

# Setup descriptors
# Use defaults for quick test
H_config = TurboSOAPSpecie(rcut=4.0, nmax=7)
O_config = TurboSOAPSpecie(rcut=5.0, nmax=7)

#Calculate the descriptor
turbosoap_desc = TurboSOAP(lmax=10, species = {'H': H_config,  'O': O_config})
soap = turbosoap_desc.create(water)
print(soap.shape)

# More control on each species
H_config = TurboSOAPSpecie(rcut=5.0, nmax=8, buffer=0.5, atom_sigma_r = 0.3, atom_sigma_t = 0.5,
        atom_sigma_r_scaling = 0.05, atom_sigma_t_scaling = 0.025)
O_config = TurboSOAPSpecie(rcut=5.0, nmax=8, buffer=0.5, atom_sigma_r = 0.2, atom_sigma_t = 0.3,
        atom_sigma_r_scaling = 0.03, atom_sigma_t_scaling = 0.020)
turbosoap_desc = TurboSOAP(lmax=10, species = {'H': H_config,  'O': O_config})

soap = turbosoap_desc.create(water)
print(soap.shape)

#Re-using oxygen configuration with slight modification
H_config = TurboSOAPSpecie(rcut=5.0, nmax=8, buffer=0.7, atom_sigma_r = 0.4, atom_sigma_t = 0.5,
        atom_sigma_r_scaling = 0.1, atom_sigma_t_scaling = 0.025)
import copy
O_config_new = copy.copy(O_config)
O_config_new.nmax = 6
turbosoap_desc = TurboSOAP(lmax=10, species = {'H': H_config,  'O': O_config})

soap = turbosoap_desc.create(water)
print(soap.shape)
```

# Installation
The latest stable release is available through pip: (add the --user flag if root access is not available)

```sh
pip install turbosoap-dscribe
```
