import os
import ase.io
from ase.data import chemical_symbols

TURBOGAP_CMD = "turbogap"

def write_turbogap_input_files(filename, system_filename, system, desc, atomic_numbers_to_indices, idx=None):
    ase.io.write(system_filename, system)

    def write_iterable(f, name, iterable):
        f.write(f"{name} = ")
        for it in iterable:
            f.write(f"{it} ")
        f.write("\n")
    species = len(atomic_numbers_to_indices)*[None]
    assert len(species) == desc['num_species']
    for k,v in atomic_numbers_to_indices.items():
        species[v] = chemical_symbols[k]

    with open(filename, 'w') as f: 
        f.write(f"input_file = {system_filename}\n")
        f.write(f"num_species = {desc['num_species']}\n")
        write_iterable(f, "species", species )
        rcuts = desc['rcut_hard']
        buffers = rcuts - desc['rcut_soft']
        write_iterable(f, "rcut", rcuts)
        write_iterable(f, "buffer", buffers)
        write_iterable(f, "atom_sigma_r", desc['atom_sigma_r'])
        write_iterable(f, "atom_sigma_t", desc['atom_sigma_t'])
        write_iterable(f, "atom_sigma_r_scaling", desc['atom_sigma_r_scaling'])
        write_iterable(f, "atom_sigma_t_scaling", desc['atom_sigma_t_scaling'])
        write_iterable(f, "amplitude_scaling", desc['amplitude_scaling'])
        write_iterable(f, "global_scaling", desc['global_scaling'])
        write_iterable(f, "n_max", desc['nmax'])
        f.write(f"l_max = {desc['lmax']}\n")
        write_iterable(f, "nf", desc['nf'])
        write_iterable(f, "central_weight", desc['central_weight'])
        if idx is not None:
            f.write(f"which_atom = {idx+1}\n")

        f.write(f"ase_format = .true.\nscaling_mode = {desc['scaling_mode']}\nbasis = {desc['basis']}\ntiming = .true.\n")

def run_turbogap(system, soap, idx=None, turbogap_cmd=TURBOGAP_CMD):
    write_turbogap_input_files("input", "input.xyz", system, soap, idx=idx)
    subprocess.run([turbogap_cmd], shell=True, check=True)
    result = np.loadtxt("soap.dat", skiprows=1)
    os.remove("input")
    os.remove("soap.dat")
    os.remove("input.xyz")
    return result