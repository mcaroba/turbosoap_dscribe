from dscribe.utils.geometry import get_adjacency_matrix
import numpy as np
import scipy.sparse
from math import sqrt, acos, atan2
import turbosoap_ext

#Symbol is 
class TurboSOAPSpecie:
    """Helper class for all the TurboSOAP per-species parameters
    """
    def __init__(self, rcut, nmax = 8, buffer = 0.5,
        atom_sigma_r = 0.5, atom_sigma_t = 0.5,
        atom_sigma_r_scaling = 0., atom_sigma_t_scaling = 0.,
        radial_enhancement = 0, amplitude_scaling = 1.0,
        central_weight = 1., global_scaling = 1., nf = 4.):
        """
        *************************************************************************************************
        Args:
            rcut (float):                    A cutoff for local region in angstroms. Should be bigger
                                             than the buffer.

            nmax (int):                      Number of radial basis functions.

            buffer (float):                  Width of buffer region where atomic density field will decay
                                             to zero

            atom_sigma_r (float):            Width of radial sigma (in Angstrom) at the origin

            atom_sigma_r_scaling (float):    Radial-scaling parameter for radial sigma (dimensionless)

            atom_sigma_t (float):            Width of angular sigma (in Angstrom) at the origin

            atom_sigma_t_scaling (float):    Radial-scaling parameter for angular sigma (dimensionless)

            radial_enhancement (integer):    Distant atomic densities get weighted by the integral of
                                             a Gaussian located at the same position times the radial
                                             coordinated raised to this power. It must be 0, 1 or 2

            amplitude_scaling (float):       The atomic densities get weighted, according to the position
                                             of their atomic centers, by a radial function that decays to
                                             zero at the cutoff according to a power law, where this
                                             parameter is the exponent

            global_scaling (float):          Multiplicative factor for the whole atomic density fieldç
                                             corresponding to this species

            central_weight (float):          Scaling factor for the atomic density field corresponding to
                                             the central atom, when the centra latom is of this species

            nf (float):                      Decay parameter of the exponential within the buffer region
        *************************************************************************************************
        """
        if rcut <= buffer:
            raise ValueError(
                "Rcut should be bigger than the buffer region for species {symbol}"
            )
        if buffer <= 0.:
            raise ValueError(
                "The buffer region should be greater than 0 for species {symbol}"
            )
        self.rcut = rcut
        if nmax < 1 or nmax > 12:
            raise ValueError
        self.nmax = nmax
        self.buffer = buffer
        self.atom_sigma_r = atom_sigma_r
        self.atom_sigma_t = atom_sigma_t
        self.atom_sigma_r_scaling = atom_sigma_r_scaling
        self.atom_sigma_t_scaling = atom_sigma_t_scaling
        self.radial_enhancement = radial_enhancement
        self.amplitude_scaling = amplitude_scaling
        self.global_scaling = global_scaling
        self.central_weight = central_weight
        self.nf = nf

default_turbosoap_specie = TurboSOAPSpecie(5.0, 8)

def prepare_turbosoap_configuration(lmax = 8, species):
    """Checks configuration parameters of each species and prepares them for fortran interface
    Order is significant. Same species given in different order correspond
    to different descriptor.
    Args:
        lmax (int): The maximum degree of spherical harmonics.
        species (iterable): List of TurboSOAPSpecie 
    Return:
        Configuration values and np.ndarray which are in correct format 
        for turbosoap fortran interface
    """
    config = {}
    #Global values
    if lmax < 0:
        raise ValueError(f"lmax cannot be negative. lmax={lmax}")
    if lmax > 25:
        raise ValueError(f"Maximum lmax is 25. lmax={lmax}")
    config['lmax'] = lmax
    config['basis'] = 'poly3gauss'
    config['scaling_mode'] = 'polynomial'
    num_species = len(species)
    config['num_species'] = num_species

    #Per species values
    config['rcut_hard'] = np.empty(num_species, dtype=float)
    config['rcut_soft'] = np.empty(num_species, dtype=float)
    config['nmax'] = np.empty(num_species, dtype=int)
    config['atom_sigma_r'] = np.empty(num_species, dtype=float)
    config['atom_sigma_r_scaling'] = np.empty(num_species, dtype=float)
    config['atom_sigma_t'] = np.empty(num_species, dtype=float)
    config['atom_sigma_t_scaling'] = np.empty(num_species, dtype=float)
    config['amplitude_scaling'] = np.empty(num_species, dtype=float)
    config['radial_enhancement'] = np.empty(num_species, dtype=int)
    config['nf'] = np.empty(num_species, dtype=float)
    config['global_scaling'] = np.empty(num_species, dtype=float)
    config['central_weight'] = np.empty(num_species, dtype=float)
    for i,specie in enumerate(species):
        if specie.rcut <= 1.0:
            raise ValueError(
            f"Rcut should be bigger than 1 angstrom for species {specie.rcut}"
            )
        config['rcut_hard'][i] = specie.rcut
        if specie.buffer < 0.0 or specie.buffer > specie.rcut:
            raise ValueError(
                f"Buffer must be positive and smaller than rcut. Buffer: {specie.buffer}"
            )
        config['rcut_soft'][i] = specie.rcut - specie.buffer
        if specie.nmax < 1:
            raise ValueError(
            f"Must have at least one radial basis function. nmax={specie.nmax}"
            )
        if specie.nmax > 12:
            raise ValueError(
            f"Too many basis functions. Orthonormalization is numerically unstable." 
            "nmax={specie.nmax}."
            )
        config['nmax'][i] = specie.nmax
        config['atom_sigma_r'][i] = specie.atom_sigma_r
        config['atom_sigma_r_scaling'][i] = specie.atom_sigma_r_scaling
        config['atom_sigma_t'][i] = specie.atom_sigma_t
        config['atom_sigma_t_scaling'][i] = specie.atom_sigma_t_scaling
        config['amplitude_scaling'][i] = specie.amplitude_scaling
        if specie.radial_enhancement < 0 or specie.radial_enhancement > 2:
            raise ValueError(
                f"Radial enhancement must be 0, 1 or 2. Radial enhancement: {specie.radial_enhancement}."
            )
        config['radial_enhancement'][i] = specie.radial_enhancement
        config['nf'][i] = specie.nf
        config['global_scaling'][i] = specie.global_scaling
        config['central_weight'][i] = specie.central_weight


    n = sum(config['nmax'])
    config['num_components'] = n*(n+1)//2 * (config['lmax']+1)

    return config

def calculate_turbosoap_descriptor(config, system, periodic, 
    atomic_numbers_to_indices, atomic_numbers_to_rcuts):
    n_sites = len(system)
    rcuts = [atomic_numbers_to_rcuts[ n ] for n in system.numbers]
    species = [atomic_numbers_to_indices[ n ] for n in system.numbers]
    f_species = np.empty((1,n_sites), dtype=int, order='F')
    f_species[0,:] = np.array(species, dtype=int) + 1 #fortran indexing from 1
    assert len(atomic_numbers_to_indices) == len(config['rcut_hard'])
    max_rcut = max(rcuts)
    n_species = len(atomic_numbers_to_indices)
    species_multiplicity = np.ones(n_sites, dtype=int)
    positions = system.positions
    if periodic:
        system = get_extended_system(system, max_rcut, return_cell_indices=False)
    adj_m = get_adjacency_matrix(max_rcut, positions, system.positions)
    adj_l = get_adjacency_list_rcut(adj_m, rcuts)

    n_neigh = np.empty(n_sites, dtype=int)
    n_atom_pairs = 0
    for i,l in enumerate(adj_l):
        num_neigh = len(l)
        n_neigh[i] = num_neigh
        n_atom_pairs += num_neigh

    rjs = np.empty(n_atom_pairs, dtype=float)
    thetas = np.empty(n_atom_pairs, dtype=float)
    phis = np.empty(n_atom_pairs, dtype=float)
    mask = np.zeros((n_atom_pairs, n_species), order='F', dtype=np.int8)
    idx = 0
    for nl in adj_l:
        rjs[idx] = 0.0
        thetas[idx] = 0.0
        phis[idx] = 0.0
        mask[idx, species[nl[0]]] = True
        idx += 1
        cpos = positions[nl[0]]
        for n in nl[1:]:
            p = system.positions[n] - cpos
            d = p*p
            d = sqrt(d[0] + d[1] + d[2])
            rjs[idx] = d
            thetas[idx] = acos(p[2]/d)
            phis[idx] = atan2(p[1], p[0])
            mask[idx, species[n]] = True
            idx += 1
    assert idx == n_atom_pairs

    do_timing = False
    do_derivatives = False 
    n_soap = config['num_components']
    soap_m = np.zeros((n_soap, n_sites), order='F')
    soap_cart_der = np.zeros((3, n_soap, n_atom_pairs), order='F')
    #soap_cart_der = np.empty((1,1,1))
    if False:
        print(soap_m.shape)
        print(n_neigh)
        print(n_species)
        print(f_species)
        print(species_multiplicity)
        print(n_atom_pairs)
        print(mask)
        print(rjs)
        print(thetas)
        print(phis)
        print("Dictionary info:")
        for k,v in config.items():
            print(f"{k}: {v}")
            
    turbosoap_ext.soap_desc.get_soap(n_sites, n_neigh, n_species, f_species, 
            species_multiplicity, n_atom_pairs,
            mask, rjs, thetas, phis, config['nmax'], config['lmax'], 
            config['rcut_hard'], config['rcut_soft'], config['nf'], config['global_scaling'], 
            config['atom_sigma_r'], config['atom_sigma_r_scaling'], 
            config['atom_sigma_t'], config['atom_sigma_t_scaling'],
            config['amplitude_scaling'], config['radial_enhancement'], config['central_weight'], 
            config['basis'], config['scaling_mode'], do_timing,
            do_derivatives, soap_m, soap_cart_der)
    soap_m = np.ascontiguousarray(soap_m.T)
    return soap_m

#Multiple rcut aware version of adjacency list construction.
#Also aware of soap.f90 convention where central atom is at the center 
#Original in dscribe/utils/geometry.py
def get_adjacency_list_rcut(adjacency_matrix, rcuts):
    if type(adjacency_matrix) != scipy.sparse.coo_matrix:
        adjacency_matrix = adjacency_matrix.tocoo()
    assert adjacency_matrix.shape[0] == len(rcuts)
    adjacency_list = [[i] for i in range(adjacency_matrix.shape[0])]
    for i,j,d in zip(adjacency_matrix.row, adjacency_matrix.col, adjacency_matrix.data):
        if d <= rcuts[i] and i != j:
            adjacency_list[i].append(j)
    return adjacency_list
