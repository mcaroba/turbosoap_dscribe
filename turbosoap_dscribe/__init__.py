# -*- coding: utf-8 -*-
"""
Copyright (c) 2018-2020 by Miguel A. Caro and others.

Licensed under the Creative Commons Attribution-NonCommercial-ShareAlike
4.0 International Public License (the "Licence"); 
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode.

Read the README file for up-to-date information on how to appropriately give attribution to
the code author(s). If you want to use the TurboGAP code or parts thereof for commercial purposes, contact the
copyright holder, Miguel A. Caro (mcaroba@gmail.com).
"""

from dscribe.utils.geometry import get_adjacency_matrix, get_extended_system
import numpy as np
import scipy.sparse
from math import sqrt, acos, atan2
import turbosoap_ext

#TurboSOAPSpecie is used to define per-species parameters
#prepare_turbosoap_configuration will then compile TurboSOAPSpecie
#into an format suitable for fortran interface.
#For consistency all the error checking is done in prepare_turbosoap_configuration
#where all the information is available.

def prepare_turbosoap_configuration(species, lmax = 8):
    """Checks configuration parameters of each species and prepares them for fortran interface
    Order is significant. Same species given in different order correspond
    to different descriptor.
    Args:
        lmax (int): The maximum degree of spherical harmonics.
        species (iterable): List of TurboSOAPSpecie 
    Returns:
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
            f"Rcut should be bigger than 1 angstrom for species {specie.debug_name}"
            )
        config['rcut_hard'][i] = specie.rcut
        if specie.buffer < 0.0 or specie.buffer > specie.rcut:
            raise ValueError(
                f"Buffer must be positive and smaller than rcut for species {specie.debug_name}. Buffer: {specie.buffer}"
            )
        config['rcut_soft'][i] = specie.rcut - specie.buffer
        if specie.nmax < 1:
            raise ValueError(
            f"Must have at least one radial basis function for species {specie.debug_name}. nmax={specie.nmax}"
            )
        if specie.nmax > 12:
            raise ValueError(
            f"Too many basis functions for species {specie.debug_name}. Orthonormalization is numerically unstable." 
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
    rcuts = [atomic_numbers_to_rcuts[ n ] for n in system.numbers]
    max_rcut = max(rcuts)
    if periodic:
        original_system = system
        system = get_extended_system(system, max_rcut, return_cell_indices=False)
        rcuts = [atomic_numbers_to_rcuts[ n ] for n in system.numbers]
    n_sites = len(system)
    species = [atomic_numbers_to_indices[ n ] for n in system.numbers]
    f_species = np.empty((1,n_sites), dtype=int, order='F')
    f_species[0,:] = np.array(species, dtype=int) + 1 #fortran indexing from 1
    assert len(atomic_numbers_to_indices) == len(config['rcut_hard'])
    n_species = len(atomic_numbers_to_indices)
    species_multiplicity = np.ones(n_sites, dtype=int)
    positions = system.positions
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
    soap_cart_der = np.empty((3, n_soap, n_atom_pairs), order='F')
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
    if periodic:
        #We "map" results back to original systems. We assume that the order of the atoms is unchanged i.e. the original unit cells atoms are first
        num_original_sites = len(original_system)
        soap_m = soap_m[0:num_original_sites]
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
