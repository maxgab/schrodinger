#!/usr/bin/env python3
"""testing environment"""
import numpy as np
import pytest
import index

_TEST_DIRECTORIES = ["inf_well", "asym_well", "double_well_cub", "fin_well",
                     "harm_osz", "double_well_linear"]

@pytest.mark.parametrize("inp_dir", _TEST_DIRECTORIES)
def test_main(inp_dir):
    """test function for solver output"""
    path = 'solutions/' + inp_dir
    index.main(path + '/schrodinger.inp')

    test_pot = np.genfromtxt('output/potential.dat')
    test_eigen = np.genfromtxt('output/eigenvalues.dat')
    test_wave = np.genfromtxt('output/wavefuncs.dat')
    test_expv = np.genfromtxt('output/expvalues.dat')

    solution_pot = np.genfromtxt(path + '/potential.dat')
    solution_eigen = np.genfromtxt(path + '/eigenvalues.dat')
    solution_wave = np.genfromtxt(path + '/wavefuncs.dat')
    solution_expv = np.genfromtxt(path + '/expvalues.dat')

    assert  np.allclose(test_pot, solution_pot, 0.1)
    assert  np.allclose(test_eigen, solution_eigen, 0.1)
    assert  np.allclose(test_wave, solution_wave, 0.5)
    assert  np.allclose(test_expv, solution_expv, 0.03)
