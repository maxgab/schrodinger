import index
import pytest
import numpy as np

_TEST_DIRECTORIES = ["inf_well", "asym_well", "double_well_cub", "fin_well", "harm_osz", "double_well_linear"]

@pytest.mark.parametrize("dir", _TEST_DIRECTORIES)
def test_main(dir):
    path = 'solutions/' + dir
    index.main(path + '/schrodinger.inp')

    test_pot = np.genfromtxt('output/potential.dat')
    test_eigen = np.genfromtxt('output/eigenvalues.dat')
    test_wave = np.genfromtxt('output/wavefuncs.dat')
    test_expv = np.genfromtxt('output/expvalues.dat')

    solution_pot = np.genfromtxt( path + '/potential.dat')
    solution_eigen = np.genfromtxt( path + '/eigenvalues.dat')
    solution_wave = np.genfromtxt( path + '/wavefuncs.dat')
    solution_expv = np.genfromtxt( path + '/expvalues.dat')
    
    assert  abs(np.all(test_pot - solution_pot)) < 0.1
    assert  abs(np.all(test_eigen - solution_eigen)) < 0.1
    assert  abs(np.all(test_wave - solution_wave))  < 0.5
    assert  abs(np.all(test_expv - solution_expv)) < 0.03