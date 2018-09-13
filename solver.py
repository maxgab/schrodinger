"""tools for interpolation and discretization, tools for solving eigenvalue problem"""
#!/usr/bin/env python3
import numpy as np
import scipy as sp
from scipy import interpolate
def interp(pot_values, interp_type):
    """interpolates given potential

    Args:
        pot_values (list): potential values in x-V(x) pairs (list)
        type (string): type of interpolation (string)
    Returns:
        list: interpolated potential function
    """
    values = list(zip(*pot_values))
    if interp_type == "linear":
        pot_interp = interpolate.interp1d(values[0], values[1], "linear")
    elif interp_type == "cspline":
        pot_interp = interpolate.CubicSpline(values[0], values[1])
    elif interp_type == "polynomial":
        pot_interp = interpolate.BarycentricInterpolator(values[0], values[1])
    else:
        pot_interp = 0
    return pot_interp

def discretize(pot_values, minval, maxval, disc_points):
    """discretizes given function

    Args:
        pot_values (function): potential values function
        minval (float): x-min point
        maxval (float): x-max point
        disc_points (int): number of discretization points
    Returns:
        list: discretized potential (list)
    """
    pot_disc = []
    pot_disc.append(np.linspace(minval, maxval, disc_points))
    pot_disc.append(pot_values(pot_disc[0]))

    pot_disc = list(zip(pot_disc[0], pot_disc[1]))
    return pot_disc

def solve_sgl(pot, mass, eigenv_range):
    """solves eigenvalue problem for discretized potential

    Args:
        pot (list): discretized potential
        mass (float): mass of particle
        range (list): first and last eigenvalue
    Returns:
        list: eigenvalues and eigenfunctions [[v, x1, x2, ...] ...]
    """
    pot = list(zip(*pot))
    distance = pot[0][1] - pot[0][0]
    factor = 1/(mass*(distance)**2)
    diag = factor + pot[1]
    offdiag = [-0.5*factor]*(len(pot[1]) - 1)
    result = sp.linalg.eigh_tridiagonal(diag, offdiag, False, 'i', (int(eigenv_range[0])-1,
                                                                    int(eigenv_range[1])-1))
    eigenval = [result[0], np.c_[pot[0], result[1]]]
    return eigenval

def normalize(vec):
    """normalizes a given wavefunction

    Args:
        wavefunc (list): wavefunction
    Returns:
        numpy array: normalized wavefunction [x][wfx]
    """
    vec = np.array(vec)
    dist = abs(vec[0][1] - vec[0][0])
    norm_sum = 0
    for item in vec[1]:
        norm_sum = norm_sum + abs(item)**2
    norm = (norm_sum * dist)**(0.5)

    vec[1] = vec[1] / norm
    return vec

def expected(func):
    """calculates expected value of operator

    Args:
        func (list): observalble list([list(x), list(val)])
    Returns:
        numpy array: expected value
    """
    func = np.array(func)
    dist = abs(func[0][1] - func[0][0])
    exp = (func[1]**2)*func[0]
    exp_sum = sum(exp)
    return exp_sum*dist

def uncertainty(func):
    """calculates heisenberg uncertainity of operator

    Args:
        func (list): normalized observable function [list(x), list(val)]
    Returns:
        numpy array: uncertainity
    """
    func = np.array(func)
    exp = expected(func)
    func[0] = func[0]**2
    exp_sq = expected(func)
    return (exp_sq - exp**2)**(0.5)
