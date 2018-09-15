"""tools for interpolation and discretization, tools for solving eigenvalue problem"""
#!/usr/bin/env python3
import numpy as np
import scipy as sp
from scipy import interpolate
def discretize(values, interp_type, minval, maxval, disc_points):
    """interpolates and discretizes given potential

    Args:
        pot_values (list): potential values [[x], [V(x)]]
        type (string): type of interpolation
        minval (float): x-min point
        maxval (float): x-max point
        disc_points (int): number of discretization points
    Returns:
        list: discretized potential [[x], [V(x)]]
    """

    if interp_type == "linear":
        pot_values = interpolate.interp1d(values[0], values[1], "linear")
    elif interp_type == "cspline":
        pot_values = interpolate.CubicSpline(values[0], values[1])
    elif interp_type == "polynomial":
        pot_values = interpolate.BarycentricInterpolator(values[0], values[1])
    else:
        pot_values = 0

    pot_disc = []
    pot_disc.append(np.linspace(minval, maxval, disc_points))
    pot_disc.append(pot_values(pot_disc[0]))

    return pot_disc

def solve_sgl(pot, mass, eigenv_range):
    """solves eigenvalue problem for discretized potential

    Args:
        pot (list): discretized potential [[x], [V(x)]]
        mass (float): mass of particle
        range (list(float)): first and last eigenvalue [min, max]
    Returns:
        list: eigenvalues and eigenfunctions [[eigenv], [[x1-wf1, x1-wf2, ...], ...]]
    """
    distance = abs(pot[0][1] - pot[0][0])
    factor = 1/(mass*(distance)**2)
    diag = factor + pot[1]
    offdiag = [-0.5*factor]*(len(pot[1]) - 1)
    (eigenval, eigenfunc) = sp.linalg.eigh_tridiagonal(diag, offdiag, False, 'i',
                                                       (int(eigenv_range[0])-1,
                                                        int(eigenv_range[1])-1))
    eigenfunc = eigenfunc.T
    for index, item in enumerate(eigenfunc):
        eigenfunc[index] = normalize(item, distance)

    return [eigenval, eigenfunc]

def normalize(vec, dist):
    """normalizes a given vector

    Args:
        vec (list): vector
        dist (float): distance between x points
    Returns:
        numpy array: normalized wavefunction [x][wfx]
    """
    vec = np.array(vec)
    norm_sum = 0
    for item in vec:
        norm_sum += item**2
    norm = (norm_sum * dist)**(0.5)

    vec = vec / norm
    return vec

def expected(func, dist):
    """calculates expected value of operator

    Args:
        func (list): vector [[x], [val]]
        dist (float): distance between x points
    Returns:
        numpy array: expected value
    """
    func = np.array(func)
    exp = (func[1]**2)*func[0]
    exp_sum = sum(exp)
    return exp_sum*dist

def uncertainty(func):
    """calculates heisenberg uncertainty

    Args:
        func (list): normalized observable function [[x], [val])]
    Returns:
        numpy array: uncertainity
    """
    func = np.array(func)
    dist = abs(func[0][1] - func[0][0])
    exp = expected(func, dist)
    func[0] = func[0]**2
    exp_sq = expected(func, dist)
    return (exp_sq - exp**2)**(0.5)
