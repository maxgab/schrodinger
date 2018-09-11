"""tools for solving eigenvalue problem"""
#!/usr/bin/env python3
import scipy as sp
import numpy as np
def eigenval_sgl(pot, mass, eigenv_range):
    """solves eigenvalue problem for discretized potential
    input:
        pot: discretized potential
        mass: mass of particle
        range: first and last eigenvalue (list)
    output:
        list of eigenvalues and eigenfunctions [[v, x1, x2, ...] ...]
    """
    pot = list(zip(*pot))
    distance = pot[0][1] - pot[0][0]
    factor = 1/(mass*distance**2)
    diag = factor + pot[1]
    offdiag = [-0.5*factor]*(len(pot[1]) - 1)
    result = sp.linalg.eigh_tridiagonal(diag, offdiag, False, 'i', (int(eigenv_range[0])-1,
                                                                    int(eigenv_range[1])-1))
    eigenval = [result[0], np.c_[pot[0], result[1]]]
    return eigenval

def normalize(wavefunc):
    """normalizes a given wavefunction
    input:
        wavefunc: wavefunction (list[x][wfx])
    output:
        normalized wavefunction (list[x][wfx])
    """
    dist = abs(wavefunc[0][1] - wavefunc[0][0])
    norm_sum = 0
    for item in wavefunc[1]:
        norm_sum = norm_sum + abs(item)**2
    norm = (norm_sum * dist)**(0.5)
    wavefunc[1] = wavefunc[1] / norm
    return wavefunc

def expected(obs_func):
    """calculates expected value of operator
    input:
        observalble (list[x][wfx])
    output:
        expected value
    """
    dist = abs(obs_func[0][1] - obs_func[0][0])
    obs_sum = 0
    for index, item in enumerate(obs_func[1]):
        obs_sum = obs_sum + item**2 * obs_func[0][index]
    return dist * obs_sum

def uncertainity(obs_func):
    """calculates heisenberg uncertainity of operator
    input:
        normalized observable function (list)
    output:
        uncertainity
    """
    obs = expected(obs_func)
    print("-------")
    print(obs)
    obs_func[0] = [x**2 for x in obs_func[0]]
    obs_sq = expected(obs_func)
    print(obs_sq)
    return (obs_sq - obs**2)**(0.5)
