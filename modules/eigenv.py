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
