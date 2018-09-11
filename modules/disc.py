"""tools for interpolation and discretization"""
#!/usr/bin/env python3
import numpy as np
from scipy import interpolate
def interp(pot_values, interp_type):
    """interpolates given potential
    input:
        pot_values: potential values in x-V(x) pairs (list)
        type: type of interpolation
    output:
        interpolated potential or false if failure
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
    input:
        pot_values: potential values function
        minval: x-min point
        maxval: x-max point
        disc_points: number of discretization points
    output:
        discretized potential (list)
    """
    pot_disc = []
    pot_disc.append(np.linspace(minval, maxval, disc_points))
    pot_disc.append(pot_values(pot_disc[0]))

    pot_disc = list(zip(pot_disc[0], pot_disc[1]))
    return pot_disc
