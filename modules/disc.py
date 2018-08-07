"""tools for interpolation and discretization"""
#!/usr/bin/env python3
import numpy as np
from scipy import interpolate
def interp(pot_values, interp_type):
    """interpolates given potential
    input:
        pot_values: potential values in x-V(x) pairs
        type: type of interpolation
    return:
        interpolated potential or false if failure
    """
    values = list(zip(*pot_values))
    if interp_type == "linear":
        pot_interp = interpolate.interp1d(values[0], values[1], "linear")
    elif interp_type == "cspline":
        pot_interp = interpolate.interp1d(values[0], values[1], "cubic")
    elif interp_type == "polynomial":
        pot_interp = interpolate.BarycentricInterpolator(values[0], values[1])
    else:
        pot_interp = 0
    return pot_interp

def discretize(pot_values, minval, maxval, disc_points):
    """discretizes given function"""
    pot_disc = []
    pot_disc.append(np.linspace(minval, maxval, disc_points))
    pot_disc.append(pot_values(pot_disc[0]))

    #import matplotlib.pyplot as plt
    #plt.plot(pot_disc[0], pot_disc[1])
    #plt.show()

    pot_disc = list(zip(pot_disc[0], pot_disc[1]))
    return pot_disc
