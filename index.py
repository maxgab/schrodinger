#!/usr/bin/env python3
"""solves the schrodinger equation for given parameters"""
import os
import numpy as np
import solver as sv
import visualizer as vz
import custom as ct

def main(inp_file):
    """Main function, includes top level code"""

    #check if schrodinger.inp exists
    inp_file = ct.file_input(inp_file)

    #check if output files exist
    outp_files = ["eigenvalues.dat", "wavefuncs.dat", "potential.dat", "expvalues.dat"]
    outp_files = ["output/" + i for i in outp_files]
    os.makedirs(os.path.dirname("output/"), exist_ok=True)
    for file_name in outp_files:
        file = open(file_name, "a")
        file.close()

    #read from schrodinger.inp
    with open(inp_file, "r") as fp:
        txt = fp.readlines()
        inp_data = {
            "mass": float(np.loadtxt(txt[0:1])),
            "range": list(np.loadtxt(txt[1:2])),
            "eigenvalues": list(np.loadtxt(txt[2:3])),
            "interpolation": str(np.loadtxt(txt[3:4], "str")),
            "points": float(np.loadtxt(txt[4:5]))
        }
    #pot_values = list(np.loadtxt(txt[5:]))
    pot_values = np.loadtxt(txt[5:]).T
    #pot_values = list(zip(*pot_values))

    #interpolate and discretize potantial
    pot_values = sv.discretize(pot_values, inp_data["interpolation"], inp_data["range"][0],
                               inp_data["range"][1], inp_data["range"][2])

    #solve eigenvalue problem
    (eigenvalues, eigenfuncs) = sv.solve_sgl(pot_values, inp_data["mass"], inp_data["eigenvalues"])
    np.savetxt("output/potential.dat", list(zip(pot_values[0], pot_values[1])))
    np.savetxt("output/eigenvalues.dat", eigenvalues)
    np.savetxt("output/wavefuncs.dat", np.c_[pot_values[0], eigenfuncs.T])

    #additional values
    expvalues = []
    dist = abs(pot_values[0][1] - pot_values[0][0])
    for item in eigenfuncs:
        expected = sv.expected([pot_values[0], item], dist)
        uncert = sv.uncertainty([pot_values[0], item])
        expvalues.append([expected, uncert])

    np.savetxt("output/expvalues.dat", expvalues)
    return inp_file

if __name__ == "__main__":
    inp_file = main("schrodinger.inp")
    vz.display(inp_file)
