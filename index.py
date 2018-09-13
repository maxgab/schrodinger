#!/usr/bin/env python3
"""solves the schrodinger equation for given parameters in /input/schrodinger.inp"""
import os
import numpy as np
import solver as sv
import visualizer as vz

def main():
    """Main function, includes top level code"""

    #check if schrodinger.inp exists
    inp_file = "schrodinger.inp"
    while os.path.isfile(inp_file) == 0:
        inp_file = input("Cannot find file schrodinger.inp. Please input file location: ")

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
    pot_values = list(np.loadtxt(txt[5:]))

    #interpolate and discretize potantial
    pot_values = sv.interp(pot_values, inp_data["interpolation"])
    pot_values = sv.discretize(pot_values, inp_data["range"][0],
                               inp_data["range"][1], inp_data["range"][2])

    #solve eigenvalue problem
    eigenvalues = sv.solve_sgl(pot_values, inp_data["mass"], inp_data["eigenvalues"])
    np.savetxt("output/potential.dat", pot_values)
    np.savetxt("output/eigenvalues.dat", eigenvalues[0])
    np.savetxt("output/wavefuncs.dat", eigenvalues[1])

    #additional values
    expvalues = []
    funcs = list(zip(*eigenvalues[1]))
    for item in funcs[1:]:
        norm = sv.normalize([funcs[0], item])
        expected = sv.expected(norm)
        uncert = sv.uncertainty(norm)
        expvalues.append([expected, uncert])

    np.savetxt("output/expvalues.dat", expvalues)

    vz.display()

main()
