#!/usr/bin/env python3
"""This visualizer is used to plot the results created by the solver"""
import matplotlib.pyplot as plt
import numpy as np

def display():
    """Main function includes the top level code """

    s1 = "schrodinger.inp"
    with open(s1, "r") as fp:
    # opens and reads the input file
        t1 = fp.readlines()
        s1 = {"x":list(np.loadtxt(t1[2:3]))}
    #extracts the information on how many eigenvalues are wanted


    a1 = int(((s1["x"])[1])-(s1["x"])[0]) + 1
    #saves the number of eigenvalues that are wanted

    x1 = np.genfromtxt('output/potential.dat', usecols=(0))
    y1 = np.genfromtxt('output/potential.dat', usecols=(1))
    x2 = np.genfromtxt('output/expvalues.dat', usecols=(0))
    #creates zeros as x-values for the corresponding eigenvalues
    y2 = np.genfromtxt('output/eigenvalues.dat', usecols=(0))
    x3 = np.genfromtxt('output/expvalues.dat', usecols=(1))
    y3 = np.genfromtxt('output/eigenvalues.dat', usecols=(0))
    #Reads and saves the x and y values for potential erwartungswerte und eigenvalues

    plt.subplot(121)
    #creates the first subplot
    plt.plot(x1, y1, color='black', label='Potential')
    plt.scatter(x2, y2, color='green', linewidth='2.0', marker='x', label='Expected values')
    #plots the arrays in to one plot

    wf = np.genfromtxt('output/wavefuncs.dat')
    x4 = wf[:, 0]


    co = 0


    for i in range(1, a1 + 1):

        energy = y2[i-1]
        currentcolumn = wf[:, i]
        for j in range(0, currentcolumn.size):
            currentcolumn[j] = 3*currentcolumn[j] + energy
        plt.plot(x4, currentcolumn, color='blue' if co % 2 == 0 else 'red')
        co = co + 1
    #reads and saves the x- and y-values for the wave functions
    #creaes plot for the waves
    #defines the x-axis
    plt.xlabel('x[Bohr]')
    #names the y-axis
    plt.ylabel('Energy[Hardtree]')

    for i in range(0, a1):
        plt.axhline(y=y2[i], linewidth=0.1, color='black')
    #sdets the title
    plt.title('Potential, Expected values, Wavefunctions')

    plt.ylim(ymax=y2[a1-1]+1)
    #schneidet noch unten ab und nicht oben
    #creats a legend
    plt.legend()


    plt.subplot(122)
    #creates sublot nr2
    plt.scatter(x3, y3, color='fuchsia', marker='+')
    #plots the array into the second plot
    plt.title('$\sigma_x$')
    #creates an title with greeek letters

    plt.show()
    #shows the plot
