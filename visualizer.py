#!/usr/bin/env python3
"""This visualizer is used to plot the results created by the solver"""
import matplotlib.pyplot as plt
import numpy as np
def display():
    """Main function includes the top level code """
    s1 = "schrodinger.inp"
    with open(s1, "r") as fp:
        t1 = fp.readlines()
        s1 = {"x":list(np.loadtxt(t1[2:3]))}
    a1 = int(((s1["x"])[1])-(s1["x"])[0]) + 1
    px = np.genfromtxt('output/potential.dat', usecols=(0))
    py = np.genfromtxt('output/potential.dat', usecols=(1))
    ex = np.genfromtxt('output/expvalues.dat', usecols=(0))
    ev = np.genfromtxt('output/eigenvalues.dat', usecols=(0))
    ey = np.genfromtxt('output/expvalues.dat', usecols=(1))
    a2 = np.min(py)
    wf = np.genfromtxt('output/wavefuncs.dat')
    x4 = wf[:, 0]
    ei = np.max(ev) - np.min(ev)
    yx = ev[a1-1] + ei/4
    yn = a2 - ei/10
    sf = 1 + (-4/7) * a1 + (40/7)

    plt.subplot(121)
    plt.plot(px, py, color='black', label='Potential')
    plt.scatter(ex, ev, color='green', linewidth='2.0', marker='x', label='Expected values')
    co = 0
    for i in range(1, a1 + 1):
        en = ev[i-1]
        cc = wf[:, i]
        for j in range(0, cc.size):
            cc[j] = sf * cc[j] + en
        plt.plot(x4, cc, color='blue' if co % 2 == 0 else 'red')
        co = co + 1
    plt.xlabel('x[Bohr]')
    plt.ylabel('Energy[Hardtree]')
    for i in range(0, a1):
        plt.axhline(y=ev[i], linewidth=0.1, color='black')
    plt.title('Potential, Expected values, Wavefunctions')
    plt.ylim(ymax=yx)
    plt.ylim(ymin=yn)
    plt.legend()
    plt.subplot(122)
    plt.scatter(ey, ev, color='fuchsia', marker='+')
    for i in range(0, a1):
        plt.axhline(y=ev[i], linewidth=0.1, color='black')
    plt.ylim(ymax=yx)
    plt.ylim(ymin=yn)
    plt.yticks([])
    plt.xlabel('x[Bohr]')
    plt.title('$\sigma_x$')
    plt.show()
