#!/usr/bin/env python3
"""This visualizer is used to plot the results created by the solver"""
import matplotlib.pyplot as plt
import numpy as np

def main():
    """Main function includes the top level code """

S1 = "input/schrodinger.inp"
with open(S1, "r") as fp:
# opens and reads the input file
    T1 = fp.readlines()
    S1 = {"x":list(np.loadtxt(T1[2:3]))}
#extracts the information on how many eigenvalues are wanted


A1 = int(((S1["x"])[1])-(S1["x"])[0]) + 1
#saves the number of eigenvalues that are wanted

X1 = np.genfromtxt('output/potential.dat', usecols=(0))
Y1 = np.genfromtxt('output/potential.dat', usecols=(1))
X2 = np.zeros(A1)
#creates zeros as x-values for the corresponding eigenvalues
Y2 = np.genfromtxt('output/eigenvalues.dat', usecols=(0))
#X3 = np.genfromtxt('erwartungswerte.dat', usecols=(0))
#Y3 = np.genfromtxt('erwartungswerte.dat', usecols=(1))
#Reads and saves the x and y values for potential erwartungswerte und eigenvalues

plt.subplot(121)
#creates the first subplot
plt.plot(X1, Y1, color='black', label='Potential')
plt.scatter(X2, Y2, color='green', linewidth='2.0', marker='x', label='Eigenvalues')
#plots the arrays in to one plot

WF = np.genfromtxt('output/wavefuncs.dat')
X4 = WF[:, 0]


CO = 0


for i in range(1, A1 + 1):

    energy = Y2[i-1]
    currentColumn = WF[:, i]
    for j in range(0, currentColumn.size):
        currentColumn[j] = 3*currentColumn[j] + energy
    plt.plot(X4, currentColumn, color='blue' if CO % 2 == 0 else 'red')
    CO = CO + 1
#reads and saves the x- and y-values for the wave functions
#creaes plot for the waves
#defines the x-axis
plt.xlabel('x[Bohr]')
#names the y-axis
plt.ylabel('Energy[Hardtree]')

for i in range(0, A1):
    plt.axhline(y = Y2[i], linewidth = 0.1, color = 'black')
#sdets the title
plt.title('Potential, Eigenvalues, Wavefunctions')

plt.ylim(ymax=Y2[A1-1]+1) 
#schneidet noch unten ab und nicht oben
#creats a legend
plt.legend()


#plt.subplot(122)
#creates sublot nr2
#plt.scatter(X3, Y3, color='fuchsia', marker='+')
#plots the array into the second plot
#plt.title('$\sigma_x$')
#creates an title with greeek letters

plt.show()
#shows the plot
