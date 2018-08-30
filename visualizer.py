#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np


"""visualizes the results from solver.py given in the output folder"""
def main():
	"""Main function, includes top level code"""
    #reads the files from the solvers and puts it in to an array
    x1, y1 = np.loadtxt('potential.dat', delimiter=',', unpack=True)
    x2, y2 = np.loadtxt('eigenvalues.dat', delimiter=',', unpack=True)
    x3, y3 = np.loadtxt('wavefuncs.dat', delimiter=',', unpack=True)
    

	#plots the arrays in to one plot	                 
    plt.plot(x1,y1, color='black', label='Potential')
    plt.scatter(x2,y2, color='green', linewidth='2.0', marker='x' label='Eigenvalues')
    plt.plot(x3,y3, color='blue', label='Wavefunctions')
    
	#defines the x-axis
    plt.xlabel('x[Bohr]') 
	#names the y-axis
    plt.ylabel('Energy[Hardtree]')
	#sdets the title
    plt.title('Potential, Eigenvalues, Wavenfunctions')
	#creats a legend
    plt.legend()
	#shows the plot
    plt.show()

main()
                     
                     
    
