#!/usr/bin/env python

import f
import J
import chi
import entropy
import newton
import printFile
import nan
import model
import solver

def read_spectral(fileName):
    import numpy as np
    import sys

    omega = []
    A = []
    try:
        ifile = open(fileName, "r")
    except:
        sys.exit(fileName + " does not exist. ")
    for index, string in enumerate(ifile):
        a = string.split()
        omega.append(float(a[0]))
        A.append(float(a[1]))
    ifile.close()
    return omega, np.asarray(A)

def readFiles(Greal, Gimag):
    import os
    import sys
    import numpy as np
    
    omega_n = []
    G_real = []
    G_imag = []
    
    try:
        ifile = open(Greal, "r")
    except:
        sys.exit(Greal + " does not exist. ")
    for index, string in enumerate(ifile):
        a = string.split()
        omega_n.append(float(a[0]))
        G_real.append(float(a[1]))
    ifile.close()

    try:
        ifile = open(Gimag, "r")
    except:
        sys.exit(Gimag + " does not exist. ")
    for index, string in enumerate(ifile):
        a = string.split()
        G_imag.append(float(a[1]))
    ifile.close()

    G_real = np.asarray(G_real)
    G_imag = np.asarray(G_imag)
    
    return omega_n, G_real, G_imag

def norm(vector):
    import numpy as np
    return np.sqrt(vector.dot(vector))

def main():
    import os
    import sys
    import numpy as np

    if (len(sys.argv) == 1):
        print "a0 = sys.argv[1], b0 = sys.argv[2]. alpha = a0*exp(-i*b0). "
        a0 = 10000
        b0 = 0.05
        print "The default value of a0 =", a0, " and b0 = ", b0, " has been chosen. "
    if (len(sys.argv) == 2):
        a0 = float(sys.argv[1])
        b0 = 0.05
    if (len(sys.argv) == 3):
        a0 = float(sys.argv[1])
        b0 = float(sys.argv[2])
    
    Greal = "G_real_rotated.txt"
    Gimag = "G_imag_rotated.txt"
    omega_n, G_real_rotated, G_imag_rotated = readFiles(Greal, Gimag)
    Niom = len(omega_n)

    if (not os.path.exists("model.txt")):
        mu = 0
        sigma = 1
        model.generateModel(mu, sigma)

    if (not os.path.exists("A_initial.txt")):
        os.system("cp model.txt A_initial.txt")
    omega, A_initial = read_spectral("A_initial.txt")
    Nomega = len(omega)
    D = np.zeros(Nomega)
    for i in range(Nomega):
        D[i] = A_initial[i]

    K_real_rotated = np.zeros((Niom, Nomega))
    K_imag_rotated = np.zeros((Niom, Nomega))
    ifile = open("K_real_rotated.txt", "r")
    for i in range(Niom):
        for j in range(Nomega):
            string = ifile.readline()
            a = string.split()
            K_real_rotated[i,j] = float(a[2])
    ifile.close()
    ifile = open("K_imag_rotated.txt", "r")
    for i in range(Niom): 
        for j in range(Nomega):
            string = ifile.readline()
            a = string.split()
            K_imag_rotated[i,j] = float(a[2])
    ifile.close()
    
    singular_values = []
    ifile = open("s.txt", "r")
    for index, string in enumerate(ifile):
        singular_values.append(float(string))
    ifile.close()
    if (len(singular_values) != Niom):
        sys.exit("Number of singular values should be equal to the number of Matsubara frequencies. ")
    LambdaInverse = np.zeros((Niom, Niom))
    for i in range(Niom):
        LambdaInverse[i,i] = 1.0/singular_values[i]

    if (True):
        alpha = []
        error = []
        for i in range(30):
            alpha.append(a0*np.exp(-i*b0))
        
        ofile = open("alpha.txt", "a")
        for i in range(len(alpha)):
            #A_updated = newton.newton(alpha[i], G_real_rotated, G_imag_rotated, K_real_rotated, K_imag_rotated, omega, A_initial, D, LambdaInverse)
            A_updated = solver.solver(alpha[i], G_real_rotated, G_imag_rotated, K_real_rotated, K_imag_rotated, omega, A_initial, D, LambdaInverse)
            if (nan.array_isnan(A_updated)):
                omega, A_initial = read_spectral("A_initial.txt")
                continue
            output = "A_updated_alpha_" + str(alpha[i]) + ".txt"
            printFile.printFile(omega, A_updated, output)
            os.system("cp " +  output + " A_initial.txt")
            ofile.write(str(alpha[i]) + "\n")
            diff = norm(A_updated - A_initial)
            error.append(diff)
            print "alpha = ", alpha[i], ", error = ", diff
            A_initial = A_updated
        ofile.close()
        printFile.printFile(alpha, error, "error_alpha.txt")
        os.system("python pltfiles.py A_initial.txt")

    return 0

main()
