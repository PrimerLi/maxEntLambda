import chi
import entropy
import numpy as np

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


def readA(fileName):
    omega = []
    A = []
    ifile = open(fileName, "r")
    for i, string in enumerate(ifile):
        a = string.split()
        omega.append(float(a[0]))
        A.append(float(a[1]))
    ifile.close()
    return omega, A

def main():
    import os
    import sys
    import numpy.linalg
    
    Greal = "G_real_rotated.txt"
    Gimag = "G_imag_rotated.txt"
    omega_n, G_real_rotated, G_imag_rotated = readFiles(Greal, Gimag)
    Niom = len(omega_n)

    singular_values = []
    ifile = open("s.txt", "r")
    for index, string in enumerate(ifile):
        singular_values.append(float(string))
    ifile.close()

    dimension = len(singular_values)
    LambdaInverse = np.zeros((dimension, dimension))
    for i in range(dimension):
        LambdaInverse[i,i] = 1.0/singular_values[i]

    alpha = []
    ifile = open("alpha.txt", "r")
    for i, string in enumerate(ifile):
        alpha.append(float(string))
    ifile.close()

    D = []
    if (not os.path.exists("model.txt")):
        sys.exit("model.txt does not exist. ")
    ifile = open("model.txt", "r")  
    for index, string in enumerate(ifile):
        D.append(float(string.split()[1]))
    ifile.close()
    D = np.asarray(D)
    Nomega = len(D)

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
    
    spectrals = []
    probability = []
    chi_values = []
    entropy_values = []
    for i in range(len(alpha)):
        print alpha[i]
        fileName = "A_updated_alpha_" + str(alpha[i]) + ".txt"
        omega, A = readA(fileName)
        spectrals.append(np.asarray(A))
        chi_values.append(chi.chi(G_real_rotated, G_imag_rotated, K_real_rotated, K_imag_rotated, A, omega, LambdaInverse))
        entropy_values.append(entropy.entropy(omega, A, D))
        probability.append(np.exp(alpha[i]*entropy_values[i] - 0.5*chi_values[i]))

    spectral_mean = np.zeros(len(spectrals[0]))
    for i in range(1, len(spectrals)):
        spectral_mean[:] = spectral_mean[:] + ((-alpha[i] + alpha[i-1])*probability[i]*spectrals[i])[:]
    s = 0.0
    for i in range(1, len(alpha)):
        s = s + (-alpha[i] + alpha[i-1])*probability[i]
    print "s = ", s
    ofile = open("Chi_log_alpha.txt", "w")
    for i in range(len(alpha)):
        ofile.write(str(np.log(alpha[i])) + "   " + str(chi_values[i]) + "\n")
    ofile.close()
    ofile = open("entropy_log_alpha.txt", "w")
    for i in range(len(alpha)):
        ofile.write(str(np.log(alpha[i])) + "    " + str(entropy_values[i]) + "\n")
    ofile.close()
    ofile = open("P_log_alpha.txt", "w")
    for i in range(len(alpha)):
        ofile.write(str(np.log(alpha[i])) + "   " + str(probability[i]) + "\n")
    ofile.close()
    ofile = open("bryan.txt", "w")
    for i in range(len(omega)):
        ofile.write(str(omega[i]) + "    " + str(spectral_mean[i]/s) + "\n")
    ofile.close()
    return 0

main()
