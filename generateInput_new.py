def gauss(x, mu, sigma):
    import numpy as np
    return 1.0/np.sqrt(2*np.pi*sigma**2)*np.exp(-(x-mu)**2/(2*sigma**2))

def main():
    import os
    import sys
    import numpy as np
    import numpy.linalg 

    ifile = open("Niom", "r")
    string = ifile.readline()
    Niom = int(string)
    ifile.close()

    ifile = open("param_CDMFT", "r")
    string = ifile.readline()
    a = string.split(",")
    beta = float(a[0])
    string = ifile.readline()
    string = string.replace(",", " ")
    a = string.split()
    Nbin = int(a[2])
    ifile.close()

    xiom = np.zeros(Niom)
    Giom_cc_real = []
    Giom_cc_imag = []
    Giom_ff_real = []
    Giom_ff_imag = []
    ifile = open("grtot_iom", "r")
    for i in range(Nbin):
        string = ifile.readline()
        for row in range(2):
            for col in range(2):
                string = ifile.readline()
                for nw in range(Niom):
                    string = ifile.readline()
                    string = string.replace("(", "")
                    string = string.replace(")", "")
                    string = string.replace(",", "    ")
                    a = string.split()
                    frequency = float(a[0])
                    real = float(a[1])
                    imag = float(a[2])
                    xiom[nw] = frequency
                    if (row == 0 and col == 0):
                        Giom_cc_real.append(real)
                        Giom_cc_imag.append(-imag)
                    if (row == 1 and col == 1):
                        Giom_ff_real.append(real)
                        Giom_ff_imag.append(-imag)
    ifile.close()

    def covariance(listA, listB):
        import sys
        import numpy as np
        if (len(listA) != len(listB)):
            print "Dimension incompatible in covariance. "
            sys.exit(-1)

        length = len(listA)
        czero = 0.0*1j
        def mean(array):
            s = 0.0*1j
            for i in range(len(array)):
                s = s + array[i]
            return s/float(len(array))
        meanA = mean(listA)
        meanB = mean(listB)
        s = czero
        for i in range(len(listA)):
            s = s + np.conj(listA[i] - meanA)*(listB[i] - meanB)
        return s/float(length*(length-1))

    def mean(array):
        s = 0
        length = len(array)
        for i in range(length):
            s = s + array[i]
        return s/float(length)

    def printFile(x, y, fileName):
        ofile = open(fileName, "w")
        for i in range(len(x)):
            ofile.write(str(x[i]) + "    " + str(y[i]) + "\n")
        ofile.close()

    Giom_cc_real_mean = np.zeros(Niom) 
    Giom_cc_imag_mean = np.zeros(Niom)
    Giom_ff_real_mean = np.zeros(Niom)
    Giom_ff_imag_mean = np.zeros(Niom)
    for nw in range(Niom):
        tempListReal = []
        tempListImag = []
        for i in range(Nbin):
            tempListReal.append(Giom_cc_real[nw + i*Niom])
            tempListImag.append(Giom_cc_imag[nw + i*Niom])
        Giom_cc_real_mean[nw] = mean(tempListReal)
        Giom_cc_imag_mean[nw] = mean(tempListImag)
        tempListReal = []
        tempListImag = []
        for i in range(Nbin):
            tempListReal.append(Giom_ff_real[nw + i*Niom])
            tempListImag.append(Giom_ff_imag[nw + i*Niom])
        Giom_ff_real_mean[nw] = mean(tempListReal)
        Giom_ff_imag_mean[nw] = mean(tempListImag)
    printFile(xiom, Giom_cc_real_mean, "G_cc_real.txt")
    printFile(xiom, Giom_cc_imag_mean, "G_cc_imag.txt")
    printFile(xiom, Giom_ff_real_mean, "G_ff_real.txt")
    printFile(xiom, Giom_ff_imag_mean, "G_ff_imag.txt")
            
    covariance_matrix_real = np.zeros((Niom, Niom))
    covariance_matrix_imag = np.zeros((Niom, Niom))
    for nw in range(Niom):
        for nwp in range(Niom):
            tempListA = []
            tempListB = []
            for i in range(Nbin):
                tempListA.append(Giom_cc_real[nw + i*Niom] + Giom_cc_imag[nw+i*Niom]*1j)
                tempListB.append(Giom_cc_real[nwp + i*Niom] + Giom_cc_imag[nwp+i*Niom]*1j)
            temp = covariance(tempListA, tempListB)
            covariance_matrix_real[nw, nwp] = temp.real
            covariance_matrix_imag[nw, nwp] = temp.imag            

    covariance_matrix_real_ff = np.zeros((Niom, Niom))
    covariance_matrix_imag_ff = np.zeros((Niom, Niom))
    for nw in range(Niom):
        for nwp in range(Niom):
            tempListA = []
            tempListB = []
            for i in range(Nbin):
                tempListA.append(Giom_cc_real[nw + i*Niom] + Giom_ff_imag[nw+i*Niom]*1j)
                tempListB.append(Giom_cc_real[nwp + i*Niom] + Giom_ff_imag[nwp+i*Niom]*1j)
            temp = covariance(tempListA, tempListB)
            covariance_matrix_real_ff[nw, nwp] = temp.real
            covariance_matrix_imag_ff[nw, nwp] = temp.imag
                
    nl = Niom
    Nuse = Niom
    run = 2
    isign = 1
    ofile = open("inputData.txt", "w")
    ofile.write("Results nl\n")
    ofile.write(str(nl) + "  " + str(Nuse) + "  " + str(run) + "  " + str(beta) + "  " + str(isign) + "\n")
    ofile.write("Giom(n)\n")
    for i in range(Niom):
        ofile.write(str(xiom[i]) + "  " + str(Giom_cc_real_mean[i]) + "\n")
    for i in range(Niom):
        ofile.write(str(xiom[i]) + "  " + str(Giom_cc_imag_mean[i]) + "\n")
    ofile.write("Gij-GiGj\n")
    for i in range(Niom):
        for j in range(i, Niom):
            ofile.write(str(i+1) + "   " + str(j+1) + "   " + str(covariance_matrix_real[i,j]) + "\n")
    ofile.close()

    nl = Niom
    Nuse = Niom
    run = 2
    isign = 1
    ofile = open("inputData_ff.txt", "w")
    ofile.write("Results nl\n")
    ofile.write(str(nl) + "  " + str(Nuse) + "  " + str(run) + "  " + str(beta) + "  " + str(isign) + "\n")
    ofile.write("Giom(n)\n")
    for i in range(Niom):
        ofile.write(str(xiom[i]) + "  " + str(Giom_ff_real_mean[i]) + "\n")
    for i in range(Niom):
        ofile.write(str(xiom[i]) + "  " + str(Giom_ff_imag_mean[i]) + "\n")
    ofile.write("Gij-GiGj\n")
    for i in range(Niom):
        for j in range(i, Niom):
            ofile.write(str(i+1) + "   " + str(j+1) + "   " + str(covariance_matrix_real_ff[i,j]) + "\n")
    ofile.close()

    def printMatrix(matrix, output):
        ofile = open(output, "w")
        shape = matrix.shape
        row = shape[0]
        col = shape[1]
        for i in range(row):
            for j in range(col):
                ofile.write(str(i+1) + "   " + str(j+1) + "  " + str(matrix[i,j]) + "\n")
        ofile.close()

    printMatrix(covariance_matrix_real, "CM_cc_real.txt")      
    printMatrix(covariance_matrix_imag, "CM_cc_imag.txt")
    printMatrix(covariance_matrix_real_ff, "CM_ff_real.txt")
    printMatrix(covariance_matrix_imag_ff, "CM_ff_imag.txt")


    if (not os.path.exists("input")):
        os.mkdir("input")
    os.system("cp G_cc_real.txt G_cc_imag.txt CM_cc_real.txt CM_cc_imag.txt ./input")
    os.system("cp G_ff_real.txt G_ff_imag.txt CM_ff_real.txt CM_ff_imag.txt ./input")

    return 0

main()
