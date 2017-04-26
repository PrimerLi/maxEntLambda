# maxEntLambda
max entropy method in python, covariance matrix is diagonalized and the Green function and kernel matrix are rotated into the eigen-space.
A deteailed description of the program is in the file algorithm.pdf. 

Guidelines about how to run the program. 
1. What does the programs accomplish? 

    The program takes as input the Matsubara imaginary frequency Green's function and the covariance matrix, and output the spectral function of the corresponding Green's function. 
  
2. How to run the program?

    First, make sure the input files are present. The program needs three input files. File 1 should contain the real part of the Matsubara Green's function; file 2 should contain the imaginary part of the Matsubara Green's function, and file 3 should contain the real part of the covariance matrix. The imaginary part of the covariance matrix is not needed, as explained the algorithm.pdf. 
    
    Second, run the command "python readmc.py" in your terminal. This step reads in the three input files, and generates the files that will be used in step 3. 
    
    Third, run the command "python main.py (optional parameters)" in your terminal. You can specifiy how to anneal the result through the optional parameters, or you can use the default values set in the program. The final spectral function is contained in the file A_initial.txt. You can delete all the other intermediate txt files that are generated when the program is running. 
