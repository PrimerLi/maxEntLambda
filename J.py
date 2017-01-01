def J(alpha, A, omega, K_real_rotated, K_imag_rotated, LambdaInverse):
    import numpy as np
    
    Nomega = len(omega)
    result = np.zeros((Nomega, Nomega))
    domega = omega[1] - omega[0]
    
    result = -np.transpose(K_real_rotated).dot(LambdaInverse).dot(K_real_rotated)
    result = result -np.transpose(K_imag_rotated).dot(LambdaInverse).dot(K_imag_rotated)
    for i in range(Nomega):
        result[i,i] = result[i,i] - alpha*domega/A[i]
    return result
