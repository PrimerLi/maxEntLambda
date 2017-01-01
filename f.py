def f(alpha, G_real_rotated, G_imag_rotated, K_real_rotated, K_imag_rotated, A, D, omega, LambdaInverse):
    import numpy as np

    domega = omega[1] - omega[0]
    Nomega = len(omega)
    result = np.zeros(Nomega)
    for nw in range(Nomega):
        result[nw] = -alpha*domega*(1 + np.log(A[nw]/D[nw]))

    vector = G_real_rotated - K_real_rotated.dot(A)
    result = result + np.transpose(K_real_rotated).dot(LambdaInverse).dot(vector)

    vector = G_imag_rotated - K_imag_rotated.dot(A)
    result = result + np.transpose(K_imag_rotated).dot(LambdaInverse).dot(vector)

    return result
