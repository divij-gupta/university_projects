import numpy as np

def gaussian_func(x, A, mu, sig):

    norm = A / (sig * np.sqrt(2*np.pi))
    exponent = -0.5 * ((x - mu) / sig)**2

    return norm * np.exp(exponent)

def chebyshev_1D(value, degree):

    if degree == 0:
        return np.ones(value.shape)

    if degree == 1:
        return value

    return 2 * value * chebyshev_1D(value, degree-1) - chebyshev_1D(value, degree-2)

def find_pairs(degree):

    pairs = []
    i, j = 0, 0

    while i <= degree:
        while j <= degree:
            current_sum = i + j
            if current_sum <= degree:
                pairs.append([j, i])
            j += 1
        i += 1
        j = 0

    # ordering for degree==3 for example is : [[0, 0], [1, 0], [2, 0], [3, 0], [0, 1], [1, 1], [2, 1], [0, 2], [1, 2], [0, 3]]
    # (this determines the order of best fit coefficients in fitted plane - doesn't affect fit result!)

    return pairs
