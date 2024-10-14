import numpy as np
import scipy
import math
import scipy.special as sp

def pdf(s, b):

    return (1/math.factorial(b)) * np.power(s+b, b) * np.exp(-(s+b))

# def log_pdf(s, b):
#
#     return b * np.log(s + b) - (s + b)

# def find_alpha_other(s_up, b):
#
#     # Logarithm of the gamma function
#     log_gamma_b_plus_1 = np.log(sp.gamma(b + 1))
#
#     # Logarithm of the upper incomplete gamma function
#     log_upper_incomplete_gamma = np.log(sp.gammainc(b + 1, s_up + b) * sp.gamma(b + 1))
#
#     # Logarithm of the integral from 0 to infinity
#     log_integral_inf = log_gamma_b_plus_1
#
#     # Ratio of the two integrals in the logarithmic domain
#     log_ratio = log_upper_incomplete_gamma - log_integral_inf
#
#     # Exponentiate the result to get the ratio
#     ratio = np.exp(log_ratio)
#
#     return ratio

def integrate(upper, b, diff=0.00001):

    no_points = int(upper/diff) + 1

    x_i = np.linspace(0, upper, no_points)[:-1]
    f_i = pdf(x_i, b)

    diff = np.mean(np.diff(x_i))

    return np.sum(f_i*diff)

def find_alpha(s_up, b):

    numerator = scipy.integrate.quad(pdf, 0, s_up, args=(b))[0]
    denominator = scipy.integrate.quad(pdf, 0, np.inf, args=(b))[0]

    return numerator / denominator

# def find_alpha_manual(s_up, b):
#
#     numerator = integrate(s_up, b)
#     denominator = integrate(10, b)
#
#     return numerator / denominator

def optimise(s_up, b, alpha):

    return np.abs(find_alpha(s_up, b) - alpha)

b_test = 17
alpha_test = 0.9
s_up_guess = 5.5

result = scipy.optimize.minimize(optimise, np.array([s_up_guess]), args=(b_test, alpha_test), bounds=[(0, None)])
print(result)

print(find_alpha(12, b_test))

# print(optimise(result.x, b_test, alpha_test))

# print(find_alpha(10, 35))
# print(find_alpha_other(1100, 1565))
# print(find_alpha_manual(670/750, 1))
