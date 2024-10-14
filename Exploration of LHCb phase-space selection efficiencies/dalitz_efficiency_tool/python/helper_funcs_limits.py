import numpy as np

def within_kinematic_range(m12_sq, m13_sq, masses):

    M, m1, m2, m3 = masses

    E1 = 1 / (2 * np.sqrt(m12_sq)) * (m12_sq + m1 ** 2 - m2 ** 2)
    E3 = 1 / (2 * np.sqrt(m12_sq)) * (M ** 2 - m12_sq - m3 ** 2)

    m13_sq_min = (E1 + E3) ** 2 - np.power(np.sqrt(np.power(E1, 2) - m1 ** 2) + np.sqrt(np.power(E3, 2) - m3 ** 2), 2)
    m13_sq_max = (E1 + E3) ** 2 - np.power(np.sqrt(np.power(E1, 2) - m1 ** 2) - np.sqrt(np.power(E3, 2) - m3 ** 2), 2)

    return (m13_sq_min < m13_sq) & (m13_sq < m13_sq_max)

def kinematic_limits_mask(x, y, decay_masses, exclude_boundary=False):

    xx, yy = np.meshgrid(x, y)
    x_diff = np.mean(np.diff(x)) / 2
    y_diff = np.mean(np.diff(y)) / 2

    bin_vals_top_left = np.array([xx - x_diff, yy + y_diff])
    bin_vals_top_right = np.array([xx + x_diff, yy + y_diff])
    bin_vals_bottom_left = np.array([xx - x_diff, yy - y_diff])
    bin_vals_bottom_right = np.array([xx + x_diff, yy - y_diff])
    bin_vals_all = [bin_vals_top_left, bin_vals_top_right, bin_vals_bottom_left, bin_vals_bottom_right]

    masks_all = []
    for bin_vals in bin_vals_all:
        mask_temp = within_kinematic_range(bin_vals[0], bin_vals[1], decay_masses)
        masks_all.append(mask_temp)

    masks_all = np.dstack(masks_all)

    if exclude_boundary:
        return np.all(masks_all, axis=2)
    else:
        return np.any(masks_all, axis=2)
