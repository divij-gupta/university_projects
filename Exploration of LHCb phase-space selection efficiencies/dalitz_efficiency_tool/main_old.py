
# decay_plot_1_test, bin_mids_x_test, bin_mids_y_test = uniform_bins_plot(
#     masses_Ks_Pip_2, masses_Ks_Pim_2, [bin_mids_m12[12:26], bin_mids_m13[12:26]], 'flat1_test')
# params_test, param_errs_test = fit_lin_plane(bin_mids_x_test, bin_mids_y_test, decay_plot_1_test, 50,
#                                              'flat1_test')
# print()
# print(params_test)
# print(param_errs_test)
#
# print()

# masses_Ks_Pip_1_quad_rej, masses_Ks_Pim_1_quad_rej = artificial_rejection_quad(masses_Ks_Pip_1, masses_Ks_Pim_1)
# decay_plot_1_quad_rej, _, _ = uniform_bins_plot(
#     masses_Ks_Pip_1_quad_rej, masses_Ks_Pim_1_quad_rej, [bin_mids_m12, bin_mids_m13], 'flat1_quad_rej', D0_decay_names)
# inv_mass_projection(masses_Ks_Pip_1_quad_rej, f'{D0_decay_names[1]} {D0_decay_names[2]}', 50, 'Ks_Pip_flat1_quad_rej')
# inv_mass_projection(masses_Ks_Pim_1_quad_rej, f'{D0_decay_names[1]} {D0_decay_names[3]}', 50, 'Ks_Pim_flat1_quad_rej')
#
# eff_1_quad_rej, eff_errs_1_quad_rej = efficiency_plot(decay_plot_1, decay_plot_1_quad_rej, bin_mids_m12, bin_mids_m13,
#                                                       50, 'flat1_true_quad_rej', D0_decay_names)
# eff_1_params_quad_rej_lin, eff_1_param_errs_quad_rej_lin = (
#     fit_lin_plane(bin_mids_m12, bin_mids_m13, eff_1_quad_rej, eff_errs_1_quad_rej, 50, 'eff_1_quad_rej_lin',
#                   D0_decay_names, 'efficiency'))
# print(eff_1_params_quad_rej_lin)
# print(eff_1_param_errs_quad_rej_lin)
# eff_1_params_quad_rej_quad, eff_1_param_errs_quad_rej_quad = (
#     fit_quad_plane(bin_mids_m12, bin_mids_m13, eff_1_quad_rej, eff_errs_1_quad_rej, 50, 'eff_1_quad_rej_quad',
#                    D0_decay_names, 'efficiency'))
# print(eff_1_params_quad_rej_quad)
# print(eff_1_param_errs_quad_rej_quad)

#-----------------------------------------------------------------------------------------------------------------------
