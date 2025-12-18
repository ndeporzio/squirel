import numpy as np
import matplotlib.pyplot as plt
from classy import Class
x = None
sigma_array = [x]
# print(sigma_array != None)
print(sigma_array != None and len(sigma_array) > 0)
# # Initialize CLASS
# cosmo = Class()

# # Use default parameters
# cosmo.set({'output':'tCl,pCl,lCl,mPk',
#             'P_k_max_1/Mpc':10,
#             'l_max_scalars':2500,
#             'lensing': 'yes',
#             'write_background':'yes',            # Write background parameter table
#             'background_verbose': 3,
#             'thermodynamics_verbose': 1,
#             'perturbations_verbose': 1,
#             'input_verbose': 1,
#             'transfer_verbose': 1,
#             'primordial_verbose': 1,
#             'harmonic_verbose': 1,
#             'fourier_verbose': 1,
#             'lensing_verbose': 1,
#             'output_verbose': 1,
#         })

# nu_ur_array = [
#     3.044,
#     2.0308,
#     1.0176,
#     0.00441
# ]
# N_mnu = 1
# M_mnu = 0.06
# omega_m = 0.1431354439      
# m_ncdm_str = ','.join([str(M_mnu/N_mnu) for x in range(N_mnu)])
# sigma = 0.6
# cosmo.set({
#     'N_ncdm':N_mnu,
#     'N_ur':nu_ur_array[N_mnu],
#     'omega_m':omega_m,
#     'background_verbose': 0,
#     })
# if N_mnu>0:
#     cosmo.set({
#         'm_ncdm':m_ncdm_str,
#         'T_ncdm':0.7161790659871354,
#         'ncdm_quadrature_strategy':'0',
#         'ncdm_maximum_q':'20',
#         'ncdm_N_momentum_bins':'20',
#         'ncdm_a':'20',
#         'ncdm_psd_parameters':'3,'+str(sigma)+',0', #RD distribution (sigma ignored)
#     })
# cosmo.compute()
# print(cosmo.Neff()-3.044)

# # # print([str(M_mnu/N_mnu) for x in range(N_mnu)]) #remove brackets
# # # # same as above but returning '0.02,0.02,0.02' for N_mnu=3
# # # m_ncdm_str = ','.join([str(M_mnu/N_mnu) for x in range(N_mnu)])
# # # print(m_ncdm_str)
# # print(cosmo.Neff())
# # print(cosmo.Om_ncdm(0)*cosmo.h()**2)
# # if N_mnu>0:
# #     cosmo.set({
# #         'm_ncdm':m_ncdm_str,
# #         'ncdm_quadrature_strategy':'0',
# #         'ncdm_maximum_q':'20',
# #         'ncdm_N_momentum_bins':'20',
# #         'ncdm_a':'20',
# #         'ncdm_psd_parameters':'2,'+str(sigma)+',0',
# #     })
# # cosmo.compute()
# # print(cosmo.Neff())

# # #

# # sigma = 0.5
# # if sigma != None and not 0.24 <= sigma <= 0.77:
# #     print('true')
# # else:
# #     print('false')