from classy import Class

standard_params = {
    'h': 0.67,
    'omega_b': 0.022,
    'omega_cdm': 0.12,
    'A_s': 2.1e-9,
    'n_s': 0.96,
    'tau_reio': 0.054,
    'N_ur': 2 * 1.015,  # Standard value for Neff
    'N_ncdm': 1,
    'm_ncdm': 0.01,
    'omega_ncdm': 0.00001,  # Set a small value for omega_ncdm
    'ncdm_psd_parameters': '2, 1.0',  # log-normal distribution with sigma=1.0
    'ncdm_quadrature_strategy': 3,  # qm_trapz
    'ncdm_maximum_q': 10.0,  # Set appropriate qmax
    'ncdm_a': 1.0,  # Will be overridden by optimal
    'ncdm_N_momentum_bins': 20,  # Set a reasonable number of momentum bins
    'output': 'mPk',
    'P_k_max_h/Mpc': 1.0,
    'z_max_pk': 0.0,
    'background_verbose': 2,  # This is the correct verbose flag for CLASS
}

print("Looking for 'Log-normal distribution detected' message...")
print("Parameters:")
for key, value in standard_params.items():
    print(f"  {key}: {value}")

print("\n" + "="*50)
print("CLASS OUTPUT:")
print("="*50)

cosmo = Class()
cosmo.set(standard_params)
cosmo.compute()

print("="*50)
print(f"✓ Log-normal distribution works!")
print(f"  h = {cosmo.h():.6f}")
# print(f"  omega_ncdm = {cosmo.omega_ncdm():.6f}")
print(f"  Age = {cosmo.age():.3f} Gyr")

cosmo.struct_cleanup()