import qforte as qf
import numpy as np
import os

# sys_str = 'h8'
# sys_str = 'h2be'
# sys_str = 'h2o'
# sys_str = 'c6h6'

sys_str = 'o2'

if(sys_str == 'h8'):
    geom = [
        ('H', (0., 0., 1.00)),
        ('H', (0., 0., 2.00)),
        ('H', (0., 0., 3.00)),
        ('H', (0., 0., 4.00)),
        ('H', (0., 0., 5.00)),
        ('H', (0., 0., 6.00)),
        ('H', (0., 0., 7.00)),
        ('H', (0., 0., 8.00)),
        # ('H', (0., 0., 9.00)),
        # ('H', (0., 0., 10.00))
        ]

elif(sys_str == 'h2be'):
    geom = [
        ('Be', (0., 0., 0.00)), 
        ('H', (0., 0., -1.00)),
        ('H', (0., 0., +1.00)),
        ]   

elif(sys_str == 'h2o'):
    geom = [
        ('O', (0., 0., 0.00)), 
        ('H', (0., 0., -1.00)),
        ('H', (0., 0., +1.00)),
        ]   

elif(sys_str == 'o2'):
    geom = [
        ('O', (0., 0., -0.50)), 
        ('O', (0., 0., +0.50)),
        ]  

elif(sys_str == 'n2'):
    geom = [
        ('N', (0., 0., -1.00)), 
        ('N', (0., 0., +1.00)),
        ]  

elif(sys_str == 'c6h6'):
    geom = [
        ('C',   ( 0.000000,    1.396792,    0.000000)), 
        ('C',   ( 0.000000,   -1.396792,    0.000000)),
        ('C',   ( 1.209657,    0.698396,    0.000000)),
        ('C',   (-1.209657,   -0.698396,    0.000000)),
        ('C',   (-1.209657,    0.698396,    0.000000)),
        ('C',   ( 1.209657,   -0.698396,    0.000000)),
        ('H',   ( 0.000000,    2.484212,    0.000000)),
        ('H',   ( 2.151390,    1.242106,    0.000000)),
        ('H',   (-2.151390,   -1.242106,    0.000000)),
        ('H',   (-2.151390,    1.242106,    0.000000)),
        ('H',   ( 2.151390,   -1.242106,    0.000000)),
        ('H',   ( 0.000000,   -2.484212,    0.000000)),
                ]       


timer = qf.local_timer()

timer.reset()

if(sys_str == 'c6h6'):
    symm_str = 'c1'
    fdocc = 0
    fuocc = 0

    basis_set = 'cc-pvdz'
    avas_atoms_and_atomic_orbs = ['C 2pz']

    mol = qf.system_factory(
        build_type='pyscf', 
        symmetry=symm_str,
        mol_geometry=geom, 
        basis=basis_set, 
        run_fci=True, 
        use_avas=True, #                     <=====
        avas_atoms_or_orbitals=avas_atoms_and_atomic_orbs,
        run_ccsd=False,
        store_mo_ints=True,
        build_df_ham=False,
        num_frozen_uocc = fuocc,
        num_frozen_docc = fdocc,
        build_qb_ham = True,
        )
else:

    mol = qf.system_factory(
        build_type='psi4', 
        symmetry='d2h',
        mol_geometry=geom, 
        basis='sto-6g',
        run_fci=1)

timer.record("Setup")


# ===> Parameters <===

r_g_opt_thresh = 1.0e-9
pool_type = 'SD'
# e_opt_thresh = 1.0e-12


noise_list = [1.0e-6, 1.0e-5, 1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1]
# noise_list = [1.0e-8]

# Number of repeats
R = 10

max_diis = 12
opt_maxiter = 50

# ===> Parameter Strings <===
ndiis_str = f'{max_diis}'
gtol_str = f'{r_g_opt_thresh:1.3e}'


# ===> VQE <===
timer.reset()

alg_vqe = qf.UCCNVQE(
    mol,
    computer_type = 'fci',
    diis_max_dim= max_diis,
    print_summary_file = True,
    )


for j, noise_factor in enumerate(noise_list):
    print(f"\n\nRunning with noise_factor = {noise_factor:1.3e}\n\n")
    noise_str = f'{noise_factor:1.3e}'    

    for r in range(R):
        print(f"Repeat {r+1} of {R}")
        r_str = str(r).lower()
        
        try:
            alg_vqe.run(
                pool_type=pool_type,
                opt_thresh = r_g_opt_thresh,
                opt_maxiter = opt_maxiter,
                noise_factor = noise_factor,
                optimizer = 'BFGS',
                )
        
        except Exception as e:
            print(f"VQE run failed: {e}")
            continue  # Skip to next repeat
        

        new_summary_name = f"vqe_{sys_str}_pool_{pool_type}_gtol_{gtol_str}_ndiis_{ndiis_str}_updt_BFGS_noise_{noise_str}_r_{r_str}.dat"

        timer.record("VQE FCI")

        runs_data_dir = "vqe_noisy_runs_data"
        # Ensure the directory exists
        os.makedirs(runs_data_dir, exist_ok=True)
        new_summary_path = os.path.join(runs_data_dir, new_summary_name)

        if os.path.exists("summary.dat"):
            os.rename("summary.dat", new_summary_path)
            print(f"\n\nRenamed summary.dat to {new_summary_path}\n\n")

            # Prepend a line to the new summary file
            header_line = f'#Efci:{mol.fci_energy:+12.10f}'
            with open(new_summary_path, "r") as original:
                data = original.read()
            with open(new_summary_path, "w") as modified:
                modified.write(header_line + data)

        print(timer)

        print(f'\n\n Efci:   {mol.fci_energy:+12.10f}')