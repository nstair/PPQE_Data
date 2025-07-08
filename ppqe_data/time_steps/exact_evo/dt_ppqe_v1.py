import qforte as qf
import numpy as np
import os

sys_str = 'h8'
# sys_str = 'h2be'
# sys_str = 'h2o'
# sys_str = 'n2'

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

elif(sys_str == 'n2'):
    geom = [
        ('N', (0., 0., -1.00)), 
        ('N', (0., 0., +1.00)),
        ]   


timer = qf.local_timer()

timer.reset()

mol = qf.system_factory(
    build_type='psi4', 
    symmetry='d2h',
    mol_geometry=geom, 
    basis='sto-6g',
    run_fci=1)

timer.record("Psi4 Setup")


# ===> Parameters <===
e_opt_thresh = 1.0e-8
r_g_opt_thresh = 1.0e-6
pool_type = 'SD'
noise_factor = 0.0e-4  # Set to zero for exact residuals
dt = 0.05
use_dt_from_l1_norm = False 

# update_type = 'jacobi_like'
update_type = 'two_level_rotation'
max_diis = 12
tord = np.inf
opt_maxiter = 100




# ===> Parameter Strings <===
dt_str = f"{dt:1.3e}"
ndiis_str = f'{max_diis}'
etol_str = f'{e_opt_thresh:1.3e}'
tord_str = str(tord) 
if(update_type == 'jacobi_like'):   
    update_str = 'jl'
elif(update_type == 'two_level_rotation'):
    update_str = 'tlr'

new_summary_name = f"ppqe_{sys_str}_dt_{dt_str}_pool_{pool_type}_etol_{etol_str}_ndiis_{ndiis_str}_updt_{update_str}.dat"


# ===> PPQE <===
timer.reset()

alg_ppqe = qf.UCCNPPQE(
    mol,
    computer_type = 'fci',
    diis_max_dim= max_diis,
    print_summary_file = True,
    )

alg_ppqe.run(
    pool_type=pool_type,
    opt_thresh = r_g_opt_thresh,
    opt_e_thresh = e_opt_thresh,
    opt_maxiter = opt_maxiter,
    noise_factor = noise_factor,
    time_step = dt,
    use_dt_from_l1_norm = use_dt_from_l1_norm,
    optimizer = 'rotation',
    ppqe_trotter_order = np.inf,
    update_type = update_type,  
    )

timer.record("PPQE FCI")

# ===> Summary <===
# if os.path.exists("summary.dat"):
#     os.rename("summary.dat", new_summary_name)
#     print(f"\n\nRenamed summary.dat to {new_summary_name}\n\n")

if os.path.exists("summary.dat"):
    os.rename("summary.dat", new_summary_name)
    print(f"\n\nRenamed summary.dat to {new_summary_name}\n\n")

    # Prepend a line to the new summary file
    header_line = f'#Efci:{mol.fci_energy:+12.10f}'
    with open(new_summary_name, "r") as original:
        data = original.read()
    with open(new_summary_name, "w") as modified:
        modified.write(header_line + data)

print(timer)
print(f'\n\n Efci:   {mol.fci_energy:+12.10f}')