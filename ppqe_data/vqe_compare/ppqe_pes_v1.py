import qforte as qf
import numpy as np
import os

# sys_str = 'n2'
# sys_str = 'c2'
sys_str = 'o2'

# sys_str = 'h8'
# sys_str = 'h2be'
# sys_str = 'h2o'
# sys_str = 'c6h6'

tord = 2

# update_type = 'jacobi_like'
update_type = 'two_level_rotation'

def get_geom(sys_str, r):
    if(sys_str == 'c2'):
        return [
            ('C', (0., 0.,  0.0)),
            ('C', (0., 0.,  r))
        ]
    
    if(sys_str == 'o2'):
        return [
            ('O', (0., 0.,  0.0)),
            ('O', (0., 0.,  r))
        ]
    else:
        raise ValueError("Unknown system string")



timer = qf.local_timer()


# ===> Parameters <===
noise_factor = 0.0e-6  # Set to zero for exact residuals

e_opt_thresh = 1.0e-8

r_g_opt_thresh = 1.0e-6
pool_type = 'SD'



dt = 0.05
use_dt_from_l1_norm = False 

max_diis = 12
opt_maxiter = 50




# ===> Parameter Strings <===
# dt_str = f"{dt:1.3e}"
ndiis_str = f'{max_diis}'
etol_str = f'{e_opt_thresh:1.3e}'
nois_str = f'{noise_factor:1.3e}'
tord_str = str(tord) 
if(update_type == 'jacobi_like'):   
    update_str = 'jl'
elif(update_type == 'two_level_rotation'):
    update_str = 'tlr'


use_dt_from_l1_norm_str = str(use_dt_from_l1_norm).lower()

# ===> PPQE <===


dt_vals_used = []

Npts = 20

r_vec = np.linspace(0.4, 1.9, Npts)

r_vec = [1.1]

for i, r in enumerate(r_vec):

    geom = get_geom(sys_str, r)  

    timer.reset()

    mol = qf.system_factory(
        build_type='psi4', 
        symmetry='d2h',
        mol_geometry=geom, 
        basis='sto-6g',
        run_fci=1)

    timer.record("Setup")

    alg_ppqe = qf.UCCNPPQE(
    mol,
    computer_type = 'fci',
    diis_max_dim= max_diis,
    print_summary_file = True,
    )

    timer.record("Setup UCCNPPQE")
    
    timer.reset()

    alg_ppqe.run(
        pool_type=pool_type,
        opt_thresh = r_g_opt_thresh,
        opt_e_thresh = e_opt_thresh,
        opt_maxiter = opt_maxiter,
        noise_factor = noise_factor,
        time_step = dt,
        use_dt_from_l1_norm = use_dt_from_l1_norm,
        optimizer = 'rotation',
        ppqe_trotter_order = tord,
        update_type = update_type,  
        )
    
    dt_used = alg_ppqe._dt
    
    dt_str = f"{dt_used:1.3e}"
    r_str = f"{r:1.3e}"

    new_summary_name = f"ppqe_{sys_str}_r_{r_str}_dt_{dt_str}_inv_hnrm_{use_dt_from_l1_norm_str}_pool_{pool_type}_etol_{etol_str}_ndiis_{ndiis_str}_updt_{update_str}_tord_{tord_str}_sig_{nois_str}.dat"

    timer.record("PPQE FCI")


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