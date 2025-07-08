import qforte as qf
import numpy as np
import os

# sys_str = 'n2'

# sys_str = 'h8'
# sys_str = 'h2be'
# sys_str = 'h2o'
sys_str = 'c6h6'

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
r_g_opt_thresh = 1.0e-6
pool_type = 'SD'
noise_factor = 0.0e-6  # Set to zero for exact residuals


# ===> VQE <===
timer.reset()

alg_pqe = qf.UCCNVQE(
    mol,
    computer_type = 'fci'
    )

alg_pqe.run(
    opt_thresh=r_g_opt_thresh, 
    pool_type=pool_type,
    noise_factor=noise_factor,
)

timer.record("VQE FCI")

print(timer)
print(f'\n\n Efci:   {mol.fci_energy:+12.10f}')