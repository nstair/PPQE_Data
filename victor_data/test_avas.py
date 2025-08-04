import qforte as qf

# cyclobutadiene
# geom = [
#     ('C', ( 0.0000,  1.2070,  0.0000)),
#     ('C', ( 1.2070,  0.0000,  0.0000)),
#     ('C', ( 0.0000, -1.2070,  0.0000)),
#     ('C', (-1.2070,  0.0000,  0.0000)),
#     ('H', ( 0.0000,  2.1470,  0.0000)),
#     ('H', ( 2.1470,  0.0000,  0.0000)),
#     ('H', ( 0.0000, -2.1470,  0.0000)),
#     ('H', (-2.1470,  0.0000,  0.0000)),
# ]

# butadiene
geom = [
    ('C', (-1.2090,  0.0000,  0.0000)),
    ('C', ( 0.0000,  0.0000,  0.0000)),
    ('C', ( 1.3390,  0.0000,  0.0000)),
    ('C', ( 2.5480,  0.0000,  0.0000)),
    ('H', (-1.7490,  0.9430,  0.0000)),
    ('H', (-1.7490, -0.9430,  0.0000)),
    ('H', ( 0.5280,  0.9430,  0.0000)),
    ('H', ( 0.5280, -0.9430,  0.0000)),
    ('H', ( 3.0880,  0.9430,  0.0000)),
    ('H', ( 3.0880, -0.9430,  0.0000)),
]

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
