# LMPTorch

Reference: Gang Seob Jung, Hunjoo Myung, and Stephan Irle. "Artificial Neural Network Potentials for Mechanics and Fracture Dynamics of Two-Dimensional Crystals." Machine Learning: Science and Technology (2023)

GS JUNG@ORNL made the codes based on fix_python.* and fix_dftb.* (Reference: Gang Seob Jung, Stephan Irle, and Bobby Sumpter. "Dynamic aspects of graphene deformation and fracture from approximate density functional theory." Carbon (2022)): https://github.com/gsjung0419/DFTBP

HJ MYUNG@KISTI helped to make the code work with a python model pointer. 

0. Intro

 -. The codes allow LAMMPS communicate with python codes (e.g., Pytorch based NNP;torchani and ASE library for ab-initio calculations, QE, NWchem, VASP and other NNPs, e.g., SchNet)

 -. The example of TorchANI is provided (https://github.com/gsjung0419/LammpsANI)

 -. Other examples will be provided soon in following studies

1. Requirements

 -. LAMMPS 29OCT20

 -. MPI (MPICH or OPENMPI)

 -. mpi4py

2. Installation

 -. copy the files (fix_python_torch.cpp and fix_python_torch.h) to lammps/src

 -. Edit lib/python/Makefile.lammps

   e.g., python_SYSINC = -I/home/$USER/anaconda3/include/python3.8
 
   e.g., python_SYSLIB = -L/home/$USER/anaconda3/lib -lpython3.8 -lnsl -ldl -lreadline -ltermcap -lpthread -lutil -lm

 -. make yes-python (in lammps/src)

 -. make mode=shared mpi (Please check whether a shared library of lammps is generated, e.g., liblammps.so)

 -. make install-python



3. LAMMPS 2025 status

This branch has been updated enough to compile with LAMMPS 22Jul2025.  The
current tested style name is:

```
fix python/torch
```

The old alias `fix python` is intentionally not registered here because recent
LAMMPS versions already provide their own built-in `fix python` style.

The 2025 port includes the following compatibility updates:

-. updated fix metadata for the newer LAMMPS API (`virial_global_flag`,
   `energy_global_flag`, `thermo_energy`)

-. added `extvector=1`, required by LAMMPS 2025 when `vector_flag` is set

-. replaced the old Python pointer/evaluation helpers with the current LAMMPS
   Python utility wrappers

-. keeps a persistent Python `Model` object and passes it to the callback as the
   third argument

4. Smoke-tested environment

The following local stack was tested:

```
LAMMPS: 22 Jul 2025
Packages: PYTHON, COLVARS, DFTBP
TorchANI: gsjung0419/torchani_gs branch v, version 2.2.0+gs0
Example: GS_LMPANI Benchmark_Size/10nm_10nm.data
Run: run 1
MPI: 1 rank and 4 ranks
```

Example command used for the 4-rank smoke test:

```
PYTHONPATH=/path/to/GS_LMPANI \
KMP_DUPLICATE_LIB_OK=TRUE \
OMP_NUM_THREADS=1 \
mpirun -np 4 lmp -in input.in
```

The test completed with LAMMPS using a `2 by 2 by 1 MPI processor grid` for
3680 atoms.  This verifies that the LAMMPS fix, embedded Python callback,
mpi4py gather, and TorchANI model path can proceed for a short run.

5. TorchANI compatibility note

The current GS_LMPANI example uses legacy TorchANI APIs, including:

```
torchani.neurochem.Constants
torchani.AEVComputer(**consts)
ANIModel([C_network])
nn._atomic_energies(...)
```

Official TorchANI 2.7/2.8 no longer provides the same legacy API surface.  A
future port to official TorchANI 2.8+ will likely need coordinated updates in
the Python model and force callback, not only this LAMMPS fix.  In particular,
the newer API uses `load_aev_constants_and_symbols`,
`AEVComputer.from_constants`, symbol-mapped network containers, and atomic
energy output through the model forward path.

For now, use `gsjung0419/torchani_gs` branch `v` when reproducing the legacy
GS_LMPANI workflow:

```
pip install "git+ssh://git@github.com/gsjung0419/torchani_gs.git@v"
```

Please report any bug/comments to jungg@ornl.gov or gs4phone@gmail.com
