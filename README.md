Reference: Jung, Gang Seob, Hoon Joo Myung, and Stephan Irle. "Artificial Neural Network Potentials for Mechanics and Fracture Dynamics of Two-Dimensional Crystals." Machine Learning: Science and Technology (2023)

GS JUNG@ORNL made the codes based on fix_dftb.cpp and fix_dftb.h (Reference: Jung, Gang Seob, Stephan Irle, and Bobby Sumpter. "Dynamic aspects of graphene deformation and fracture from approximate density functional theory." Carbon (2022)) Check https://github.com/gsjung0419/DFTBP
HJ MYUNG@KISTI helped to make the code work with a python model pointer. 

0. Intro

 -. The codes allow LAMMPS communicate with python codes (e.g., Pytorch based NNP;torchani and ASE library for ab-initio calculations, QE,NWchem,VASP and other NNP;SchNet)

 -. The example of TorchANI is provided (https://github.com/gsjung0419/LammpsANI)

 -. Other examples will be provided soon in following studies

1. Requirements

 -. LAMMPS 29OCT20

 -. MPI (MPICH or OPENMPI)

 -. mpi4py

2. Installation

 -. copy the files (fix_python_torch.cpp and fix_python_torch.h) to lammps/src

 -. Edit lib/python/Makefile.lammps

   e.g., python_SYSINC = -I/home/8gj/anaconda3/include/python3.8
 
   e.g., python_SYSLIB = -L/home/8gj/anaconda3/lib -lpython3.8 -lnsl -ldl -lreadline -ltermcap -lpthread -lutil -lm

 -. make yes-python (in lammps/src)

 -. make mode=shared mpi (Please check whether a shared library of lammps is generated, e.g., liblammps.so)

 -. make install-python



Please report any bug/commetns to jungg@ornl.gov or gs4phone@gmail.com
