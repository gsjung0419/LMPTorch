/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://lammps.sandia.gov/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Richard Berger (Temple U)
------------------------------------------------------------------------- */

#include "fix_python_torch.h"

#include "error.h"
#include "lmppython.h"
#include "python_compat.h"
#include "python_utils.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "memory.h"
#include "error.h"
#include "utils.h"

#include <iostream>
#include <cstring>
#include <Python.h>   // IWYU pragma: export

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixPythonTorch::FixPythonTorch(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  //Initiation of local vector
  local_flag=1;
  virial_global_flag=1;
  thermo_virial=1;
  energy_global_flag=1;
  thermo_energy=1;
  vector_flag=1;
  extvector=1;
  size_vector=10;
  
  vector_local = new double[10];
  for(int i =0;i<10;i++){
    vector_local[i]=0.0;
  }
  
  if (narg != 6) error->all(FLERR,"Illegal fix python/torch command");

  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  if (nevery <= 0) error->all(FLERR,"Illegal fix python/torch command");


  // ensure Python interpreter is initialized
  python->init();

  if (strcmp(arg[4],"post_force") == 0) {
    selected_callback = POST_FORCE;
  } else if (strcmp(arg[4],"end_of_step") == 0) {
    selected_callback = END_OF_STEP;
  } else {
    error->all(FLERR,"Unsupported callback name for fix python/torch");
  }

  // get Python function
  PyUtils::GIL lock;

  PyObject * pyMain = PyImport_AddModule("__main__");

  if (!pyMain) {
    PyUtils::Print_Errors();
    error->all(FLERR,"Could not initialize embedded Python");
  }

  char * fname = arg[5];
  pFunc = PyObject_GetAttrString(pyMain, fname);

  if (!pFunc) {
    PyUtils::Print_Errors();
    error->all(FLERR,"Could not find Python function");
  }

  PyObject * module = PyImport_ImportModule("model");
  if (!module ){
	PyUtils::Print_Errors();
	error->all(FLERR,"Could not initialize embedded Python");
  }
  	pModel = PyObject_CallMethod((PyObject *)module,"Model", NULL);
  	if (!pModel) {
    	PyUtils::Print_Errors();
    	error->all(FLERR,"Could not find Python class, \"Model\"");
	}

  lmpPtr = PyCapsule_New((void *)lmp, nullptr, nullptr);
}

/* ---------------------------------------------------------------------- */

FixPythonTorch::~FixPythonTorch()
{
  PyUtils::GIL lock;
  Py_CLEAR(lmpPtr);
  Py_CLEAR(pFunc);
  Py_CLEAR(pModel);
  delete [] vector_local;
//std::cout<<"Torch Debug: FixPythonTorch destructor called: "<<std::endl;  
}

/* ---------------------------------------------------------------------- */

int FixPythonTorch::setmask()
{
  int mask = 0;
  //mask |= INITIAL_INTEGRATE;
  //mask |= FINAL_INTEGRATE;
  mask |= PRE_REVERSE;
  mask |= POST_FORCE;
  mask |= MIN_POST_FORCE;

  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPythonTorch::init()
{
  // error checks
  // for now, assume nlocal will never change
  
  //std::cout<<"Torch Debug: FixPythonTorch init() called: "<<std::endl;  
}

/* ---------------------------------------------------------------------- */

void FixPythonTorch::init_list(int /*id*/, NeighList * /*ptr*/)
{
  // list = ptr;
}


/* ---------------------------------------------------------------------- */

void FixPythonTorch::setup(int vflag)
{

  post_force(vflag);
  //std::cout<<"Torch Debug: FixPythonTorch setup() called: "<<std::endl;  
}

/* ---------------------------------------------------------------------- */

void FixPythonTorch::min_setup(int vflag)
{
  post_force(vflag);
  //std::cout<<"Torch Debug: FixPythonTorch min_setup() called: "<<std::endl;  
}

/* ---------------------------------------------------------------------- */

void FixPythonTorch::setup_pre_reverse(int eflag, int vflag)
{
  pre_reverse(eflag,vflag);
  //std::cout<<"Torch Debug: FixPythonTorch setup_pre_reverse() called: "<<std::endl;  
}

void FixPythonTorch::initial_integrate(int /*vflag*/) {}

/* ---------------------------------------------------------------------- */

void FixPythonTorch::pre_reverse(int eflag, int /*vflag*/)
{
  eflag_caller = eflag;
  //std::cout<<"Torch Debug: FixPythonTorch pre_reverse() called: "<<std::endl;    
}


void FixPythonTorch::post_force(int vflag)
{
  //if (update->ntimestep % nevery != 0) return;

  int eflag = eflag_caller;
  ev_init(eflag,vflag);
  
  PyUtils::GIL lock;

  //PyObject * arglist = Py_BuildValue("(Oi)", ptr, vflag);
  PyObject * arglist = Py_BuildValue("(OiO)", (PyObject *)lmpPtr, vflag, (PyObject *)pModel); // Model object pointer is the third parameter

  PyObject * result = PyObject_CallObject((PyObject*)pFunc, arglist);

  if (!result) {
    Py_DECREF(arglist);
    PyUtils::Print_Errors();
    error->all(FLERR,"Fix python/torch post_force() method failed");
  }

  Py_CLEAR(result);
  Py_DECREF(arglist);

  nnani_energy = vector_local[9];
  virial[0]=vector_local[0];
  virial[1]=vector_local[4];
  virial[2]=vector_local[8];    
  virial[3]=0.5*(vector_local[1]+vector_local[3]);
  virial[4]=0.5*(vector_local[2]+vector_local[6]);
  virial[5]=0.5*(vector_local[5]+vector_local[7]);   

  //return;
  
}

/* ---------------------------------------------------------------------- */
void FixPythonTorch::min_post_force(int vflag)
{
  //To perform energy minimization, eflag should be set > 0 
  eflag_caller=1;
  post_force(vflag);
  //std::cout<<"Torch Debug: FixPythonTorch min_post_force() called: "<<std::endl;
  //std::cout<<"eflag_callerL "<<eflag_caller<<std::endl;  
  
  
}

/* ---------------------------------------------------------------------- */

void FixPythonTorch::final_integrate() {}

/* ---------------------------------------------------------------------- */


void FixPythonTorch::reset_dt()
{

}

/* ---------------------------------------------------------------------- */

double FixPythonTorch::compute_scalar()
{
  return nnani_energy;
}

/* ---------------------------------------------------------------------- */

double FixPythonTorch::compute_vector(int n)
{
  return vector_local[n];
}

/* ---------------------------------------------------------------------- */

void FixPythonTorch::write_restart(FILE *fp){
#define RESTART_ITEM 6
  double buf[RESTART_ITEM+1];
  //if(comm->me == 0){
  //check MPI?
    buf[0]=nnani_energy;
    for (int i=0;i<6;i++){
      buf[i+1]=virial[i];
    } 
    //}
}

/* ---------------------------------------------------------------------- */

void FixPythonTorch::restart(char *buf){
  double *restored = (double *)buf;
  nnani_energy=restored[0];
  for (int i=0;i<6;i++){
    virial[i]=restored[i+1];
  }
  delete [] restored;
}

/* ---------------------------------------------------------------------- */
double FixPythonTorch::memory_usage(){
  //std::cout<<"Python Torch Debug: memory_usage() called: "<<std::endl;  
  
  double bytes = 0.0;
  bytes += 10 * sizeof(double);
  return bytes;
}
