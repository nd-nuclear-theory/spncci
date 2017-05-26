"""runsu3rme.py

  Example run script to calculate RMEs of relative unit tensors and
  the symplectic operators (A/B/Nintr) in SU(3)-NCSM many-body basis.

  Language: Python 3

  A. E. McCoy
  Department of Physics, University of Notre Dame

  - 12/1/16 (aem): Created.
  - 5/26/17 (mac):
      + Rename from generate_lsu3shell_rmes.py to runsu3rme01.py.
      + Move task handler out to spncci.py.

"""

import glob
import os
import sys

import mcscript
import spncci

mcscript.init()

##################################################################
# directory configuration
##################################################################

interaction_directory = os.environ["SPNCCI_INTERACTION_DIR"]
unit_tensor_directory = os.environ["SPNCCI_LSU3SHELL_DIR"]
interaction_filename_template = os.path.join(interaction_directory,"JISP16_Nmax20","JISP16_Nmax20_hw{:2.1f}_rel.dat")
unit_tensor_filename_template = os.path.join("lsu3shell_Z{nuclide[0]:02d}_N{nuclide[1]:02d}_{Nsigma_max:02d}")  

##################################################################
# build task list
##################################################################

task_list = [
    {
        "nuclide" : (3,3),
        "Nmax" : Nmax,
        "Nstep" : 2,
        "N1v" : 1,
        "Nsigma_0" : 11,
        "Nsigma_max" : Nmax,
        "J0" : -1,  # -1 for no restriction (needed for spncci); 0 for only Hamiltonian like operators
        "unit_tensor_filename_template" : unit_tensor_filename_template
    }
    for Nmax in mcscript.utils.value_range(0,6,2)
]

################################################################
# run control
################################################################

def task_descriptor(task):
    """"""
    return ("Z{nuclide[0]:d}-N{nuclide[1]:d}-Nmax{Nmax:02d}".format(**task))

def task_pool(task):
    """"""
    return ("Nmax{Nmax:02d}".format(**task))

##################################################################
# task control
##################################################################

mcscript.task.init(
    task_list,
    task_descriptor=task_descriptor,
    task_pool=task_pool,
    phase_handler_list=[
        spncci.do_generate_lsu3shell_rmes
        ],
    # Note: change to mcscript.task.archive_handler_hsi for tape backup
    archive_phase_handler_list=[mcscript.task.archive_handler_generic]
    )

################################################################
# termination
################################################################

mcscript.termination()
