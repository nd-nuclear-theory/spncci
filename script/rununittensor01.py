"""rununittensor01.py

  Example run script to generate relative unit tensors and symplectic
  operators (Arel/Brel/Nrel) in LSU3shell operator format.

  Language: Python 3

  A. E. McCoy and M. A. Caprio
  Department of Physics
  University of Notre Dame

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
# build task list
##################################################################

N1v_set = [0,1]
Nsigmamax_set = mcscript.utils.value_range(0,11,1)

task_list = [
    {
        "N1v" : N1v,
        "Nsigmamax" : Nsigmamax,
        "Nstep" : 2
    }
    for N1v in N1v_set
    for Nsigmamax in Nsigmamax_set
]

################################################################
# run control
################################################################

def task_descriptor(task):
    """"""
    return spncci.relative_operator_descriptor(task)

def task_pool(task):
    """"""
    return ("Nsigmamax{Nsigmamax:02d}".format(**task))

##################################################################
# task control
##################################################################

mcscript.task.init(
    task_list,
    task_descriptor=task_descriptor,
    task_pool=task_pool,
    phase_handler_list=[
        spncci.do_generate_relative_operators
        ],
    # Note: change to mcscript.task.archive_handler_hsi for tape backup
    archive_phase_handler_list=[mcscript.task.archive_handler_generic]
    )

################################################################
# termination
################################################################

mcscript.termination()
