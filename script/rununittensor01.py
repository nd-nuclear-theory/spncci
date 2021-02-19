"""rununittensor01.py

  Example run script to generate relative unit tensors and symplectic
  operators (Arel/Brel/Nrel) in LSU3shell operator format.

  Example invocation (under csh):

    foreach n (00 01 02 03 04 05 06 07 08 09 10 11)
      qsubm unittensor01 long 999 --pool="Nsigmamax${n}" --num=2
    end

  Then manually save the results:

    cp -rv rununittensor01/results ${SPNCCI_OPERATOR_DIR}/rununittensor01

  (Of course, the shorthand ${SPNCCI_OPERATOR_DIR} only works if
  SPNCCI_OPERATOR_DIR is a single directory, not a list of directories.)

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
Nsigma_max_set = mcscript.utils.value_range(0,11,1)

task_list = [
    {
        "N1v" : N1v,
        "Nsigma_max" : Nsigma_max,
        "Nstep" : 2
    }
    for N1v in N1v_set
    for Nsigma_max in Nsigma_max_set
]

################################################################
# run control
################################################################

def task_descriptor(task):
    """"""
    return su3rme.relative_operator_descriptor(task)

def task_pool(task):
    """"""
    return ("Nsigmamax{Nsigma_max:02d}".format(**task))

##################################################################
# task control
##################################################################

mcscript.task.init(
    task_list,
    task_descriptor=task_descriptor,
    task_pool=task_pool,
    phase_handler_list=[
        su3rme.do_generate_relative_operators
        ],
    # Note: change to mcscript.task.archive_handler_hsi for tape backup
    archive_phase_handler_list=[mcscript.task.archive_handler_generic]
    )

################################################################
# termination
################################################################

mcscript.termination()
