"""runsu3rme01.py

  Example run script to calculate RMEs of relative unit tensors and
  the symplectic operators (A/B/Nintr) in SU(3)-NCSM many-body basis.

  Example invocation (under csh):

    foreach n (00 02 04 06 08 10)
      qsubm su3rme01 long 999 --pool="Nsigmamax${n}" --ranks=24
    end

  Then manually save the results:

    cd runsu3rme01/results
    cp *.tgz ${SPNCCI_LSU3SHELL_DIR}

  Language: Python 3

  A. E. McCoy and M. A. Caprio
  Department of Physics, University of Notre Dame

  - 12/1/16 (aem): Created.
  - 5/26/17 (mac):
      + Rename from generate_lsu3shell_rmes.py to runsu3rme01.py.
      + Move task handler out to spncci.py.
  - 5/30/17 (mac): Split out creation of unit tensors from evaluation of rmes.

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

Nsigma_max_list = mcscript.utils.value_range(0,10,2)

task_list = [
    {
        "nuclide" : (3,3),
        ## "Nmax" : Nmax,
        "Nstep" : 2,
        "N1v" : 1,
        "Nsigma_0" : 11,
        "Nsigma_max" : Nsigma_max,
        "J0" : -1,  # -1 for no restriction (needed for spncci); 0 for only Hamiltonian like operators
        "su3rme_descriptor_template": spncci.su3rme_descriptor_template_Nsigmamax
    }
    for Nsigma_max in Nsigma_max_list
]

################################################################
# run control
################################################################

def task_descriptor(task):
    """"""
    return (spncci.su3rme_descriptor_template_Nsigmamax.format(**task))

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
        spncci.do_generate_lsu3shell_rmes,
        spncci.retrieve_su3rme_files  # for debugging: to test of file retrieval
        ],
    # Note: change to mcscript.task.archive_handler_hsi for tape backup
    archive_phase_handler_list=[mcscript.task.archive_handler_generic]
    )

################################################################
# termination
################################################################

mcscript.termination()
