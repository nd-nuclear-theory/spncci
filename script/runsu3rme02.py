"""runsu3rme02.py

  Parallel timing testbed for su3rme code, based on Nv=0 spaces.

  CAVEAT: Pools are based on MPI allocation.  Runtimes within a pool
  are inhomogeneous (and increasing).  So the mcscript.task graceful
  termination based on walltime will not prevent walltime timeouts.

  Example invocation (under csh):

    foreach nslaves (01 04 16 64)
      @ nranks = ${nslaves} + 1
      qsubm su3rme02 long 2880 --pool="nslaves${nslaves}" --ranks=${nranks}
    end

    foreach nslaves (01 04 16 64)
      @ nranks = ${nslaves} + 1
      qsubm su3rme02 long 2880 --pool="nslaves${nslaves}" --ranks=${nranks}
    end


  Language: Python 3

  A. E. McCoy and M. A. Caprio
  Department of Physics, University of Notre Dame

  - 6/3/17 (mac): Created, based on runsu3rme01.

"""

import glob
import os
import sys

import mcscript
import spncci

mcscript.init()

##################################################################
# data file search paths
##################################################################

spncci.operator_subdirectory_list += ["rununittensor01"]

##################################################################
# build task list
##################################################################

Nsigma_max_list = mcscript.utils.value_range(0,6,2)
nslaves_list = [2**k for k in mcscript.utils.value_range(0,6,2)]
task_list = [
    {
        "nuclide" : (2,2),
        ## "Nmax" : None,
        "Nstep" : 2,
        "N1v" : 0,
        "Nsigma_0" : 11,
        "Nsigma_max" : Nsigma_max,
        "J0" : -1,  # -1 for no restriction (needed for spncci); 0 for only Hamiltonian like operators
        "su3rme_descriptor_template": spncci.su3rme_descriptor_template_Nsigmamax,

        # for timing pool
        "nslaves" : nslaves
    }
    for Nsigma_max in Nsigma_max_list
    for nslaves in nslaves_list
]


## def do_runs(task):
##     """ Do runs iterating over Nmax."""
## 
##     for Nsigma_max in Nsigma_max_list:
##         task["Nsigma_max"] = Nsigma_max
##         spncci.do_generate_lsu3shell_rmes(task)

################################################################
# run control
################################################################

def task_descriptor(task):
    """"""
    return (spncci.su3rme_descriptor_template_Nsigmamax.format(**task))

def task_pool(task):
    """"""
    return ("nslaves{nslaves:02d}".format(**task))

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
