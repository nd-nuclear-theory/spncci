""" runmfdn01.py

    See runmfdn.txt for description.

    Mark A. Caprio
    University of Notre Dame

    - 12/14/16 (mac): Created.
    - 12/29/16 (mac): Complete run.  Add full run list.
      
"""

import mcscript
import spncci

# initialize mcscript
mcscript.init()

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
        "Nsigma_ex_max" : 4,
        "num_eigenvalues" : 10,
        "J_values" : [0,1]
    }
    for Nmax in mcscript.utils.value_range(0,4,2)
]

################################################################
# run control
################################################################

def task_descriptor(task):
    """"""
    return ("Nmax{Nmax:02d}".format(**task))

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
    phase_handler_list=[spncci.generate_lsu3shell_rmes,spncci.call_spncci]
    )

################################################################
# termination
################################################################

mcscript.termination()
