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
        "N1v" : 1
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
    phase_handler_list=[spncci.task_handler_relative_tensor_rmes]
    )

################################################################
# termination
################################################################

mcscript.termination()
