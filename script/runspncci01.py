""" runmfdn01.py

    See runmfdn.txt for description.

    Mark A. Caprio
    University of Notre Dame

    - 12/14/16 (mac): Created.
    - 12/29/16 (mac): Complete run.  Add full run list.
      
"""

import mcscript
import spncci
import os
# initialize mcscript
mcscript.init()

##################################################################
# build task list
##################################################################
interaction_directory = os.path.join(os.environ["HOME"],"nuclthy","data","interaction","rel")
interaction_filename_template = os.path.join(interaction_directory,"JISP16_Nmax20","JISP16_Nmax20_hw{:2.1f}_rel.dat")
unit_tensor_directory_template = os.path.join(
        "/afs/crc.nd.edu/group/nuclthy/data/lsu3shell-data/unit_tensors",
        "lsu3shell_{Nsigma_ex_max:02d}"
        )


task_list = [
    {
        "nuclide" : (3,3),
        "Nmax" : Nmax,
        "Nstep" : 2,
        "N1v" : 1,
        "Nsigma_0" : 11,
        "Nsigma_ex_max" : 2,
        "num_eigenvalues" : 10,
        "J0" : 0,
        "J_range" : [1,3,2], #min, max, step
        "hw_range" : [10,30,2.5], # min, max, step
        "interaction_filename_template" :interaction_filename_template,
        "unit_tensor_directory" : unit_tensor_directory_template,
        "observables" : ["r2intr"]
    }
    for Nmax in mcscript.utils.value_range(2,20,2)
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
    phase_handler_list=[
        spncci.do_full_spncci_run
        ]
    )

################################################################
# termination
################################################################

mcscript.termination()
