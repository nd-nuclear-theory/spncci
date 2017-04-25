"""generate_lsu3shell_rmes.py

  Calculate each of the relative unit tensors and 
  Brel and Nintr used for obtaining LGI expansion

  Arguments:
    Z N twice_Nsigma_0 Nmax Nstep N1B

  Note: Requires utils subpackage from mcscript.  This should be
  restructured to load it from within the mcscript package once the
  overhaul of the mcscript package structure is complete.

  Language: Python 3

  A. E. McCoy
  Department of Physics, University of Notre Dame

  12/1/16 (aem): Created.
"""

import os
import glob
import sys
import utils
import spncci
import mcscript
mcscript.init()

##################################################################
# directory configuration
##################################################################

interaction_directory = os.environ["SPNCCI_INTERACTION_DIR"]
unit_tensor_directory = os.environ["SPNCCI_LSU3SHELL_DIR"]
interaction_filename_template = os.path.join(interaction_directory,"JISP16_Nmax20","JISP16_Nmax20_hw{:2.1f}_rel.dat")
unit_tensor_filename_template = os.path.join("lsu3shell_Z{nuclide[0]:02d}_N{nuclide[1]:02d}_{Nsigma_ex_max:02d}")  
# TODO label by N and Z; no longer "directory"

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
        "Nsigma_ex_max" : Nmax,
        "num_eigenvalues" : 10,
        "J0" : 0,
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
# get relative unit tensor rmes 
##################################################################
def do_generate_lsu3shell_rmes(task):
    """"""
    spncci.generate_lsu3shell_rmes(task)
    ## archiving directory of rmes 
    unit_tensor_directory=task["unit_tensor_filename_template"].format(**task)
    unit_tensor_directory_archive_filename = "{}.tgz".format(unit_tensor_directory)
    mcscript.call(["tar","-zcvf",unit_tensor_directory_archive_filename, "lsu3shell_rme"])
    # mcscript.call(["rm","-r","lsu3shell_rme"])
##################################################################
# task control
##################################################################

mcscript.task.init(
    task_list,
    task_descriptor=task_descriptor,
    task_pool=task_pool,
    phase_handler_list=[
        do_generate_lsu3shell_rmes
        ],
    # Note: change to mcscript.task.archive_handler_hsi for tape backup
    archive_phase_handler_list=[mcscript.task.archive_handler_generic]
    )

################################################################
# termination
################################################################

mcscript.termination()