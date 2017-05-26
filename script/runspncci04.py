""" runspncci04.py

    spncci demonstration run

    Environment:
      SPNCCI_INTERACTION_DIR -- base directory for relative lsjt interaction files
      SPNCCI_LSU3SHELL_DIR -- base directory for unit tensor rme files

    Ex (NDCRC nuclthy):
      setenv SPNCCI_INTERACTION_DIR /afs/crc.nd.edu/group/nuclthy/data/interaction/rel
      setenv SPNCCI_LSU3SHELL_DIR /afs/crc.nd.edu/group/nuclthy/data/spncci/lsu3shell

    Ex (NERSC m2032):
      setenv SPNCCI_INTERACTION_DIR /project/projectdirs/m2032/data/interaction/rel
      setenv SPNCCI_LSU3SHELL_DIR /project/projectdirs/m2032/data/spncci/lsu3shell

    Mark A. Caprio
    University of Notre Dame

    + 1/8/17 (mac): Created.
    + 2/23/17 (mac): Implement basic scripting.

"""

import mcscript
import os
import spncci

# initialize mcscript
mcscript.init()

##################################################################
# directory configuration
##################################################################

interaction_directory = os.environ["SPNCCI_INTERACTION_DIR"]
unit_tensor_directory = os.environ["SPNCCI_LSU3SHELL_DIR"]
interaction_filename_template = os.path.join(interaction_directory,"JISP16_Nmax20","JISP16_Nmax20_hw{:2.1f}_rel.dat")
# interaction_filename_template = os.path.join(interaction_directory,"coulomb_Nmax20_rel.dat")

unit_tensor_directory_template = os.path.join(unit_tensor_directory,"lsu3shell_Z{nuclide[0]:02d}_N{nuclide[1]:02d}_{Nsigma_max:02d}")

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
        "Nsigma_max" : Nsigma_max,
        "num_eigenvalues" : 10,
        # "J0" : 0,
        "J_range" : (1,3,2), #min, max, step
        "hw_range" : (20,20,2.5), # min, max, step
        "interaction_directory" : interaction_directory,
        "interaction_filename_template" : interaction_filename_template,
        "use_coulomb" : False,
        "observables" : [("r2intr",0),("Qintr",2)],
        "unit_tensor_directory" : unit_tensor_directory_template
    }
    for Nsigma_max in mcscript.utils.value_range(0,6,2)  # CAVEAT: Nmax0 requires special treatment for num eigenvectors
    for Nmax in mcscript.utils.value_range(Nsigma_max,20,2)
]

################################################################
# run control
################################################################

def task_descriptor(task):
    """"""
    return ("Z{nuclide[0]:d}-N{nuclide[1]:d}-Nsigmaexmax{Nsigma_max:02d}-Nmax{Nmax:02d}".format(**task))

def task_pool(task):
    """"""
    return ("{Nsigma_max:02d}-{Nmax:02d}".format(**task))


##################################################################
# task control
##################################################################

mcscript.task.init(
    task_list,
    task_descriptor=task_descriptor,
    task_pool=task_pool,
    phase_handler_list=[
        spncci.do_full_spncci_run
        ],
    # Note: change to mcscript.task.archive_handler_hsi for tape backup
    archive_phase_handler_list=[mcscript.task.archive_handler_generic]
    )

################################################################
# termination
################################################################

mcscript.termination()
