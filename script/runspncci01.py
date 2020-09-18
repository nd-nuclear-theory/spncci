""" runspncci01.py

    spncci demonstration run

    DEPRECATED -- outdated syntax

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
import spncci
import os
# initialize mcscript
mcscript.init()

##################################################################
# directory configuration
##################################################################

interaction_directory = os.environ["SPNCCI_INTERACTION_DIR"]
unit_tensor_directory = os.environ["SPNCCI_LSU3SHELL_DIR"]
interaction_filename_template = os.path.join(interaction_directory,"JISP16_Nmax20","JISP16_Nmax20_hw{:2.1f}_rel.dat")
unit_tensor_directory_template = os.path.join(unit_tensor_directory,"lsu3shell_{Nsigma_ex_max:02d}")  # TODO label by N and Z; no longer "directory"
spncci.seed_subdirectory_list += [
"runaem0031", #6Li
"runaem0079", #7Be
"runaem0080" #7Li
] 

##################################################################
# build task list
##################################################################

task_list = [
    {
        "nuclide" : (3,3),
        "Nmax" : Nmax,
        "Nstep" : 2,
        "N1v" : 1,
        "Nsigma_max" : Nsigma_ex_max,
        "J_range" : (1,3,2), #min, max, step
        "hw_range" : (20,20,2.5), # min, max, step
        "seed_descriptor_template" : spncci.seed_descriptor_template_Nsigmamax,
        "hyperblocks_dir" : None,
        "interaction" : "JISP16",
        "interaction_filename_template" :interaction_filename_template,
        "use_coulomb" : True,
        "coulomb_filename" : "coulomb_Nmax20_steps500_rel.dat",
        "observables" : [("r2intr",0),("Qintr",2),("Qpintr",2),("Qnintr",2)],
        "num_eigenvalues" : 5,
        "truncation_filename": None,
        "transformation_filename": None,
        # eigensolver convergence parameters
        "eigensolver_num_convergence" : 2*5+1, # docs for Spectra SymEigsSolver say to take "ncv>=2*nev" and want smaller than dim of matrix 
        "eigensolver_max_iterations" : 1000, # at least 100 times num_eigenvalues, possible more
        "eigensolver_tolerance" : 1e-6 
    }
    for Nsigma_ex_max in mcscript.utils.value_range(0,6,2)  # CAVEAT: Nmax0 requires special treatment for num eigenvectors
    for Nmax in mcscript.utils.value_range(Nsigma_ex_max,20,2)
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
