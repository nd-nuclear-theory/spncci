""" runspncci01.py

    spncci demonstration run

    Environment:
      See spncci.py header

    Mark A. Caprio and Anna E. McCoy
    University of Notre Dame

    + 1/8/17 (mac): Created.
    + 2/23/17 (mac): Implement basic scripting.
    + 9/24/20 (aem): Updated to reflect current scripting

"""

import mcscript
import seeds
import spncci
import os
# initialize mcscript
mcscript.init()

##################################################################
# data file search paths
##################################################################
## spncci.py searches for interaction file in subdirectorys 
## spncci.interaction_subdirectory_list in directory SPNCCI_INTERACTION_DIR
spncci.interaction_subdirectory_list += ["",
"JISP16_Nmax20",
"Daejeon16_Nmax40",
"N3LO"
]

#Template for specific interaction filename
interaction_filename_template = "JISP16_Nmax20_hw{hw:2.1f}_rel.dat"

## subdirectories containing seeds for spncci recurrence 
## in directory SPNCCI_SEED_DIR
seeds.seed_subdirectory_list += [
     "runaem0031", #6Li
     "runaem0079", #7Be
     "runaem0080"  #7Li
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
        "Nsigma_max" : Nsigma_max,
        "J_range" : (1,3,2), #min, max, step
        "hw_range" : (20,20,2.5), # min, max, step
        "seed_descriptor_template" : spncci.seed_descriptor_template_Nsigmamax,
        ## If hyperblocks_dir set to none, temporary directory created in run directory
        "hyperblocks_dir" : None,
        # Tag for interaction
        "interaction" : "JISP16",
        # Interaction filename template
        "interaction_filename_template" :interaction_filename_template,
        "use_coulomb" : True,
        "coulomb_filename" : "coulomb_Nmax20_steps500_rel.dat",
        # "observables" : [(<Observable name>,J)]
        "observables" : [("r2intr",0),("Qintr",2),("Qpintr",2),("Qnintr",2)],
        "num_eigenvalues" : 5,
        #parameters for spncci basis truncation
        "truncation_filename": None,
        "transformation_filename": None,
        # eigensolver convergence parameters
        "eigensolver_num_convergence" : 2*5+1, # docs for Spectra SymEigsSolver say to take "ncv>=2*nev" and want smaller than dim of matrix 
        "eigensolver_max_iterations" : 100, # at least 100 times num_eigenvalues, possible more
        "eigensolver_tolerance" : 1e-6 
    }
    for Nsigma_max in mcscript.utils.value_range(0,6,2)  
    for Nmax in mcscript.utils.value_range(Nsigma_max,20,2)
]

################################################################
# run control
################################################################

def task_descriptor(task):
    """"""
    coulomb = int(task["use_coulomb"])
    return ("Z{nuclide[0]:d}-N{nuclide[1]:d}-{interaction:s}-{coulomb:1d}-Nsigmamax{Nsigma_max:02d}-Nmax{Nmax:02d}".format(coulomb=coulomb,**task))

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
