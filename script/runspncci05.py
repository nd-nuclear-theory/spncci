"""runspncci05.py

    spncci demonstration run

    See spncci/script/spncci.py for proper settings for environment
    variables.

    Mark A. Caprio
    University of Notre Dame

    + 1/8/17 (mac): Created.
    + 2/23/17 (mac): Implement basic scripting.
    + 5/30/17 (mac): Update handling of SU(3) RME storage.

"""

import mcscript
import os
import su3rme
import seeds
import spncci

# initialize mcscript
mcscript.init()

##################################################################
# data file search paths
##################################################################

spncci.seed_subdirectory_list += [
     "runaem0031"  # 6Li -- binary
    ]

spncci.interaction_subdirectory_list += ["","JISP16_Nmax20"]

spncci.truncation_subdirectory = "6Li"
##################################################################
# build task list
##################################################################

# interaction filename template
#
# can use dummy variable hw in format specification
interaction_filename_template = "JISP16_Nmax20_hw{hw:2.1f}_rel.dat"

task_list = [
    {
        # space parameters
        "nuclide" : (3,3),
        "Nmax" : Nmax,
        "Nstep" : 2,
        "N1v" : 1,
        "Nsigma_max" : Nsigma_max,

        # su3rme parameters
        # "J0" : 0,
        "su3rme_descriptor_template" : su3rme.su3rme_descriptor_template_Nsigmamax,
        "seed_descriptor_template" : seeds.seed_descriptor_template_Nsigmamax,
        # "hyperblocks_dir" examples : "/tmp","/scratch/hyperblocks", None
        "hyperblocks_dir" : None,

        # spncci parameters
        "J_range" : (1,3,2), #min, max, step
        "hw_range" : (20,20,2.5), # min, max, step
        "interaction" : "JISP16",
        "interaction_filename_template" : interaction_filename_template,
        "use_coulomb" : False,
        "observables" : [("r2intr",0),("Qintr",2)],
        "num_eigenvalues" : 10,
        "coulomb_filename" : "coulomb_Nmax20_steps500_rel.dat",
        "truncation_filename": None,
        "transformation_filename": None,

        # eigensolver convergence parameters
        "eigensolver_num_convergence" : 2*10,  # docs for Spectra SymEigsSolver say to take "ncv>=2*nev"
        "eigensolver_max_iterations" : 100*10,  # at least 100 times num_eigenvalues, possible more
        "eigensolver_tolerance" : 1e-8
    }
    for Nsigma_max in mcscript.utils.value_range(0,6,2)  # CAVEAT: Nmax0 requires special treatment for num eigenvectors
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
