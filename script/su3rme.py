"""su3rme.py -- define scripting for lsu3shell generation of su3rmes

    Environment variables:

        SPNCCI_PROJECT_ROOT_DIR -- root directory under which
            lsu3shell and spncci codes are found

        SPNCCI_OPERATOR_DIR -- base directory for relative operator
            files (colon-delimited search path); operator files will
            be sought within subdirectories named in
            spncci.operator_subdirectory_list

        SPNCCI_SU3RME_DIR -- base directory for su3rme files
            (colon-delimited search path); su3rme files will be sought
            within subdirectories named in su3rme.su3rme_subdirectory_list

        SU3SHELL_DATA -- directory needed for running lsu3shell.  
            within directory there must be a subdirectory rme.

        Ex (personal workstation):
            setenv SPNCCI_PROJECT_ROOT_DIR "${HOME}/code"
            setenv SPNCCI_OPERATOR_DIR ${HOME}/data/spncci/operator
            setenv SPNCCI_SU3RME_DIR ${HOME}/data/spncci/su3rme

        Ex (NDCRC nuclthy):
            setenv SPNCCI_PROJECT_ROOT_DIR "${HOME}/code"
            setenv SPNCCI_OPERATOR_DIR ${HOME}/data/spncci/operator:/afs/crc.nd.edu/group/nuclthy/data/spncci/operator
            setenv SPNCCI_SU3RME_DIR ${HOME}/data/spncci/su3rme:/afs/crc.nd.edu/group/nuclthy/data/spncci/su3rme


        Ex (NERSC m2032):
            setenv SPNCCI_PROJECT_ROOT_DIR "${HOME}/code"
            setenv SPNCCI_OPERATOR_DIR ${HOME}/data/spncci/operator:/project/projectdirs/m2032/data/spncci/operator
            setenv SPNCCI_SU3RME_DIR /project/projectdirs/m2032/data/spncci/su3rme
            setenv SPNCCI_SU3RME_DIR ${SPNCCI_SU3RME_DIR}:${SCRATCH}/data/spncci/su3rme-expanded:${CSCRATCH}/data/spncci/su3rme-expanded
            setenv PYTHONPATH ${SPNCCI_PROJECT_ROOT_DIR}/spncci/script:${PYTHONPATH}

        Example of pre-expanded rme files (see script/su3rme-untar.csh):

            edison09:/global/cscratch1/sd/mcaprio/data/spncci/su3rme-expanded% ls runmac0408/
            su3rme-Z03-N03-Nsigmamax00-Nstep2/  su3rme-Z03-N03-Nsigmamax04-Nstep2/  su3rme-Z03-N03-Nsigmamax08-Nstep2/
            su3rme-Z03-N03-Nsigmamax02-Nstep2/  su3rme-Z03-N03-Nsigmamax06-Nstep2/

            edison09:/global/cscratch1/sd/mcaprio/data/spncci/su3rme-expanded% ls runmac0408/su3rme-Z03-N03-Nsigmamax00-Nstep2/
            Arel.rme                  relative_unit_000007.rme  relative_unit_000020.rme  relative_unit_000033.rme  relative_unit_000046.rme
            ...


    You will also need the directory containing the present script
    file (spncci.py) to be in your Python path:

        setenv PYTHONPATH ${SPNCCI_PROJECT_ROOT_DIR}/spncci/script:${PYTHONPATH}

    Task parameters:

        # space parameters
        nuclide (tuple of int): (N,Z)
        Nmax (int): oscillator Nmax for many-body basis
        Nstep (int): step in N for many-body basis (1 or 2)
        N1v (int): valence shell oscillator N
        ## Nsigma_0 (float): U(1) quantum number of lowest configuration
        ##     (in general can be half integer since contains zero-point offset)
        Nsigma_max (int): maximum oscillator exitation for LGIs

        # su3rme parameters
        J0 (int): restriction on J0 for unit tensors (normally -1 to include
            all J0 as needed for spncci recurrence)
        "su3rme_descriptor_template" (str): template for string describing SU(3)-NCSM
            space truncation used in SU3RME calculation
        "su3rme_mode" (str): mode argument for SU3RME

        # eigenspace parameters -- for spncci runs ONLY
        num_eigenvalues (int): number of eigenvalues to compute in each J space
        J_range (tuple of float/int): (min,max,step) for J spaces
        hw_range (tuple of float): (min,max,step) for oscillator hw
        ...
          
  Language: Python 3

  A. E. McCoy and M. A. Caprio
  University of Notre Dame and TRIUMF

  SPDX-License-Identifier: MIT

  2/16/18 (aem): Created with code extracted from spncci.py.

"""
import glob
import mcscript
import os

import lsu3shell_operator as lsu3
################################################################
# global configuration
################################################################

# environment configuration variables
project_root = os.environ["SPNCCI_PROJECT_ROOT_DIR"]
operator_directory_list = os.environ["SPNCCI_OPERATOR_DIR"].split(":")
operator_subdirectory_list = []
su3rme_directory_list = os.environ["SPNCCI_SU3RME_DIR"].split(":")
su3rme_subdirectory_list = []

# executable files
# ... from lsu3shell
su3rme_executable = os.path.join(project_root,"lsu3shell","programs","tools","SU3RME_MPI")
su3basis_executable = os.path.join(project_root,"lsu3shell","programs","tools","ncsmSU3xSU2IrrepsTabular")
## su3rme_executable = spncci.get_lsu3shell_executable("SU3RME_MPI")
## su3basis_executable = spncci.get_lsu3shell_executable("ncsmSU3xSU2IrrepsTabular")

# ... from spncci
generate_lsu3shell_model_space_executable = os.path.join(project_root,"spncci","programs","unit_tensors","generate_lsu3shell_model_space")
u3s_subspace_lister_executable=os.path.join(project_root,"spncci","programs","lgi","get_u3s_subspaces")
# TODO (mac): update paths to proper installation directory...
## generate_lsu3shell_model_space_executable = spncci.get_spncci_executable("generate_lsu3shell_model_space")
## u3s_subspace_lister_executable = spncci.get_spncci_executable("get_u3s_subspaces")



################################################################
# LSU3Shell basis setup
################################################################
def generate_model_space_files(task):
    """
    Create LSU3shell model space files for SU3RME:
        model_space_bra.dat
        model_space_ket.dat
    
    If no model space file is given (for bra or ket) in task, then generate 
    model space for full space using generate_lsu3shell_model_space_executable.

    Otherwise, copy truncated model space to model space input files: 
        model_space_bra.dat
        model_space_ket.dat
    """
    model_space_executable_command_line=[
        generate_lsu3shell_model_space_executable,
        "{nuclide[0]:d}".format(**task),
        "{nuclide[1]:d}".format(**task),
        "{Nsigma_max:d}".format(**task),
        "{Nstep:d}".format(**task)
    ]

    if task["model_space_file_bra"]==None:
        print("generating model_space.dat")
        
        mcscript.call(command_line,mode=mcscript.CallMode.kSerial)
        mcscript.call(["cp","model_space.dat","model_space_bra.dat"])

    else:
        mcscript.call(["cp",task["model_space_file_bra"], "model_space_bra.dat"])

    if task["model_space_file_ket"]==None:
        mcscript.call(command_line,mode=mcscript.CallMode.kSerial)
        mcscript.call(["cp","model_space.dat","model_space_ket.dat"])

    else:
        mcscript.call(["cp",task["model_space_file_ket"], "model_space_ket.dat"])



def generate_basis_table(task):
    """Create SU(3)-NCSM basis table for bra and ket modle spaces.

    Invokes ncsmSU3xSU2IrrepsTabular.

    """
    model_space_filename = "model_space_ket.dat".format(**task)
    basis_listing_filename = "lsu3shell_basis_ket.dat"

    command_line=[su3basis_executable,model_space_filename,basis_listing_filename]
    mcscript.call(
        command_line,
        mode=mcscript.CallMode.kSerial
    )


    model_space_filename = "model_space_bra.dat".format(**task)
    basis_listing_filename = "lsu3shell_basis_bra.dat"

    command_line=[su3basis_executable,model_space_filename,basis_listing_filename]
    mcscript.call(
        command_line,
        mode=mcscript.CallMode.kSerial
    )

################################################################
# relative operator SU(3) RME construction
################################################################
# su3rme descriptor string: describes the many-body space on which SU(3) 
#   RMEs are calculated.  It is used, e.g., in the filename of the archive file
#   in which the SU(3) RMEs are stored.
# descriptor string for straightforward case of pure Nsigmamax truncation
su3rme_descriptor_template_Nsigmamax = "Z{nuclide[0]:02d}-N{nuclide[1]:02d}-Nsigmamax{Nsigma_max:02d}-Nstep{Nstep:d}"


def setup_su3shell_directories(task):
    """
    Make sure temporary work directory exists
    """
    su3rme_data_directory=os.environ["SU3SHELL_DATA"]
    directory_name=su3rme_data_directory+"/rme"
    mcscript.call(["mkdir","-p",directory_name])

def calculate_rmes(task):
    """ Invoke lsu3shell SU3RME code to calculate rmes of relative unit
    tensors and symplectic raising/lowering/N operators in SU(3)-NCSM basis.

    Invokes SU3RME_MPI.
    """
    if ("su3rme_mode" not in task):
        task["su3rme_mode"] = "text"

    # call SU3RME
    command_line = [
        su3rme_executable,
        "model_space_bra.dat",
        "model_space_ket.dat",
        "relative_operators.dat",
        task["su3rme_mode"]
    ]
    mcscript.call(
        command_line,
        mode=mcscript.CallMode.kHybrid
    )


def save_su3rme_files(task):
    """Create archive of SU(3) RMEs of relative operators.

    Some auxiliary files (e.g., the list of operators) are saved as well.

    Manual follow-up: The rme files are bundled into tgz files and saved to
    the current run's results directory.  They should then be moved
    (manually) to the directory specified in SPNCCI_LSU3SHELL_DIR, for
    subsequent use.

    """

    su3rme_descriptor = task["su3rme_descriptor_template"].format(**task)

    # select files to save
    archive_file_list = [
        "model_space_ket.dat",
        "model_space_bra.dat",
        "relative_operators.dat",
        "lsu3shell_basis.dat",
        "relative_unit_tensor_labels.dat"
    ]
    archive_file_list += glob.glob('*.rme')

    # generate archive
    archive_filename = f"su3rme-{su3rme_descriptor}.tgz"
    mcscript.call(["tar", "-zcvf", archive_filename] + archive_file_list)

    # save independent copy of basis listing outside tarball for easy inspection
    basis_filename = f"lsu3shell_basis_{su3rme_descriptor}.dat"
    mcscript.call(["cp", "lsu3shell_basis.dat", basis_filename])

    # move archive to results directory (if in multi-task run)
    if (mcscript.task.results_dir is not None):
        mcscript.call(
            [
                "mv",
                "--verbose",
                basis_filename,archive_filename,
                f"--target-directory={mcscript.task.results_dir}"
            ]
        )


def retrieve_su3rme_files(task):
    """ Retrieve archive of relative operator SU(3) RME files.

    (1) Directory is symlinked as a subdirectory named lsu3shell_rme, or...
    (2) Files are retrieved into a subdirectory named lsu3shell_rme.
    """
    # identify su3rme data directory
    su3rme_descriptor = task["su3rme_descriptor_template"].format(**task)
    directory_name = mcscript.utils.search_in_subdirectories(
        su3rme_directory_list,
        su3rme_subdirectory_list,
        f"su3rme-{su3rme_descriptor}",
        error_message="Data directory for SU(3) RMEs not found",
        fail_on_not_found=False
    )
    archive_filename = mcscript.utils.search_in_subdirectories(
        su3rme_directory_list,
        su3rme_subdirectory_list,
        f"su3rme-{su3rme_descriptor}.tgz",
        error_message="Archive file for SU(3) RMEs not found",
        fail_on_not_found=False
    )

    if (directory_name is not None):

        # remove any existing symlink or data directory
        #
        # Notes: On a symlink to a directory: rmdir fails; rm or "rm -r"
        # removes symlink.  But "rm -r" will also work if tar file had
        # been directly expanded before and needs to be replaced by a
        # symlink.
        if (os.path.exists("lsu3shell_rme")):
            mcscript.call(["rm","-r","lsu3shell_rme"])

        # link to data su3rme directory
        mcscript.call(
            [
                "ln",
                "-s",
                directory_name,
                "lsu3shell_rme"
            ]
        )

    elif (archive_filename is not None):

        # set up data directory
        if (not os.path.exists("lsu3shell_rme")):
            mcscript.utils.mkdir("lsu3shell_rme")

        # extract archive contents
        mcscript.call(
            [
                "tar",
                "-xvf",
                archive_filename,
                "--directory=lsu3shell_rme"
            ]
        )

    else:
        raise(mcscript.exception.ScriptError("Cannot find SU(3) RME data"))


def generate_lsu3shell_rmes(task):
    """
    Generates LSU3Shell rmes for operators given in relative_operators.dat
    and saves them to file

    """
    # set up necessary directories for lsu3shell
    setup_su3shell_directories(task)

    # generate model space files for bra and ket
    generate_model_space_files(task)

    # generate basis listing for basis defined by bra and ket model space
    generate_basis_table(task)

    # generate operators rmes
    calculate_rmes(task)    

def do_generate_lsu3shell_generator_rmes(task,save=True):
    """
    Control code for generating RMEs in the SU(3)-NCSM basis, for relative
    unit tensors and symplectic raising/lowering/N operators.

    input:
        task : task dictionary for specific task
        save (optional, bool) : if save is True, then save files and do cleanup
            su3rmes saved to "su3rme-{descriptor}.tgz", where descriptor
            is defined by task["su3rme_descriptor_template"]
    """

    # retrieve relevant operator files
    retrieve_generator_operator_files(task)
    mcscript.call(["ln","-sf","relative_generators.dat","relative_operators.dat"])
    generate_lsu3shell_rmes(task)

    if save:
        # save results
        save_su3rme_files(task)

        # clean up working directory
        mcscript.call(["du","-hs","."])  # log working directory disk usage
        delete_filenames=glob.glob('*')
        mcscript.call(["rm"] + delete_filenames)


def do_generate_lsu3shell_unittensor_rmes(task,save=True):
    """
    Control code for generating RMEs in the SU(3)-NCSM basis, for relative
    unit tensors and symplectic raising/lowering/N operators.

    input:
        task : task dictionary for specific task
        save (optional, bool) : if save is True, then save files and do cleanup
            su3rmes saved to "su3rme-{descriptor}.tgz", where descriptor
            is defined by task["su3rme_descriptor_template"]
    """

    # retrieve relevant operator files
    retrieve_unit_operator_files(task)
    split_relative_operator_file(task)
    mcscript.call(["ln","-sf","relative_unittensors.dat","relative_operators.dat"])
    generate_lsu3shell_rmes(task)
    
    if save:
        # save results
        save_su3rme_files(task)

        # clean up working directory
        mcscript.call(["du","-hs","."])  # log working directory disk usage
        delete_filenames=glob.glob('*')
        mcscript.call(["rm"] + delete_filenames)


def do_generate_lsu3shell_rmes(task):
    ## TODO: Test updated scripting.  Use runsu3rme03
    """
    Control code for generating RMEs in the SU(3)-NCSM basis, for relative
    unit tensors and symplectic raising/lowering/N operators.

    su3rme files saved to "su3rme-{descriptor}.tgz", where descriptor
    is defined by task["su3rme_descriptor_template"]
    """
    ## retrieve relevant operator files
    retrieve_operator_files(task)
    
    ## Compute su3rmes for generators 
    do_generate_lsu3shell_generator_rmes(task,save=False)

    ## Comput su3rmes for unit operators 
    do_generate_lsu3shell_unittensor_rmes(task,save=False)

    # save results
    save_su3rme_files(task)

    # clean up working directory
    mcscript.call(["du","-hs","."])  # log working directory disk usage
    delete_filenames=glob.glob('*')
    mcscript.call(["rm"] + delete_filenames)
    

if (__name__ == "__MAIN__"):
    pass
