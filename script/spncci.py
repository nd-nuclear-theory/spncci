"""spncci.py -- define scripting for spncci runs

    Environment variables:

        SPNCCI_PROJECT_ROOT_DIR -- root directory under which
            lsu3shell and spncci codes are found

        SPNCCI_OPERATOR_DIR -- base directory for relative operator
            files (colon-delimited search path); operator files will
            be sought within subdirectories named in
            spncci.operator_subdirectory_list

        SPNCCI_SU3RME_DIR -- base directory for su3rme files
            (colon-delimited search path); su3rme files will be sought
            within subdirectories named in spncci.su3rme_subdirectory_list

        SPNCCI_INTERACTION_DIR -- base directory for relative lsjt
            interaction files (colon-delimited search path);
            interaction files will be sought within subdirectories
            named in spncci.interaction_subdirectory_list

        Ex (personal workstation):
            setenv SPNCCI_PROJECT_ROOT_DIR "${HOME}/code"
            setenv SPNCCI_OPERATOR_DIR ${HOME}/data/spncci/operator
            setenv SPNCCI_SU3RME_DIR ${HOME}/data/spncci/su3rme
            setenv SPNCCI_INTERACTION_DIR ${HOME}/data/interaction/rel

        Ex (NDCRC nuclthy):
            setenv SPNCCI_PROJECT_ROOT_DIR "${HOME}/code"
            setenv SPNCCI_OPERATOR_DIR ${HOME}/data/spncci/operator:/afs/crc.nd.edu/group/nuclthy/data/spncci/operator
            setenv SPNCCI_SU3RME_DIR ${HOME}/data/spncci/su3rme:/afs/crc.nd.edu/group/nuclthy/data/spncci/su3rme
            setenv SPNCCI_INTERACTION_DIR ${HOME}/data/interaction/rel:/afs/crc.nd.edu/group/nuclthy/data/interaction/rel

        Ex (NERSC m2032):
            setenv SPNCCI_PROJECT_ROOT_DIR "${HOME}/code"
            setenv SPNCCI_INTERACTION_DIR ${HOME}/data/interaction/rel:/project/projectdirs/m2032/data/interaction/rel
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
  Department of Physics, University of Notre Dame

  1/8/17 (aem,mac): Created with code from compute_relative_tensors_lsu3shell_rmes.py.
  2/23/17 (aem,mac): Update rme invocation and add spncci handler.
  5/26/17 (mac):
      + Fix notation "Nsigma_ex_max" to "Nsigma_max".
      + Add SPNCCI_PROJECT_ROOT_DIR environment variable.
      + Remove LSU3shell load files.
  5/20/17 (mac): Split out generation of relative operators and calculation of SU(3)
      RMEs.
  6/4/17 (mac): Add search paths for input files.
  6/27/17 (mac): Split out copy of lsu3shell basis listing from su3rme tarball.
  7/4/17 (mac): Add support for using pre-staged expanded archive of su3rmes.

"""
  
import glob
import mcscript
import os

################################################################
# global configuration
################################################################

# environment configuration variables
project_root = os.environ["SPNCCI_PROJECT_ROOT_DIR"]
interaction_directory_list = os.environ["SPNCCI_INTERACTION_DIR"].split(":")
interaction_subdirectory_list = []
operator_directory_list = os.environ["SPNCCI_OPERATOR_DIR"].split(":")
operator_subdirectory_list = []
su3rme_directory_list = os.environ["SPNCCI_SU3RME_DIR"].split(":")
su3rme_subdirectory_list = []

# executable files
# ... from lsu3shell
recoupler_executable = os.path.join(project_root,"lsu3shell","programs","upstreams","RecoupleSU3Operator")
su3rme_executable = os.path.join(project_root,"lsu3shell","programs","tools","SU3RME_MPI")
su3basis_executable = os.path.join(project_root,"lsu3shell","programs","tools","ncsmSU3xSU2IrrepsTabular")
# ... from spncci
generate_lsu3shell_model_space_executable = os.path.join(project_root,"spncci","programs","unit_tensors","generate_lsu3shell_model_space")
generate_lsu3shell_relative_operators_executable = os.path.join(project_root,"spncci","programs","unit_tensors","generate_lsu3shell_relative_operators")
generate_relative_operator_rmes_executable = os.path.join(project_root,"spncci","programs","operators","generate_relative_u3st_operators")
spncci_executable_dir = os.path.join(project_root,"spncci","programs","spncci")


################################################################
# relative operator construction (A-independent)
################################################################

def relative_operator_descriptor(task):
    """Generate descriptor string for use in relative operator archive
    filename.

    Returns:
        (string): descriptor
    """

    descriptor = "Nv{N1v:d}-Nsigmamax{Nsigma_max:02d}-Nstep{Nstep:d}".format(**task)
    return descriptor

def generate_relative_operators(task):

    """Create recoupler input files for relative unit
    tensors and symplectic raising/lowering/N operators.

    Invokes generate_lsu3shell_relative_operators.
    """

    command_line = [
        generate_lsu3shell_relative_operators_executable,
        "{Nsigma_max:d}".format(**task),
        "{Nstep:d}".format(**task),
        "{N1v:d}".format(**task),
        "-1",# All J0
        "-1"# All T0
    ]
    mcscript.call(
        command_line,
        mode=mcscript.CallMode.kSerial
    )

def read_relative_operator_basenames(task):
    """ Read list of relative operators basenames.

    Returns:
        (list) : list of relative operators basenames
    """

    relative_operator_filename = "relative_operators.dat"

    # read list of unit tensors
    relative_operator_stream = open(relative_operator_filename,mode="rt")
    relative_operator_basename_list = [
        ## line.strip()  # remove trailing newline
        (line.split())[0]  # basename is leading element on line
        for line in relative_operator_stream
    ]
    relative_operator_stream.close()
    return relative_operator_basename_list

def recouple_operators(task,relative_operator_basename_list):
    """ Invoke lsu3shell recoupler code on relative unit
    tensors and symplectic raising/lowering/N operators.

    Invokes RecoupleSU3Operator.

    Arguments:
        relative_operator_basename_list (list) : list of operator names
    """

    # iterate over unit tensors
    for basename in relative_operator_basename_list:

        # call recoupler
        command_line = [
            recoupler_executable,
            "{}.recoupler".format(basename),
            basename
        ]
        mcscript.call(
            command_line,
            mode=mcscript.CallMode.kSerial
        )


def save_operator_files(task):
    """ Create archive of relative operator files.

    Manual follow-up: The operator files are bundled into tgz files and saved to
    the current run's results directory.  They should then be moved
    (manually) to the directory specified in SPNCCI_LSU3SHELL_DIR, for
    subsequent use.
    """

    # select files to save
    archive_file_list = glob.glob('*.dat')
    archive_file_list += glob.glob('*.PN')
    archive_file_list += glob.glob('*.PPNN')

    # generate archive
    descriptor = relative_operator_descriptor(task)
    archive_filename = "relative-operators-{}.tgz".format(descriptor)
    mcscript.call(
        [
            "tar", "-zcvf", archive_filename
        ] + archive_file_list
    )

    # move archive to results directory (if in multi-task run)
    if (mcscript.task.results_dir is not None):
        mcscript.call(
            [
                "mv",
                "--verbose",
                archive_filename,
                "--target-directory={}".format(mcscript.task.results_dir)
            ]
        )

def retrieve_operator_files(task):
    """ Retrieve archive of relative operator files.
    """

    # identify archive file
    descriptor = relative_operator_descriptor(task)
    archive_filename = mcscript.utils.search_in_subdirectories(
        operator_directory_list,
        operator_subdirectory_list,
        "relative-operators-{}.tgz".format(descriptor),
        error_message="relative operator archive file not found"
    )

    # extract archive contents
    mcscript.call(
        [
            "tar", "xf", archive_filename
        ]
    )


def do_generate_relative_operators(task):
    """Control code for generating relative unit tensors and symplectic
    raising/lowering/N operators.

    """

    # generate operators and their rmes
    generate_relative_operators(task)
    relative_operator_basename_list = read_relative_operator_basenames(task)
    recouple_operators(task,relative_operator_basename_list)

    ## # do cleanup
    ## deletable_filenames=glob.glob('*.recoupler')
    ## mcscript.call(["rm"] + deletable_filenames)

    # save results
    save_operator_files(task)


################################################################
# relative operator SU(3) RME construction
################################################################

# su3rme descriptor string
#
# This describes the many-body space on which SU(3) RMEs are being
# calculated.  It is used, e.g., in the filename of the archive file
# in which the SU(3) RMEs are stored.

# descriptor string for straightforward case of pure Nsigmamax truncation
su3rme_descriptor_template_Nsigmamax = "Z{nuclide[0]:02d}-N{nuclide[1]:02d}-Nsigmamax{Nsigma_max:02d}-Nstep{Nstep:d}"

def generate_model_space_file(task):
    """Create LSU3shell model space file for SU3RME.

    Invokes generate_lsu3shell_model_space.
    """

    command_line = [
        generate_lsu3shell_model_space_executable,
        "{nuclide[0]:d}".format(**task),
        "{nuclide[1]:d}".format(**task),
        "{Nsigma_max:d}".format(**task),
        "{Nstep:d}".format(**task)
    ]
    mcscript.call(
        command_line,
        mode=mcscript.CallMode.kSerial
    )

def calculate_rmes(task):
    """ Invoke lsu3shell SU3RME code to calculate rmes of relative unit
    tensors and symplectic raising/lowering/N operators in SU(3)-NCSM basis.

    Invokes SU3RME_MPI.
    """

    model_space_filename = "model_space.dat".format(**task)

    if ("su3rme_mode" not in task):
        task["su3rme_mode"] = "text"

    # call SU3RME
    command_line = [
        su3rme_executable,
        model_space_filename,
        model_space_filename,
        "relative_operators.dat",
        task["su3rme_mode"]
    ]
    mcscript.call(
        command_line,
        mode=mcscript.CallMode.kHybrid
    )


def generate_basis_table(task):
    """Create SU(3)-NCSM basis table.

    Invokes ncsmSU3xSU2IrrepsTabular.

    Depends on model space file created by generate_lsu3shell_relative_operators.
    """

    print("{nuclide}".format(**task))
    model_space_filename = "model_space.dat".format(**task)
    basis_listing_filename = "lsu3shell_basis.dat"

    command_line=[su3basis_executable,model_space_filename,basis_listing_filename]
    mcscript.call(
        command_line,
        mode=mcscript.CallMode.kSerial
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
        "model_space.dat",
        "relative_operators.dat",
        "lsu3shell_basis.dat",
        "relative_unit_tensor_labels.dat"
    ]
    archive_file_list += glob.glob('*.rme')

    # generate archive
    archive_filename = "su3rme-{}.tgz".format(su3rme_descriptor)
    mcscript.call(
        [
            "tar", "-zcvf", archive_filename
        ] + archive_file_list
    )

    # save independent copy of basis listing outside tarball for easy inspection
    basis_filename = "lsu3shell_basis_{}.dat".format(su3rme_descriptor)
    mcscript.call(
        [
            "cp", "lsu3shell_basis.dat", basis_filename
        ]
    )

    # move archive to results directory (if in multi-task run)
    if (mcscript.task.results_dir is not None):
        mcscript.call(
            [
                "mv",
                "--verbose",
                basis_filename,archive_filename,
                "--target-directory={}".format(mcscript.task.results_dir)
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
        "su3rme-{}".format(su3rme_descriptor),
        error_message="Data directory for SU(3) RMEs not found",
        fail_on_not_found=False
    )
    archive_filename = mcscript.utils.search_in_subdirectories(
        su3rme_directory_list,
        su3rme_subdirectory_list,
        "su3rme-{}.tgz".format(su3rme_descriptor),
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

def link_su3rme_files(task):
    """ Search for expanded archive of relative operator SU(3) RME files.

    
    """


def do_generate_lsu3shell_rmes(task):
    """
    Control code for generating RMEs in the SU(3)-NCSM basis, for relative
    unit tensors and symplectic raising/lowering/N operators.
    """

    # retrieve relevant operator files
    retrieve_operator_files(task)

    # generate model space file needed by lsu3shell codes
    generate_model_space_file(task)

    # generate basis listing for basis in which rmes are calculated
    generate_basis_table(task)

    # generate operators rmes
    calculate_rmes(task)

    # save results
    save_su3rme_files(task)

    # clean up working directory
    mcscript.call(["du","-hs","."])  # log working directory disk usage
    delete_filenames=glob.glob('*')
    ##delete_filenames=glob.glob('*.rme')
    ##delete_filenames+=glob.glob('*.PN')
    ##delete_filenames+=glob.glob('*.PPNN')
    mcscript.call(["rm"] + delete_filenames)
    

################################################################
# generate SU(3)-coupled relative matrix elements of observables
################################################################

def generate_observable_rmes(task):
    """Generate relative U3ST RMEs of observable operators.
    
    This may either be by upcoupling relative RMEs or by analytic
    expressions.

    Invokes generate_relative_operator_rmes.

    Output directory:
        relative_observables

    Output filename format:
        {}_hw{:.1f}_Nmax{:02d}_u3st.dat

    """

    mcscript.utils.mkdir("relative_observables")
    os.chdir("relative_observables")
    
    # set parameters
    A = int(task["nuclide"][0]+task["nuclide"][1])
    Nmax=task["Nmax"]
    J0=0
    T0=-1
    g0=0
    J_max_jisp=4
    J_max_coulomb=21

    # generate Hamiltonian RMEs (by upcoupling)
    for hw in mcscript.utils.value_range(*task["hw_range"]):    

        # generate load file
        interaction_filename = mcscript.utils.search_in_subdirectories(
            interaction_directory_list,
            interaction_subdirectory_list,
            task["interaction_filename_template"].format(hw=hw),
            error_message="relative interaction file not found"
        )
        hamiltonian_input_lines = [
            "{}".format(hw),
            "Tintr 1.",
            "INT 1. {} {} {} {} {}".format(J_max_jisp,J0,T0,g0,interaction_filename,**task)
        ]

        if task["use_coulomb"]==True:
            coulomb_filename = mcscript.utils.search_in_subdirectories(
                interaction_directory_list,
                interaction_subdirectory_list,
                task["coulomb_filename"],
                error_message="relative interaction file not found (for Coulomb interaction)"
            )
            hamiltonian_input_lines+=["COUL 1. {} {} {} {} {}".format(J_max_coulomb,J0,T0,g0,coulomb_filename,**task)]
        hamiltonian_load_filename = "hamiltonian.load"
        mcscript.utils.write_input(hamiltonian_load_filename,hamiltonian_input_lines,verbose=True)

        # Call code to upcouple and generate input file for hamiltonian
        #
        # TODO: fix to take (N,Z) instead of (A,N1v), and remove N1v from task dictionary
        command_line = [
                generate_relative_operator_rmes_executable,
                "{}".format(A) ,   
                "{Nmax:d}".format(**task),
                "{N1v:d}".format(**task),
                "hamiltonian"
            ]
        mcscript.call(
            command_line,
            mode=mcscript.CallMode.kSerial
        )

    # generate RMEs for other observables (analytically)
    for hw in mcscript.utils.value_range(*task["hw_range"]):    
        # generate observable load files      
        for observable in task["observables"] :
            observable_name=observable[0]
            # Generate load files for other observables
            input_lines = [
                "{}".format(hw),
                "{} 1.".format(observable_name)
            ]
            load_file_name = "{}.load".format(observable_name)
            mcscript.utils.write_input(load_file_name,input_lines,verbose=True)
            # Generate observable u3st rmes 
            print("made load file")
            command_line = [
                generate_relative_operator_rmes_executable,
                "{:d}".format(A) ,   
                "{Nmax:d}".format(**task),
                "{N1v:d}".format(**task),
                "{}".format(observable_name)
            ]
            mcscript.call(
                command_line,
                mode=mcscript.CallMode.kSerial
            )

    os.chdir("..")

################################################################
# spncci execution
################################################################

def generate_spncci_control_file(task):
    """ Generate control file for spncci run.

    Output file: spncci.dat
    """

    hw_min=task["hw_range"][0]
    hw_max=task["hw_range"][1]
    hw_step=task["hw_range"][2]
    twice_J_min=int(2*task["J_range"][0])
    twice_J_max=int(2*task["J_range"][1])
    J_step=task["J_range"][2]
    J0=0#task["J0"]
    coulomb = int(task["use_coulomb"])


    # write basic parameters
    input_lines = [
        "{nuclide[0]:d} {nuclide[1]:d} {Nsigma_max:d} {Nmax:d}".format(**task),
        "{num_eigenvalues:d}".format(**task),
        "{:d} {:d} {:d}".format(twice_J_min,twice_J_max,J_step),
        "{:f} {:f} {:f}".format(hw_min,hw_max,hw_step),
        "{:s} {:d}".format(coulomb=coulomb,**task)
    ]

    # write observables
    full_observables = [("hamiltonian",0)] + task["observables"]
    for observable in full_observables:
        input_lines.append("{:s} {:d}".format(observable[0],observable[1]))

    control_filename = "spncci.dat"
    mcscript.utils.write_input(control_filename,input_lines,verbose=True)

def call_spncci(task):
    """ Carry out spncci run.
    """

    ## A = int(task["nuclide"][0]+task["nuclide"][1])  # why cast to int???
    ## twice_Nsigma_0 = int(2*task["Nsigma_0"])

    if ("spncci_variant" not in task):
        task["spncci_variant"] = "spncci"
    spncci_executable = os.path.join(spncci_executable_dir,task["spncci_variant"])

    command_line = [
        spncci_executable
    ]
    mcscript.call(
        command_line,
        mode=mcscript.CallMode.kSerial
    )

    # cleanup
    mcscript.call(["rm","-r","lsu3shell_rme","relative_observables"])

def save_spncci_results(task):
    """
    Rename and save spncci results files.
    """

    ## # log file
    ## raw_log_filename = "spncci.out"
    ## new_log_filename = os.path.join(
    ##     mcscript.task.results_dir,
    ##     "{name}-{descriptor}.out".format(name=mcscript.parameters.run.name,**task)
    ## )
    ## mcscript.call(
    ##     [
    ##         "cp",
    ##         "--verbose",
    ##         raw_log_filename,
    ##         new_log_filename
    ##     ]
    ## )

    # results file
    raw_log_filename = "spncci.res"
    new_log_filename = os.path.join(
        mcscript.task.results_dir,
        "{name}-{descriptor}.res".format(name=mcscript.parameters.run.name,**task)
    )
    mcscript.call(
        [
            "cp",
            "--verbose",
            raw_log_filename,
            new_log_filename
        ]
    )

def do_full_spncci_run(task):
    """ Carry out full task of constructing and diagonalizing
    Hamiltonian and other observables.
    """
    retrieve_su3rme_files(task)
    generate_observable_rmes(task)
    generate_spncci_control_file(task)
    call_spncci(task)
    save_spncci_results(task)

if (__name__ == "__MAIN__"):
    pass
