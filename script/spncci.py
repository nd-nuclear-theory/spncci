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
            setenv SPNCCI_SEED_DIR ${HOME}/data/spncci/seeds:/afs/crc.nd.edu/group/nuclthy/data/spncci/seeds
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
  2/15/18 (aem): Add seed calculation scripting
  2/16/18 (aem): Factored out su3rme computation into su3rme.py
"""
  
import glob
import mcscript
import os
import su3rme
################################################################
# global configuration
################################################################

# environment configuration variables
project_root = os.environ["SPNCCI_PROJECT_ROOT_DIR"]
interaction_directory_list = os.environ["SPNCCI_INTERACTION_DIR"].split(":")
interaction_subdirectory_list = []
operator_directory_list = os.environ["SPNCCI_OPERATOR_DIR"].split(":")
operator_subdirectory_list = []

# seed_directory_list = os.environ["SPNCCI_SEED_DIR"].split(":")
seed_directory_list = os.environ["SPNCCI_SEED_DIR"]
seed_subdirectory_list = []

truncation_directory = os.environ["SPNCCI_TRUNCATION_DIR"]
truncation_subdirectory = []


# ... from spncci
generate_relative_operator_rmes_executable = os.path.join(project_root,"spncci","programs","operators","generate_relative_u3st_operators")
generate_spncci_seed_files_executable = os.path.join(project_root,"spncci","programs","lgi","get_spncci_seed_blocks")
spncci_executable_dir = os.path.join(project_root,"spncci","programs","spncci")
seed_descriptor_template_Nsigmamax = "Z{nuclide[0]:02d}-N{nuclide[1]:02d}-Nsigmamax{Nsigma_max:02d}-Nstep{Nstep:d}"
################################################################
# generate SU(3)-coupled relative matrix elements of observables
################################################################
def stacksize_setup():
    """ Set OpenMP stacksize variables.

    """
    stacksize=50 #in megabytes
    # set number of threads by global qsubm depth parameter
    print("Setting OMP_STACKSIZE to {}M.".format(stacksize))
    os.environ["OMP_STACKSIZE"] = "20 M"




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
# seed execution
################################################################
def generate_spncci_seed_files(task):
    """
    Generates:
   
    "seeds/lgi_families.dat": File containing a list of lgi family labels with 
        their corresponding index 
            i  twice_Nsigma lambda_sigma mu_sigma twice_Sp twice_Sn twice_S

    For each pair of lgi families (i,j) with i>=j, there are two
    files in directory seeds containing:
        1. "seeds/seeds_00000i_00000j.rmes": All non-zero unit tensor seed blocks
          
            binary_format_code (1) and binary precision (4 or 8)
            num_rows, num_cols, rmes...

        2. "seeds/operators_00000i_00000j.dat": unit tensor labels corresponding
            to each seed block
          
            lambda0, mu0, 2S0, 2T0, etap, 2Sp, 2Tp, eta, 2S, 2T, rho0
             
    "seeds/lgi_expansions.dat": File containing expansion of lgis in lsu3shell basis.
            For each lsu3shell subspace, 
                rows, cols, rmes

    """
    # if seed directory already exist, remove and recreate fresh copy
    if (os.path.exists("seeds")):
        mcscript.call(["rm", "-r","seeds"])

    mcscript.utils.mkdir("seeds")


    # generate seed rmes
    command_line = [
        generate_spncci_seed_files_executable,
        "{nuclide[0]:d}".format(**task),
        " {nuclide[1]:d}".format(**task),
        " {Nsigma_max:d}".format(**task)
    ]
    mcscript.call(
        command_line,
        mode=mcscript.CallMode.kSerial
    )
    
    # mcscript.call(["ls"])
    # change seed directory name
    seed_descriptor = task["seed_descriptor_template"].format(**task)
    seed_directory = "seeds-{}".format(seed_descriptor)
    
    # if seed directory already exist, remove and recreate fresh copy
    if (os.path.exists(seed_directory)):
        mcscript.call(["rm", "-r",seed_directory])
    
    mcscript.call(["mv","seeds",seed_directory])


def save_seed_files(task):
    """Create archive of seed RMEs 

    Some auxiliary files (e.g., the list of operators) are saved as well.

    Manual follow-up: The rme files are bundled into tgz files and saved to
    the current run's results directory.  They should then be moved
    (manually) to the directory specified in SPNCCI_SEED_DIR, for
    subsequent use.

    """
    seed_descriptor = task["seed_descriptor_template"].format(**task)

    # generate archive
    directory="seeds-{}".format(seed_descriptor)
    archive_filename = "seeds-{}.tgz".format(seed_descriptor)


    mcscript.call(
        [
            "tar", "-zcvf", archive_filename, "--format=posix", directory
        ] 
    )

    #TODO: change to create ln to directory in results directory. 

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


    # clean up working directory
    # mcscript.call(["rm","-r", "seeds"])
    mcscript.call(["rm","-r", "lsu3shell_rme"])


def retrieve_seed_files(task):
    """ Retrieve archive of seed RME files.

    Files are retrieved into a subdirectory named seeds.
    """
    # identify seed data directory
    seed_descriptor = task["seed_descriptor_template"].format(**task)
    directory_name = mcscript.utils.search_in_subdirectories(
        seed_directory_list,
        seed_subdirectory_list,
        "seeds-{}".format(seed_descriptor),
        error_message="Data directory for seed RMEs not found",
        fail_on_not_found=False
    )
    archive_filename = mcscript.utils.search_in_subdirectories(
        seed_directory_list,
        seed_subdirectory_list,
        "seeds-{}.tgz".format(seed_descriptor),
        error_message="Archive file for seed RMEs not found",
        fail_on_not_found=False
    )

    # remove any existing symlink or data directory
    if (os.path.exists("seeds")):
        mcscript.call(["rm","-r","seeds"])
        # print("hi")
    

    if (directory_name is not None):
        print("seed directory ",directory_name)
        mcscript.call(["ln", "-s", directory_name, "seeds"])

    elif (archive_filename is not None):
        # extract archive contents
        directory_name="seeds-{}".format(seed_descriptor)
        mcscript.call(["tar","-xf",archive_filename])
        mcscript.call(["ln", "-s", directory_name, "seeds"])


    else:
        raise(mcscript.exception.ScriptError("Cannot find seed RME data"))


    # link to data seed directory
    


################################################################
# setting up lgis
################################################################

def get_lgi_file(task):
    """
    Creates symbolic link from list of lgi families to be included in basis
    to "lgi_families.dat". If no truncation file is given, symbolic link to 
    list of full space of lgi families "seeds/lgi_families.dat". 
    """
    # if link already exists, remove 

    print("lgi_families.dat exits ",os.path.exists("lgi_families.dat"))
    if (os.path.exists("lgi_families.dat")):
        mcscript.call(["rm","-r","lgi_families.dat"])
        print("removed lgi families file")

    # if no truncation filename is given, then use full Nsigma,max space
    if task["truncation_filename"]==None:
        mcscript.call(
            [
                "cp",
                "seeds/lgi_families.dat",
                "lgi_families.dat"
            ]
        )
    
    # create symbolic link to truncated list of lgi family labels
    else :
        mcscript.call(
            [
                "ln",
                "-s",
                task["truncation_filename"],
                "lgi_families.dat"
            ]
        )
            
def read_lgi_list(filename):
    """
    Read in lgi family labels from filename
    """
    lines = [line.rstrip('\n') for line in open(filename,'r')]
    for line in lines:
        lgi_labels=[[int(x) for x in line.split()][:-1] for line in lines]

    return lgi_labels

def lookup_table(lgi_labels_truncated,lgi_labels):
    """
    Create look up table between lgi_family_index in basis and 
    lgi_family_index in full space, which indexes the seed files 
    """
    index_lookup=[]
    for label in lgi_labels_truncated:
        index=lgi_labels.index(label)
        index_lookup.append(index)

    return index_lookup


def write_lookup_table(index_lookup):
    """
    Create file containing lookup table
    """
    filename="seeds/lgi_full_space_lookup_table.dat"
    with open(filename, 'w') as outstream:
        for index in range(len(index_lookup)):
            full_space_index=index_lookup[index]
            outstream.write("{} {}\n".format(index,full_space_index))

    outstream.close()


def generate_lgi_lookup_table(task):
    """
    define lgi file containing lgi families
    create look up table between lgi family indices in 
    (possibly) truncated space and full space
    """
    get_lgi_file(task)

    # Get LGI labels of full space
    lgi_labels=read_lgi_list("seeds/lgi_families.dat")

    # print(lgi_labels)
    # Get LGI labels of truncated space
    lgi_labels_truncated=read_lgi_list("lgi_families.dat")
    # print(lgi_labels_truncated)
    
    # Make look up table for related in truncated space index to full space index
    index_lookup=lookup_table(lgi_labels_truncated,lgi_labels)

    #write table to file
    write_lookup_table(index_lookup)



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

    if(task["truncation_filename"]==None):
        transform_lgi=0
    elif(task["transformation_filename"]==None):
        transform_lgi=0
    else:
        transform_lgi=1

        if (os.path.exists("lgi_transformations.dat")):
            mcscript.call(["rm","-r","lgi_transformations.dat"])

        mcscript.call(
            [
                "ln",
                "-s",
                task["transformation_filename"],
                "lgi_transformations.dat"
            ]
        )

    # write basic parameters
    input_lines = [
        "{nuclide[0]:d} {nuclide[1]:d} {Nsigma_max:d} {Nmax:d} {transform_lgi:d}".format(transform_lgi=transform_lgi,**task),
        "{num_eigenvalues:d} {eigensolver_num_convergence:d} {eigensolver_max_iterations:d} {eigensolver_tolerance:.2e}".format(**task),
        "{:d} {:d} {:d}".format(twice_J_min,twice_J_max,J_step),
        "{:f} {:f} {:f}".format(hw_min,hw_max,hw_step),
        "{interaction:s} {coulomb:d}".format(coulomb=coulomb,**task)
    ]

    # write observables
    full_observables = [("hamiltonian",0)] + task["observables"]
    for observable in full_observables:
        input_lines.append("{:s} {:d}".format(observable[0],observable[1]))

    control_filename = "spncci.dat"
    mcscript.utils.write_input(control_filename,input_lines,verbose=True)

def make_hyperblocks_dir(task):
    """
    Defines temporary director for hyperblocks 
    """
    if (os.path.exists("hyperblocks")):
        mcscript.call(["rm","-r","hyperblocks"])


    if task["hyperblocks_dir"]==None:
        mcscript.utils.mkdir("hyperblocks")

    else:
        directory=task["hyperblocks_dir"]
        mcscript.utils.mkdir(directory)

        # link to hyperblocks temporary directory
        mcscript.call(
            [
                "ln",
                "-s",
                directory,
                "hyperblocks"
            ]
        )


def call_spncci(task):
    """ Carry out spncci run.
    """

    if ("spncci_variant" not in task):
        task["spncci_variant"] = "spncci"
    spncci_executable = os.path.join(spncci_executable_dir,task["spncci_variant"])

    command_line = [spncci_executable]
    mcscript.call(
        command_line,
        mode=mcscript.CallMode.kSerial
    )

    # cleanup
    mcscript.call(["rm","-r","seeds","relative_observables","hyperblocks"])

    seed_descriptor = task["seed_descriptor_template"].format(**task)
    seed_directory = "seeds-{}".format(seed_descriptor)

    if (os.path.exists(seed_directory)):
        mcscript.call(["rm","-r",seed_directory])

def save_spncci_results(task):
    """
    Rename and save spncci results files.
    """
    # results file
    coulomb = int(task["use_coulomb"])
    results_descriptor="Z{nuclide[0]:d}-N{nuclide[1]:d}-{interaction:s}-{coulomb:1d}-Nsigmamax{Nsigma_max:02d}-Nmax{Nmax:02d}".format(coulomb=coulomb,**task)
    raw_log_filename = "spncci.res"
    new_log_filename = os.path.join(
        mcscript.task.results_dir,
        "{name}-{description}.res".format(name=mcscript.parameters.run.name,description=results_descriptor)
    )
    mcscript.call(
        [
            "cp",
            "--verbose",
            raw_log_filename,
            new_log_filename
        ]
    )

def do_seed_run(task):
    """ Generate seed files from lsu3shell files.
    """
    su3rme.retrieve_su3rme_files(task)
    generate_spncci_seed_files(task)
    save_seed_files(task)

def do_full_spncci_run(task):
    """ Carry out full task of constructing and diagonalizing
    Hamiltonian and other observables.
    """
    # stacksize_setup()
    retrieve_seed_files(task)
    generate_lgi_lookup_table(task)
    generate_observable_rmes(task)
    generate_spncci_control_file(task)
    make_hyperblocks_dir(task)
    call_spncci(task)
    save_spncci_results(task)

def test(task):
    """
    Read in lgi family labels from filename
    """
    make_hyperblocks_dir(task)
    filename="hyperblocks/testfile.txt"
    file = open(filename,'w') 
    file.write("Hello World") 
    file.write("This is our new text file") 
    file.write("and this is another line.") 
    file.write("Why? Because we can.") 
    print("done writing")
    file.close() 
    
    file = open(filename, 'r') 
    print(file.read())
    mcscript.call(["rm","-r","hyperblocks"])

if (__name__ == "__MAIN__"):
    pass
