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
import su3rme
################################################################
# global configuration
################################################################

# environment configuration variables
project_root = os.environ["SPNCCI_PROJECT_ROOT_DIR"]
interaction_directory_list = os.environ["SPNCCI_INTERACTION_DIR"].split(":")
interaction_subdirectory_list = []
su3rme_directory_list = os.environ["SPNCCI_SU3RME_DIR"].split(":")
su3rme_subdirectory_list = []

# executable files
generate_relative_operator_rmes_executable = os.path.join(project_root,"spncci","programs","operators","generate_relative_u3st_operators")
generate_spncci_seed_files_executable = os.path.join(project_root,"spncci","programs","lgi","get_spncci_seed_blocks")
spncci_executable_dir = os.path.join(project_root,"spncci","programs","spncci")
    

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
# seed execution
################################################################
def generate_spncci_seed_files(task):
    """
    Generates list of lgi families and recurrence seeds
    """
        # set up data directory
    if (not os.path.exists("seeds")):
        mcscript.utils.mkdir("seeds")

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

def read_lgi_list(filename):
    """
    """
    lines = [line.rstrip('\n') for line in open(filename,'r')]
    for line in lines:
        lgi_labels=[[int(x) for x in line.split()] for line in lines]

    return lgi_labels

def lookup_table(lgi_labels_truncated,lgi_labels):
    """
    """
    index_lookup=[]
    for label in lgi_labels_truncated:
        index=lgi_labels.index(label)
        index_lookup.append(index)

    return index_lookup


def get_lgi_file(task):
    """
    link lgi_families.dat to list of lgi families
        if none, link to list of full space of lgi families found in seeds. 
    """
    if (os.path.exists("lgi_families.dat")):
        mcscript.call(["rm","-r","lgi_families.dat"])


    if task["truncation_filename"]==None:
        mcscript.call(
            [
                "ln",
                "-s",
                "seeds/lgi_families.dat",
                "lgi_families.dat"
            ]
        )
    
    else :
        mcscript.call(
            [
                "ln",
                "-s",
                task["truncation_filename"],
                "lgi_families.dat"
            ]
        )
            

def write_lookup_table(index_lookup):
    filename="lgi_full_space_lookup_table.dat"
    with open(filename, 'w') as outstream:
        for index in range(len(index_lookup)):
            full_space_index=index_lookup[index]
            outstream.write("{} {}\n".format(index,full_space_index))

    outstream.close()


def generate_lgi_lookup_table(tasks):
    # Get LGI labels of full space
    lgi_labels=read_lgi_list("seeds/lgi_families.dat")

    # Get LGI labels of truncated space
    lgi_labels_truncated=read_lgi_list("lgi_families.dat")

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

    # write basic parameters
    input_lines = [
        "{nuclide[0]:d} {nuclide[1]:d} {Nsigma_max:d} {Nmax:d}".format(**task),
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

def do_seed_run(task):
    """ Carry out full task of constructing and diagonalizing
    Hamiltonian and other observables.
    """
    su3rme.retrieve_su3rme_files(task)
    generate_spncci_seed_files(task)
    get_lgi_file(task)
    generate_lgi_lookup_table(task)


def do_full_spncci_run(task):
    """ Carry out full task of constructing and diagonalizing
    Hamiltonian and other observables.
    """
    su3rme.retrieve_su3rme_files(task)
    generate_spncci_seed_files(task)
    get_lgi_file(task)
    generate_lgi_lookup_table(task)
    generate_observable_rmes(task)
    generate_spncci_control_file(task)
    call_spncci(task)
    save_spncci_results(task)

if (__name__ == "__MAIN__"):
    pass
