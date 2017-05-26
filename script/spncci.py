"""spncci.py -- define scripting for spncci runs


    Task parameters:

        # basic space parameters -- for lsu3shell rme evaluation and spncci:

        nuclide (tuple of int): (N,Z)
        Nmax (int): oscillator Nmax for many-body basis
        Nstep (int): step in N for many-body basis (1 or 2)
        N1v (int): valence shell oscillator N
        Nsigma_0 (float): U(1) quantum number of lowest configuration
            (in general can be half integer since contains zero-point offset)
        Nsigma_ex_max (int): maximum oscillator exitation for LGIs

        # unit tensor parameters -- for lsu3shell rme evaluation ONLY
        J0 (int): restriction on J0 for unit tensors (normally -1 to include
            all J0 as needed for spncci recurrence)
        "unit_tensor_filename_template" (str): template for unit tensor output filenames

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

"""
  
import glob
import mcscript
import os

################################################################
# global configuration
################################################################

# executable files
projects_root = os.path.join(os.environ["HOME"],"projects")
# ... from lsu3shell
recoupler_executable = os.path.join(projects_root,"lsu3shell","programs","upstreams","RecoupleSU3Operator")
su3rme_executable = os.path.join(projects_root,"lsu3shell","programs","tools","SU3RME_MPI")
su3basis_executable = os.path.join(projects_root,"lsu3shell","programs","tools","ncsmSU3xSU2IrrepsTabular")
# ... from spncci
generate_lsu3shell_relative_operators_executable = os.path.join(projects_root,"spncci","programs","unit_tensors","generate_lsu3shell_relative_operators")
generate_relative_operator_rmes_executable = os.path.join(projects_root,"spncci","programs","operators","generate_relative_u3st_operators")
spncci_executable_dir = os.path.join(projects_root,"spncci","programs","spncci")


################################################################
# relative unit tensor evaluation
################################################################

def generate_relative_operators(task):
    """Create recoupler input files for relative unit
    tensors and symplectic raising/lowering/N operators.

    Invokes generate_lsu3shell_relative_operators.
    """

    command_line = [
        generate_lsu3shell_relative_operators_executable,
        "{nuclide[0]:d}".format(**task),
        "{nuclide[1]:d}".format(**task),
        "{Nsigma_ex_max:d}".format(**task),
        "{Nstep:d}".format(**task),
        "{N1v:d}".format(**task),
        "-1",# All J0
        "-1"# All T0
    ]
    mcscript.call(
        command_line,
        mode=mcscript.CallMode.kSerial
    )


def read_unit_tensor_list(task):
    """ Read list of unit tensor basenames.

    Returns:
        (list) : list of unit tensor names
    """

    relative_operator_filename = "relative_operators.dat"

    # read list of unit tensors
    relative_operator_stream = open(relative_operator_filename,mode="rt")
    relative_operator_basename_list = [
        line.strip()  # remove trailing newline
        for line in relative_operator_stream
    ]
    relative_operator_stream.close()
    return relative_operator_basename_list

def recouple_operators(task,relative_operator_basename_list):
    """ Invoke lsu3shell recoupler code on operators.

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

def calculate_rmes(task,relative_operator_basename_list):
    """ Invoke lsu3shell SU3RME code to calculate rmes.

    Invokes SU3RME_MPI.

    Arguments:
        relative_operator_basename_list (list) : list of operator file names
    """

    model_space_filename = "model_space_{nuclide[0]:02d}_{nuclide[1]:02d}_Nmax{Nsigma_ex_max:02d}.dat".format(**task)

    # iterate over unit tensors
    # for basename in relative_operator_basename_list:

        # # generate load file
        # input_lines = [
        #     "00",   # dummy hw value
        #     "INT {}".format(basename)
        # ]
        # load_filename = "{}.load".format(basename)
        # # rme_filename="{}.rme".format(basename)
        # mcscript.utils.write_input(load_filename,input_lines,verbose=True)
    
    # print("load_files finished")
    # call SU3RME
    command_line = [
        su3rme_executable,
        model_space_filename,
        model_space_filename,
        "relative_operators.dat"
    ]
    mcscript.call(
        command_line,
        mode=mcscript.CallMode.kHybrid
    )

def generate_lsu3shell_rmes(task):
    """Generate SU(3)-NCSM RMEs for relative unit tensors and symplectic
    raising/lowering/N operators.

    Output directory: lsu3shell_rme
    """
    
    # set up data directory
    mcscript.utils.mkdir("lsu3shell_rme")
    os.chdir("lsu3shell_rme")

    # generate operators and their rmes
    generate_relative_operators(task)
    relative_operator_basename_list = read_unit_tensor_list(task)
    recouple_operators(task,relative_operator_basename_list)
    calculate_rmes(task,relative_operator_basename_list)

    # do cleanup
    delete_filenames=glob.glob('*.recoupler')
    delete_filenames+=glob.glob('*.PN')
    delete_filenames+=glob.glob('*.PPNN')
    delete_filenames+=glob.glob('*.load')
    mcscript.call(["rm"] + delete_filenames)

    # generate basis listing for basis in which rmes were calculated
    generate_basis_table(task)

    # restore working directory
    os.chdir("..")  

def generate_basis_table(task):
    """Create SU(3)-NCSM basis table.

    Invokes ncsmSU3xSU2IrrepsTabular.

    Depends on model space file created by generate_lsu3shell_relative_operators.
    """

    print("{nuclide}".format(**task))
    model_space_filename = "model_space_{nuclide[0]:02d}_{nuclide[1]:02d}_Nmax{Nsigma_ex_max:02d}.dat".format(**task)
    basis_listing_filename = "lsu3shell_basis.dat"

    command_line=[su3basis_executable,model_space_filename,basis_listing_filename]
    mcscript.call(
        command_line,
        mode=mcscript.CallMode.kSerial
    )

################################################################
# relative matrix element evalation
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
    
    A = int(task["nuclide"][0]+task["nuclide"][1])
    Nmax=task["Nmax"]
    J0=0
    T0=-1
    g0=0
    J_max_jisp=4
    J_max_coulomb=21
    coulomb_filename=os.path.join(task["interaction_directory"],"coulomb_Nmax20_rel.dat")
    # generate hamiltonian load file
    for hw in mcscript.utils.value_range(10,30,2.5):
        interaction_filename=task["interaction_filename_template"].format(hw)
        hamiltonian_input_lines = [
            "{}".format(hw),
            "Tintr 1.",
            "INT 1. {} {} {} {} {}".format(J_max_jisp,J0,T0,g0,interaction_filename,**task)
        ]

        if task["use_coulomb"]==True:
            hamiltonian_input_lines+=["INT 1. {} {} {} {} {}".format(J_max_coulomb,J0,T0,g0,coulomb_filename,**task)]


        hamiltonian_load_filename = "hamiltonian.load"
        mcscript.utils.write_input(hamiltonian_load_filename,hamiltonian_input_lines,verbose=True)

        # Call code to upcouple and generate input file for hamiltonian 
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
    #  Generate SU(3) rme files for each of the observables 
    for hw in mcscript.utils.value_range(10,30,2.5):    
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
# spncii execution
################################################################

def generate_spncci_control_file(task):
    """ Generate control file for spncci run.

    Output file: spncci.load
    """
    hw_min=task["hw_range"][0]
    hw_max=task["hw_range"][1]
    hw_step=task["hw_range"][2]
    twice_J_min=2*task["J_range"][0]
    twice_J_max=2*task["J_range"][1]
    J_step=task["J_range"][2]
    J0=0#task["J0"]

    input_lines = [
            "{} {} {}".format(twice_J_min,twice_J_max,J_step),
            "{} {} {}".format(hw_min,hw_max,hw_step),   # dummy hw value
            "hamiltonian {}".format(J0)
        ]
    for observable in task["observables"]:
        input_lines.append("{} {}".format(observable[0],observable[1]))

    load_filename = "spncci.load"
    mcscript.utils.write_input(load_filename,input_lines,verbose=True)

def call_spncci(task):
    """ Carry out spncci run.
    """
    A = int(task["nuclide"][0]+task["nuclide"][1])
    twice_Nsigma_0 = int(2*task["Nsigma_0"])

    if ("spncci_variant" not in task):
        task["spncci_variant"] = "spncci"
    spncci_executable = os.path.join(spncci_executable_dir,task["spncci_variant"])

    command_line = [
        spncci_executable,
        # TODO determine actual arguments or move into a control file
        "{A:d}".format(A=A,**task) ,   
        "{twice_Nsigma_0:d}".format(twice_Nsigma_0=twice_Nsigma_0,**task),
        "{Nsigma_ex_max:d}".format(**task),
        "{N1v:d}".format(**task),
        "{Nmax:d}".format(**task),
        "{num_eigenvalues:d}".format(**task),
        "spncci"
    ]
    mcscript.call(
        command_line,
        mode=mcscript.CallMode.kSerial
    )

def save_spncci_results(task):
    """
    Ad hoc...
    """

    # results file
    raw_results_filename = "eigenvalues_Nmax{Nmax:02d}_Nsigma_ex{Nsigma_ex_max:02d}.dat".format(**task)
    new_results_filename = os.path.join(mcscript.task.results_dir,"{name}-{descriptor}.dat".format(name=mcscript.parameters.run.name,**task))
    mcscript.call(
        [
            "cp",
            "--verbose",
            raw_results_filename,
            new_results_filename
        ]
    )

    # log file
    raw_log_filename = "spncci.out"
    new_log_filename = os.path.join(mcscript.task.results_dir,"{name}-{descriptor}.out".format(name=mcscript.parameters.run.name,**task))
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
    if ("unit_tensor_directory" in task):
        # LEGACY support for scripts with messed up name
        unit_tensor_directory=task["unit_tensor_directory"].format(**task)
    else:
        unit_tensor_directory=task["unit_tensor_directory_template"].format(**task)
    # generate_lsu3shell_rmes(task)
    ## os.symlink(unit_tensor_directory, "lsu3shell_rme")  # hard to reliably clean up?
    ## mcscript.call(["cp","--recursive",unit_tensor_directory,"lsu3shell_rme"])
    unit_tensor_directory_archive_filename = "{}.tgz".format(unit_tensor_directory)
    mcscript.call(["tar","xf",unit_tensor_directory_archive_filename])
    generate_observable_rmes(task)
    generate_spncci_control_file(task)
    call_spncci(task)
    save_spncci_results(task)

if (__name__ == "__MAIN__"):
    pass
