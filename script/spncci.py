"""spncci.py -- define scripting for spncci runs

          
  Language: Python 3

  A. E. McCoy and M. A. Caprio
  Department of Physics, University of Notre Dame

  1/8/17 (aem,mac): Created with code from compute_relative_tensors_lsu3shell_rmes.py.
  2/23/17 (aem,mac): Update rme invocation and add spncci handler.

"""
  
import glob
import os
import mcscript

################################################################
# global configuration
################################################################

# executable files
projects_root = os.path.join(os.environ["HOME"],"projects")
# ... from lsu3shell
recoupler_executable = os.path.join(projects_root,"lsu3shell","programs","upstreams","RecoupleSU3Operator")
su3rme_executable = os.path.join(projects_root,"lsu3shell","programs","tools","SU3RME")
su3basis_executable =os.path.join(projects_root,"lsu3shell","programs","tools","ncsmSU3xSU2IrrepsTabular")
# ... from spncci
generate_lsu3shell_relative_operators_executable = os.path.join(projects_root,"spncci","programs","unit_tensors","generate_lsu3shell_relative_operators")
generate_relative_operator_rmes_executable = os.path.join(projects_root,"spncci","programs","operators","generate_relative_u3st_operators")
spncci_executable = os.path.join(projects_root,"spncci","programs","spncci","spncci")
module_file=os.path.join(projects_root,"spncci","config","module-load-ndcrc.csh")
################################################################
# relative unit tensor evaluation
################################################################


# task parameters
#     nuclide (tuple of int): (N,Z)
#     Nmax (int): oscillator Nmax for many-body basis
#     Nstep (int): step in N for many-body basis (1 or 2)
#     Nsigma_0 (float): U(1) quantum number of lowest configuration
#         (in general can be half integer since contains zero-point offset)
#         (DEPRECATED) 
#     N1v (int): valence shell oscillator N


def generate_relative_operators(task):
    """ Create recoupler input files.
    """
    # call generate_lsu3shell_two_body_unit_tensors
    command_line = [
        generate_lsu3shell_relative_operators_executable,
        "{nuclide[0]:d}".format(**task),
        "{nuclide[1]:d}".format(**task),
        "{Nsigma_ex_max:d}".format(**task),
        "{Nstep:d}".format(**task),
        "{N1v:d}".format(**task)
    ]
    mcscript.call(
        command_line,
        mode=mcscript.call.serial
    )

def generate_basis_table(task):
    """Create basis table.

    Depends on model space file created by generate_lsu3shell_relative_operators.
    """

    print("{nuclide}".format(**task))
    model_space_filename = "model_space_{nuclide[0]:02d}_{nuclide[1]:02d}_Nmax{Nsigma_ex_max:02d}.dat".format(**task)
    basis_listing_filename = "lsu3shell_basis.dat"

    command_line=[su3basis_executable,model_space_filename,basis_listing_filename]
    mcscript.call(
        command_line,
        mode=mcscript.call.serial
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
            mode=mcscript.call.serial
        )

def calculate_rmes(task,relative_operator_basename_list):
    """ Invoke lsu3shell SU3RME code to calculate rmes.

    Arguments:
        relative_operator_basename_list (list) : list of operator file names
    """

    model_space_filename = "model_space_{nuclide[0]:02d}_{nuclide[1]:02d}_Nmax{Nsigma_ex_max:02d}.dat".format(**task)

    # iterate over unit tensors
    for basename in relative_operator_basename_list:

        # generate load file
        input_lines = [
            "00",   # dummy hw value
            "INT {}".format(basename)
        ]
        load_filename = "{}.load".format(basename)
        rme_filename="{}.rme".format(basename)
        mcscript.utils.write_input(load_filename,input_lines,verbose=True)
        print("load_file finished")
        # call SU3RME
        command_line = [
            su3rme_executable,
            model_space_filename,
            model_space_filename,
            load_filename,
            rme_filename
        ]
        mcscript.call(
            command_line,
            mode=mcscript.call.serial
        )

def generate_lsu3shell_rmes(task):
    """ Carry out full task of generating set of relative tensor rmes.
    """

    ## print(task)
    mcscript.utils.mkdir("lsu3shell_rme")
    os.chdir("lsu3shell_rme")
    generate_relative_operators(task)
    generate_basis_table(task)
    relative_operator_basename_list = read_unit_tensor_list(task)
    recouple_operators(task,relative_operator_basename_list)
    calculate_rmes(task,relative_operator_basename_list)
    os.chdir("..")


def generate_interaction_rmes(task):
    """ Generate u3st observable rmes.
            {}_hw{:.1f}_Nmax{:02d}_u3st.dat
    """
    mcscript.utils.mkdir("relative_observables")
    os.chdir("relative_observables")
    J0=0
    T0=0
    g0=0
    J_max=4
    A = int(task["nuclide"][0]+task["nuclide"][1])
    Nmax=task["Nmax"]
    # generate hamiltonian load file
    for hw in mcscript.utils.value_range(10,30,2.5):
        interaction_filename=task["interaction_filename_template"].format(hw)
        hamiltonian_input_lines = [
            "{}".format(hw),
            "Tintr 1.",
            "INT 1. {} {} {} {} {}".format(J_max,J0,T0,g0,interaction_filename,**task)
        ]
        hamiltonian_load_filename = "hamiltonian.load"
        mcscript.utils.write_input(hamiltonian_load_filename,hamiltonian_input_lines,verbose=True)

        command_line = [
                generate_relative_operator_rmes_executable,
                "{}".format(A) ,   
                "{Nmax:d}".format(**task),
                "{N1v:d}".format(**task),
                "hamiltonian"
            ]

        mcscript.call(
            command_line,
            mode=mcscript.call.serial
        )
      
        for observable in task["observables"] :
            # Generate load files for other observables
            input_lines = [
                "{}".format(hw),
                "{} 1.".format(observable)
            ]
            load_file_name = "{}.load".format(observable)
            mcscript.utils.write_input(load_file_name,input_lines,verbose=True)
            # Generate observable u3st rmes 
            print("made load file")
            command_line = [
                generate_relative_operator_rmes_executable,
                "{:d}".format(A) ,   
                "{Nmax:d}".format(**task),
                "{N1v:d}".format(**task),
                "{}".format(observable)
            ]
            mcscript.call(
                command_line,
                mode=mcscript.call.serial
            )

    os.chdir("..")

def generate_spncci_control_file(task):
    """ control file for spncci observables 
    """
    hw_min=task["hw_range"][0]
    hw_max=task["hw_range"][1]
    hw_step=task["hw_range"][2]
    twice_J_min=2*task["J_range"][0]
    twice_J_max=2*task["J_range"][1]
    J_step=task["J_range"][2]
    J0=task["J0"]

    input_lines = [
            "{} {} {} {}".format(J0,twice_J_min,twice_J_max,J_step),
            "{} {} {}".format(hw_min,hw_max,hw_step),   # dummy hw value
            "hamiltonian"
        ]
    for observable in task["observables"]:
        input_lines.append(observable)

    load_filename = "spncci.load"
    mcscript.utils.write_input(load_filename,input_lines,verbose=True)

def call_spncci(task):
    """ compute observable rmes in spncci basis.
    """
    A = int(task["nuclide"][0]+task["nuclide"][1])
    twice_Nsigma_0 = int(2*task["Nsigma_0"])

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
        mode=mcscript.call.serial
    )

def do_full_spncci_run(task):
    """ Carry out full task of constructing and diagonalizing hamiltonian and other observables.
    """
    unit_tensor_directory=task["unit_tensor_directory"].format(**task)
    # generate_lsu3shell_rmes(task)
    os.symlink(unit_tensor_directory, "lsu3shell_rme")
    generate_interaction_rmes(task)
    generate_spncci_control_file(task)
    call_spncci(task)


if (__name__ == "__MAIN__"):
    pass
