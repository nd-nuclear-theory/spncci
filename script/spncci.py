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
spncci_executable = os.path.join(projects_root,"spncci","programs","spncci","spncci")

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
        "{Nmax:d}".format(**task),
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
    model_space_filename = "model_space_{nuclide[0]:02d}_{nuclide[1]:02d}_Nmax{Nmax:02d}.dat".format(**task)
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

    model_space_filename = "model_space_{nuclide[0]:02d}_{nuclide[1]:02d}_Nmax{Nmax:02d}.dat".format(**task)

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

def call_spncci(task):
    """ Carry out full task of generating set of relative tensor rmes.
    """

    twice_Nsigma_0 = int(2*task["Nsigma_0"])
    command_line = [
        spncci_executable,
        # TODO determine actual arguments or move into a control file
        "{Nmax:d}".format(**task),
        "{twice_Nsigma_0:d}".format(twice_Nsigma_0=twice_Nsigma_0,**task)
    ]
    mcscript.call(
        command_line,
        mode=mcscript.call.serial
    )

    


if (__name__ == "__MAIN__"):
    pass
