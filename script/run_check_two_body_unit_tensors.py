"""run_check_two_body_unit_tensors.py

  Run check suite on two_body_unit tensor generation.

  Arguments:
    Nstep Nmax

  Note: Requires utils subpackage from mcscript.  This should be
  restructured to load it from within the mcscript package once the
  overhaul of the mcscript package structure is complete.

  Language: Python 3

  M. A. Caprio
  Department of Physics, University of Notre Dame

  9/10/16 (mac): Created.
  11/17/16 (mac): Fix up path construction for executables.

"""

import os
import sys
import utils

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
generate_lsu3shell_two_body_unit_tensors_executable = os.path.join(projects_root,"spncci","programs","unit_tensors","generate_lsu3shell_two_body_unit_tensors")
check_two_body_unit_tensors_executable = os.path.join(projects_root,"spncci","programs","unit_tensors","check_two_body_unit_tensors")

# data files
model_space_filename = "model_space.dat"
two_body_operator_filename = "two_body_operators.dat"

# run parameters
Z = 1  # need Z=1 for pn tbmes
N = 1  # need N=1 for pn tbmes
Nmax = 0  # now read from command line below
Nstep = 0  # now read from command line below
twice_Nsigma_0 = 6  # twice_Nsigma_0=6 for two-body system

################################################################
# create unit tensor operators
################################################################
def generate_basis_table():
    """Create basis table 
    """
    command_line=[su3basis_executable,model_space_filename,"lsu3shell_basis.dat"]
    utils.call(command_line)

def generate_unit_tensors():
    """ Create recoupler input files.
    """

    # call generate_lsu3shell_two_body_unit_tensors
    command_line = [
        generate_lsu3shell_two_body_unit_tensors_executable,
        str(Z),
        str(N),
        str(Nmax),
        str(Nstep)
    ]
    utils.call(command_line)

def read_unit_tensor_list():
    """ Read list of unit tensor basenames.

    Returns:
        (list) : list of unit tensor names
    """

    # read list of unit tensors
    two_body_operator_stream = open(two_body_operator_filename,mode="rt")
    two_body_operator_basename_list = [
        line.strip()  # remove trailing newline
        for line in two_body_operator_stream
    ]
    two_body_operator_stream.close()
    return two_body_operator_basename_list

def recouple_unit_tensors(two_body_operator_basename_list):
    """ Invoke lsu3shell recoupler code on unit tensors.

    Arguments:
        two_body_operator_basename_list (list) : list of unit tensor names
    """

    # iterate over unit tensors
    for basename in two_body_operator_basename_list:

        # call recoupler
        command_line = [
            recoupler_executable,
            "{}.recoupler".format(basename),
            basename
        ]
        utils.call(command_line)

def calculate_rmes(two_body_operator_basename_list):
    """ Invoke lsu3shell SU3RME code to calculate rmes.

    Arguments:
        two_body_operator_basename_list (list) : list of unit tensor names
    """

    # iterate over unit tensors
    for basename in two_body_operator_basename_list:

        # generate load file
        input_lines = [
            "00",   # dummy hw value
            "INT {}".format(basename)
        ]
        load_filename = "{}.load".format(basename)
        utils.write_input(load_filename,input_lines,silent=False)
        print("load_file finished")
        # call SU3RME
        command_line = [
            su3rme_executable,
            model_space_filename,
            model_space_filename,
            load_filename,
            basename
        ]
        utils.call(command_line)

def check_rmes():
    """ Check rmes.
    """

    # call check_two_body_unit_tensors
    command_line = [
        check_two_body_unit_tensors_executable,
        str(Z),
        str(N),
        str(Nmax),
        str(Nstep),
        str(twice_Nsigma_0)
    ]
    utils.call(command_line)



################################################################
# main program
################################################################

# command line parameters
Nmax = int(sys.argv[1])
Nstep = int(sys.argv[2])

generate_unit_tensors()
two_body_operator_basename_list = read_unit_tensor_list()
recouple_unit_tensors(two_body_operator_basename_list)
calculate_rmes(two_body_operator_basename_list)
generate_basis_table()
check_rmes()
