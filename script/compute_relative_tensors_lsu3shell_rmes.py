"""compute_relative_tensors_lsu3shell_rmes.py

  Calculate each of the relative unit tensors and 
  Brel and Nintr used for obtaining LGI expansion

  Arguments:
    Z N twice_Nsigma_0 Nmax Nstep N1B

  Note: Requires utils subpackage from mcscript.  This should be
  restructured to load it from within the mcscript package once the
  overhaul of the mcscript package structure is complete.

  Language: Python 3

  A. E. McCoy
  Department of Physics, University of Notre Dame

  12/1/16 (aem): Created.
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
su3rme_executable = os.path.join(projects_root,"lsu3shell","programs","tools","SU3RME_MPI")#changed name
su3basis_executable =os.path.join(projects_root,"lsu3shell","programs","tools","ncsmSU3xSU2IrrepsTabular")
# ... from spncci
generate_lsu3shell_relative_operators_executable = os.path.join(projects_root,"spncci","programs","unit_tensors","generate_lsu3shell_relative_operators")

# run parameters
Z = 1  # now read from command line below
N = 1  # now read from command line below
Nmax = 0  # now read from command line below
Nstep = 0  # now read from command line below
twice_Nsigma_0 = 6  #now read from command line below

################################################################
# create unit tensor operators
################################################################
def generate_basis_table():
    """Create basis table 
    """
    command_line=[su3basis_executable,model_space_filename,"lsu3shell_basis.dat"]
    utils.call(command_line)

def generate_relative_operators():
    """ Create recoupler input files.
    """

    # call generate_lsu3shell_two_body_unit_tensors
    command_line = [
        generate_lsu3shell_relative_operators_executable,
        str(Z),
        str(N),
        str(Nmax),
        str(Nstep),
        str(N1B)
    ]
    utils.call(command_line)

def read_unit_tensor_list():
    """ Read list of unit tensor basenames.

    Returns:
        (list) : list of unit tensor names
    """

    # read list of unit tensors
    relative_operator_stream = open(relative_operator_filename,mode="rt")
    relative_operator_basename_list = [
        line.strip()  # remove trailing newline
        for line in relative_operator_stream
    ]
    relative_operator_stream.close()
    return relative_operator_basename_list

def recouple_operators(relative_operator_basename_list):
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
        utils.call(command_line)

def calculate_rmes(relative_operator_basename_list):
    """ Invoke lsu3shell SU3RME code to calculate rmes.

    Arguments:
        relative_operator_basename_list (list) : list of operator file names
    """

    # iterate over unit tensors
    for basename in relative_operator_basename_list:

        # generate load file
        input_lines = [
            "00",   # dummy hw value
            "INT {}".format(basename)
        ]
        load_filename = "{}.load".format(basename)
        rme_filename="{}.rme".format(basename)
        utils.write_input(load_filename,input_lines,silent=False)
        print("load_file finished")
        # call SU3RME
        command_line = [
            su3rme_executable,
            model_space_filename,
            model_space_filename,
            load_filename,
            rme_filename
        ]
        utils.call(command_line, mode=mcscript.call.hybrid)

################################################################
# main program
################################################################

# command line parameters
Z=int(sys.argv[1])
N=int(sys.argv[2])
twice_Nsigma_0=int(sys.argv[3])
Nmax = int(sys.argv[4])
Nstep = int(sys.argv[5])
N1B = int(sys.argv[6])

# data files
model_space_filename = "model_space_{:02d}_{:02d}_Nmax{:02d}.dat".format(Z, N,Nmax)
relative_operator_filename = "relative_operators.dat"

generate_relative_operators()
relative_operator_basename_list = read_unit_tensor_list()
recouple_operators(relative_operator_basename_list)
calculate_rmes(relative_operator_basename_list)
generate_basis_table()
