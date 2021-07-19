"""su3rme.py -- define scripting for lsu3shell generation of su3rmes

    Environment variables:

        SPNCCI_PROJECT_ROOT_DIR -- root directory under which
            lsu3shell and spncci codes are found

        SPNCCI_OPERATOR_DIR -- base directory for relative operator
            files (colon-delimited search path); operator files will
            be sought within subdirectories named in
            spncci.operator_subdirectory_list

        Ex (NERSC m2032):
            export SPNCCI_PROJECT_ROOT_DIR="${HOME}/codes"
            export SPNCCI_OPERATOR_DIR=/project/projectdirs/m2032/data/spncci/operator

    You will also need the directory containing the present script
    file (lsu3shell.py) to be in your Python path:

        export PYTHONPATH=${SPNCCI_PROJECT_ROOT_DIR}/spncci/script:${PYTHONPATH}

    Task parameters:

        # space parameters
        
        Nsigma_max (int): oscillator Nmax for lgi
        Nstep (int): step in N for many-body basis (1 or 2)
        N1v (int): valence shell oscillator N
        ...
          
  Language: Python 3

  A. E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT

  7/14/21 (aem): Created with code extracted from su3rme.py.

"""
import glob
import mcscript
import os

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
recoupler_executable_dir = os.path.join(os.environ["MCSCRIPT_INSTALL_HOME"],os.environ["CRAY_CPU_TARGET"], "su3shell/bin")
recoupler_executable = os.path.join(recoupler_executable_dir,"RecoupleSU3Operator")
# recoupler_executable = os.path.join(project_root,"lsu3shell","programs","upstreams","RecoupleSU3Operator")

# ... from spncci
generate_lsu3shell_relative_operators_executable = os.path.join(project_root,"spncci","programs","unit_tensors","generate_lsu3shell_relative_operators")

################################################################
# relative operator construction (input into LSU3Shell)
################################################################
def relative_operator_descriptor(task):
    """Generate descriptor string for use in relative operator archive
    filename.

    Returns:
        (string): descriptor
    """

    descriptor = "Nv{N1v:d}-Nsigmamax{Nsigma_max:02d}-Nstep{Nstep:d}".format(**task)
    return descriptor

def generate_relative_operators_for_recoupler(task):

    """
    Generates biquad representation of relative unit tensors, Arel, Brel and Nrel
    for input into recoupler_executable.

    Invokes generate_lsu3shell_relative_operators in spncci/programs/unit_tensors/.

    Output:
        relative_operators.dat : Control file.  Each line of file is:
                basename N lambda mu 2*S
        
            + basename is name of file containing operator with tensor character N(lambda,mu)S
            + list of basenames includes all relative unit tensors, Arel, Brel and Nrel.
        
        relative_unit_00000i.recoupler : biquad representation of each unit tensor indexed by i.
            Order of unit tensors determined by function u3shell::GenerateRelativeUnitTensorLabelsU3ST

        Arel.recoupler : biquad representation of relative symplectic raising operator
        Brel.recoupler : biquad representation of relative symplectic lowering operator
        Nrel.recoupler : biquad representation of relative number operator
    """

    last_recoupler_file_exists=os.path.exists("Nrel.recoupler")
    ## if recoupler files already exist, continue to next step.  Otherwise, generate files 
    if last_recoupler_file_exists:
        print("recoupler files found")
    else:
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
    """ Recouple biquads generated by Invoke lsu3shell recoupler code on relative unit
    tensors and symplectic raising/lowering/N operators.

    Invokes RecoupleSU3Operator.

    Arguments:
        relative_operator_basename_list (list) : list of operator names
    """
    
    # iterate over unit tensors
    for basename in relative_operator_basename_list:

        ## If restarting a run 
        ## Search for .PN and .PPNN files for given base name.  
        ## If found in current working directory,go to next base name

        file_exists=os.path.exists(f"{basename}.PPNN")
        file_exists&=os.path.exists(f"{basename}.PN")
        
        if not file_exists:            
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





def split_relative_operator_file(task):
    ## split relative_operators.dat into two control files, one for generators 
    ## and one for unit tensors.

    # open up file relative_operators.dat and separate out relative_generators and 
    # relative_operators 
    relative_filename="relative_operators.dat"
 
    unittensors=[];
    generators=[];
    file=open(relative_filename,'r')
    for line in file:
    # for line in lines:
        operator=line.split()
        if(operator[0]=="Arel" or operator[0]=="Brel" or operator[0]=="Nrel"):
            generators.append(line)
        else:
            unittensors.append(line)

    ## Write operator base names and labels to corresponding file 
    generator_filename="relative_generator_operators.dat"
    with open(generator_filename, 'w') as outstream:
        for line in generators:
            outstream.write(line)
    outstream.close()

    unittensor_filename="relative_unit_operators.dat"
    with open(unittensor_filename, 'w') as outstream:
        for line in unittensors:
            outstream.write(line)
    outstream.close()

    # moving relative_operators.dat since name will be re-used as a simlink to
    # relative_generators.dat or relative_unittensors.dat
    mcscript.call(["mv","relative_operators.dat","relative_operators-old.dat"])



def save_operator_files(task):
    """ 
    Create archive of relative operator files separated into two archives
        relative-unit-operators-{descriptor}.tgz :  Contains all relative unit tensor files
        relative-generator-operators-{descrptor}.tgz : Contains symplectic generators (Arel,Brel,Nrel)
    
    Results are saved to current run's results directory and need to be (manually) moved to
    directory specified in SPNCCI_LSU3SHELL_DIR for subsequent use.
    """

    ## Split up relative_operators.dat into two files, one for unit tensors, one for generators 
    split_relative_operator_file(task)

    ## Define descriptor for archive file
    descriptor = relative_operator_descriptor(task)
    
    # select unit tensor files to save
    unit_tensor_archive_list = glob.glob("relative_unit_operators.dat")
    unit_tensor_archive_list += glob.glob("relative_unit*.PN")
    unit_tensor_archive_list += glob.glob("relative_unit*.PPNN")
    
    unit_tensor_archive_filename = f"relative-unit-operators-{descriptor}.tgz"
    mcscript.call(
        ["tar", "-zcf", unit_tensor_archive_filename] + unit_tensor_archive_list
    )

    # Select generator files to save
    generators_archive_list = glob.glob("relative_generator_operators.dat")
    generators_archive_list += glob.glob("*rel.PN")
    generators_archive_list += glob.glob("*rel.PPNN")
    
    generators_archive_filename = f"relative-generator-operators-{descriptor}.tgz"
    mcscript.call(
        ["tar", "-zcf", generators_archive_filename] + generators_archive_list
    )

    # move archive to results directory (if in multi-task run)
    if (mcscript.task.results_dir is not None):
        mcscript.call(
            [
                "mv",
                "--verbose",
                unit_tensor_archive_filename,
                "--target-directory={}".format(mcscript.task.results_dir)
            ]
        )

    if (mcscript.task.results_dir is not None):
        mcscript.call(
            [
                "mv",
                "--verbose",
                generators_archive_filename,
                "--target-directory={}".format(mcscript.task.results_dir)
            ]
        )


def generate_relative_operators(task):
    """
    Generate the second quantized relative unit tensors and generators
    needed for calculating spncci seeds using lsu3shell/tools/SU3RME 
    
    Operators are save to 
        + archive containing all relative unit tensor files and list of operator names and quantum numbers 
            relative-unit-operators-{descriptor}.tgz 
        + archive containing symplectic generators (Arel,Brel,Nrel)
            relative-generator-operators-{descrptor}.tgz 
    
    Results are saved to current run's results directory and need to be (manually) moved to
    directory specified by SPNCCI_OPERATOR_DIR for subsequent use.
    """

    ## Generate expressions for operators in terms of creation and annihilation operators 
    
    generate_relative_operators_for_recoupler(task)

    ## Recouple products of creation and annihilation operators into order needed by SU3RME
    relative_operator_basename_list=read_relative_operator_basenames(task)
    recouple_operators(task,relative_operator_basename_list)
    
    ## Save operators in tar files 
    save_operator_files(task)



def retrieve_archived_files(filename):
    """
    Retrieve archived files and untar. 
    """
    
    archive_filename = mcscript.utils.search_in_subdirectories(
        operator_directory_list,
        operator_subdirectory_list,
        filename,
        error_message=f"relative unit operator archive file {filename} not found"
    )
    print(archive_filename)
    # extract archive contents
    mcscript.call(["tar", "xf", archive_filename])


def retrieve_unit_operator_files(task):
    """
    Retrieve relative unit tensor operator files.
    """
    descriptor = relative_operator_descriptor(task)
    filename=f"relative-unit-operators-{descriptor}.tgz"
    retrieve_archived_files(filename)


def retrieve_generator_operator_files(task):
    """
    Retrieve relative generator operator files.
    """
    descriptor = relative_operator_descriptor(task)
    filename=f"relative-generator-operators-{descriptor}.tgz"
    retrieve_archived_files(filename)



def retrieve_operator_files(task):
    """ 
    Retrieve relative unit tensor and generator operator files.
    
    """
    retrieve_unit_operator_files(task)
    retrieve_generator_operator_files(task)


if (__name__ == "__MAIN__"):
    pass
