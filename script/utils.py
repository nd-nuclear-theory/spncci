""" utils -- scripting utility functions

  Language: Python 3

  M. A. Caprio
  Department of Physics, University of Notre Dame

  2/13/13 (mac): Extracted from job.py.
  5/28/13 (mac): Overhaul of subprocess invocation to use POpen.communicate().
  6/5/13 (mac): Absorbed into mcscript package.
  1/22/14 (mac): Python 3 update.
  5/14/15 (mac): 
      -- Insert "future" statements for Python 2 legacy support, conditional 
         use of decode by Python version.  
      -- Increased diagnostic output from subpath search utilities.
  6/13/16 (mac): Rename to utils.py as proper subpackage of mcscript.

TODO: replace QSUBM_MODULE_CMD with local hook; merge call_parallel and call_serial into call

"""

import glob
import math
import os
import shutil
import string
import subprocess
import sys
import time


################################################################
# module operations
################################################################

def module(args):
    """ module([arg1,arg2,...]) evaluates the code provided by the
    module command module_command python arg1 arg2 ...

    EX: module(["load", "mpich2"]) # to add mpiexec command to PATH
    """

    print (" ".join([ "module" ] + args))
    module_command = local.module_command()
    module_code_string = subprocess.check_output([ module_command, "python" ] + args)
    if (module_code_string != ""):
        print ("  Executing module code...")
        module_code = compile(module_code_string,"<string>","exec")
        eval (module_code)  # eval can crash on raw string, so compile first
    else:
        print ("  No module code to execute...")

################################################################
# subprocess execution
################################################################

# exception class for errors in script execution
class ScriptError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

def call(baseargs,input_lines=[],prefix=[],shell=False,cwd=None,print_stdout=True,check_return=True):
    """
    Invoke subprocess with simple batch input and output via
    POpen.communicate.  The subprocess arguments are obtained by
    joining the prefix list to the baseargs list.  Input lines and
    captured output are logged to stdout.

    baseargs: list of arguments for subprocess invocation
    
    input_lines: list of lines to be given to the subprocess as
    standard input (i.e., as list of strings, each of which is to be
    treated as one input line)

    print_stdout: boolean, whether or not to print subprocess's
    captured output to stdout

    check_resturn: whether or not to check subprocess's return value
    (and raise exception if nonzero)
    
    prefix: list of arguments to precede given arguments
    
    shell: pass-through to POpen, whether or not to launch subshell
    
    cwd: pass-through to POpen, current workind directory

    Returns: subprocess's captured output, as string

    Exceptions: If check_return is set and subprocess return is
    nonzero, raises a ScriptError exception.  Also raises a
    ScriptError exception if subprocess cannot be invoked in the first
    place.
    
    Limitation: All captured output is lost if process crashes.
    
    Examples:

    >>> mcscript.call(["cat"],input_lines=["a","b"]) # basic
    >>> mcscript.call(["catx"],input_lines=["a","b"]) # for execution failure
    >>> mcscript.call(["cat","badfile"],input_lines=["a","b"]) # for nonzero return
    
    """

    # set up arguments
    popenargs = prefix + baseargs 
    
    # set up input
    stdin_string = "".join([s + "\n" for s in input_lines])
    # encode as bytes
    #   else communicate complains under Python 3: 'str' does not
    #   support the buffer interface
    # caveat: bytes(stdin_string, "utf-8") fails Python 2 backward compatibility,
    # but encoding is required by Python 3
    stdin_bytes = bytes(stdin_string,"utf-8")

    # log header output
    print ("----------------------------------------------------------------")
    print ("Running %s" % str(popenargs))
    print ("Given standard input:")
    print (stdin_string)
    time_stamp()
    print ("----------------------------------------------------------------")
    sys.stdout.flush()
    
    # invoke
    try:
        # launch process
        process = subprocess.Popen(
            popenargs,
            stdin=subprocess.PIPE,     # to take input from communicate
            stdout=subprocess.PIPE,    # to send output to communicate
            stderr=subprocess.PIPE,    # separate stderr
            shell=shell,cwd=cwd,       # pass-through arguments 
            close_fds=True             # for extra neatness and protection
            )
    except OSError as err:
        print ("Execution failed:", err)
        raise ScriptError("execution failure")

    # communicate with process
    (stdout_bytes,stderr_bytes) = process.communicate(input=stdin_bytes)

    # process output
    # result of process.communicate was bytes (under Python 3)
    stdout_string = stdout_bytes.decode("utf-8")
    stderr_string = stderr_bytes.decode("utf-8")
    if (print_stdout):
        print (stdout_string)
    print (stderr_string)
    # handle return value
    returncode = process.returncode
    print ("Return code:", returncode)

    # finish logging
    print ("----------------------------------------------------------------")
    time_stamp()
    print ("----------------------------------------------------------------")
    sys.stdout.flush()  # just for good measure

    # return (or abort)
    if ( check_return and (returncode != 0) ):
        raise ScriptError("nonzero return")

    return stdout_string


def call_serial(baseargs,run,**kwargs):
    """ Invoke via mcscript.call, with serial prefix arguments.

    baseargs: list of arguments for subprocess invocation

    run: RunParameters object containing appropriate prefix

    kwargs: remaining keyword arguments for mcscrpt.call

    """

    if (self.batch_mode):
        if (self.parallel_epar != -1):
            # set configuration so serial codes will launch locally under this
            # daughter process (presumably already on compute node)
            serial_mode = "epar"
        else:
            # set configuration so serial codes will launch appropriately for
            # batch jobs at the present facility (perhaps on a separate
            # compute node if the job scripts run on dedicated nodes distinct
            # from the compute nodes)
            serial_mode = "batch"
    else:
        # local run on front end
        serial_mode = "local"

    serial_prefix = local.serial_prefix(serial_mode)
    return call(baseargs,prefix=serial_prefix,**kwargs)

def call_parallel(baseargs,run,**kwargs):
    """ Invoke via mcscript.call, with MPI parallel prefix arguments.

    baseargs: list of arguments for subprocess invocation

    run: RunParameters object containing appropriate prefix

    kwargs: remaining keyword arguments for mcscrpt.call

    """

    parallel_prefix = local.parallel_prefix(run)

    return call(baseargs,prefix=parallel_prefix,**kwargs)

################################################################
# input file generation
################################################################

def write_input(filename,input_lines=[],silent=False):
    """ Generate text file (typically an input file for a code to be
    invoked by the script), with contents given line-by-line as list
    of strings, and log to stdout.

    input_lines: list of strings for input lines

    silent: whether or not to suppress diagnostic output (e.g., for large sets of files
    """
    
    # set up input
    stdin_string = "".join([s + "\n" for s in input_lines])

    # produce diagnotic output
    if (not silent):
        print ("----------------------------------------------------------------")
        print ("Generating input file %s:" % filename)
        print (stdin_string)
        print ("----------------------------------------------------------------")

    # dump contents to file
    data_file = open(filename,"w")
    data_file.write(stdin_string)
    data_file.close()

################################################################
# timestamp utilities
################################################################

def time_stamp():
    """ Returns time stamp string.
    """
    
    return time.asctime()

def date_tag():
    """ Returns date tag string "YYMMDD".
    """
    return "%s" % time.strftime("%y%m%d")

################################################################
# string concatenation utilities
################################################################

def stringify(li,delimiter):
    """ Converts list entries to strings and joins with delimiter."""

    string_list = map(str,li)
    return delimiter.join(string_list)

def dashify(li):
    """ Converts list entries to strings and joins with "-"."""

    return stringify(li,"-")

def spacify(li):
    """ Converts list entries to strings and joins with spaces."""

    return stringify(li," ")

################################################################
# range construction utilities
################################################################

def value_range(x1,x2,dx,epsilon=0.00001):
    """ value_range(x1,x2,dx) produces float or integer results [x1, ..., x2] with step dx

    This is meant to emulate Mathematica's Range[x1,x2,x], with a "<=" upper bound, in
    contrast to Python's range, which is limited to integers and has an "<" upper bound.
    
    Limitations: Currently assumes dx is positive.  Upper cutoff is
    slightly fuzzy due to use of epsilon in floating point comparison.

    epsilon: tolerance for cutoff, as fraction of dx
    """

    value_list = []
    x=x1
    while (x < (x2 + dx * epsilon)):
        value_list.append(x)
        x += dx
    return value_list
    
def log_range(x1,x2,steps,base=2,power=1,first=True):
    """ log_range (x1,x2,steps,base=2,power=1,first=True) -> values returns a range
    starting from x1 and ending with x2 (inclusive) in logarithmic
    increments, at with steps increments per octave.  The first entry may be omitted by setting first=False.

    An octave refers to a doubling, if base=2 and power=1, or in general multiplication by base**power.

    print log_range(10,40,8)  -> 10..40 in 16 steps (as two octaves)
    print log_range(1,2,8,power=0.5) -> 1..2 in 16 steps (as two sqrt-octaves)
    """

    alpha = power**(-1) * math.log(x2/x1) / math.log(base)

    value_list = [
        x1 * base**(power*p)
        for p in value_range(0,alpha,1/steps)
        ]
    
    if (not first):
        value_list = value_list[1:]
        
    return value_list

################################################################
# arithmetic utilities
################################################################

def approx_equal(x,y,tol):
    """ approx_equal(x,y,tol) tests whether or not x==y to within tolerance tol
    """

    return (abs(y-x) <= tol)

def approx_less_equal(x,y,tol):
    """ approx_equal(x,y,tol) tests whether or not x<=y to within tolerance tol
    """

    return (x <= y+tol )

def approx_gtr_equal(x,y,tol):
    """ approx_equal(x,y,tol) tests whether or not x>=y to within tolerance tol
    """

    return (x >= y-tol )

def approx_less(x,y,tol):
    """ approx_equal(x,y,tol) tests whether or not x<y excluding values within tolerance tol
    """

    return not approx_gtr_equal(x,y,tol)

def approx_gtr(x,y,tol):
    """ approx_equal(x,y,tol) tests whether or not x>y excluding values within tolerance tol
    """

    return not approx_less_equal(x,y,tol)

################################################################
# default value utility
################################################################

def auto_value(x,default):
    """ auto_value(x,default) returns x, unless x is None, in which case it returns default
    """
    
    if x is None:
        return default
    else:
        return x

################################################################
# flow control utilities
################################################################

def ifelse(cond,a,b):
    if (cond):
        return a
    else:
        return b

################################################################
# path search utilities
################################################################

def subpath_search_filename(base_path,search_list,filename):
    """ subpath_search_filename(base_path,search_list,filename) -> name

    Searches for file base_path/subdir/filename, for subdir in
    search_list, and returns first match of existing file.  Raises
    exception on failure.
    """

    print ("Searching for file name...")
    print ("  Base path:", base_path)
    print ("  Search list:", search_list)
    print ("  Filename:", filename)

    # search in successive directories
    for subdir in search_list:
        qualified_name = os.path.join(base_path,subdir,filename)
        if (os.path.exists(qualified_name)):
            print ("  >>>>", qualified_name)
            return qualified_name

    # fallthrough
    print ("  No matching filename found...")
    raise ScriptError("no filename match on filename".format(filename))

def subpath_search_basename(base_path,search_list,basename):
    """ subpath_search_basename(base_path,search_list,basename) -> name

    Searches for file base_path/subdir/basename*, for subdir in
    search_list, and returns the base name base_path/subdir/basename
    for the first match of an existing file.  Raises
    exception on failure.
    """

    print ("Searching for file by base name...")
    print ("  Base path:", base_path)
    print ("  Search list:", search_list)
    print ("  Filename base:", basename)

    # search in successive directories file
    for subdir in search_list:
        qualified_name = os.path.join(base_path,subdir,basename)
        matches = glob.glob(qualified_name + "*")
        if (len(matches) is not 0):
            print ("  >>>>", qualified_name,len(matches))
            return qualified_name

    # fallthrough
    print ("  No matching filenames found...")
    raise ScriptError("no filename match on base filename {}".format(basename))

