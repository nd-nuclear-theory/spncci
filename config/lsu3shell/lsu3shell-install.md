Installing lsu3shell using ndconfig

  - 04/13/16 (mac): Written for use with University of Notre Dame CRC
    Linux cluster. First new build in a couple of years, with lots of
    new modules.
  - 09/5/16 (mac): Cut out old notes and notes on experimentation for 4/13 
    commit, and update config-ndcrc.mk.
  - 09/16/16 (mac): Update to lsu3shell 9/15 commit (de8ba2).  Remove instruction
    to update makefile.
  - 11/8/16 (mac): Paste in Anna's note from spncci/install_NDCRC.txt.
  - 11/17/16 (mac): Update to lsu3shell 11/13/16 commit (56981f) and
    current spncci config files (93ec11).
  - 12/22/16 (mac): Revise to serve as generic instructions across clusters,
  following shell/spncci project install.txt.
  - 01/24/17 (mac): Update path to config files (spncci/config/lsu3shell).
  - 05/22/17 (mac): Update config file description.
  - 05/28/17 (mac): Update config file symlinks.
  - 11/03/17 (mac): Update config file symlinks.
  - 06/19/20 (mac): Clean up after install at NERSC (LSU3shell 0a4db9f).
  
----------------------------------------------------------------

The following instructions assume you have already cloned spncci, so you have
access to the following directories:

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ~/code/spncci/config/ndconfig
  ~/code/spncci/config/lsu3shell
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

The compiler and library dependences for LSU3shell are similar to those for the
shell/spncci projects.  It is assumed you are already familiar with the
discussion and procedures in ndconfig/INSTALL.md.

1) Retrieving source

  Change to the directory where you want the repository to be installed,
  e.g.,

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  % cd ~/code
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  To clone a clean repository and check out the LSU3develop branch:
  
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  % git clone git://git.code.sf.net/p/lsu3shell/code lsu3shell
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  If you need the latest development version:

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  % git checkout -b LSU3develop origin/LSU3develop
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  Then change your working directory to lsu3shell for all the following steps.

2) Makefile configuration files

  We maintain the necessary config.mk file under spncci/config/lsu3shell.  (It
  could be moved to lsu3shell/devel/configs in the future.)

  First, remove the default config.mk that comes with lsu3shell.

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  % rm config.mk
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  For the makefile includes to work correctly, you first need to set up a
  symlink to the spncci config directory:

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  % ln -s ${HOME}/code/spncci/config
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  Then, make a symlink to the config file for your compiler environment, make a
  symlink to it:

  || @NDCRC: For compiling under gcc:
  || 
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ||   % ln -s config/lsu3shell/config-lsu3shell-gnu-ndcrc.mk config.mk
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  || @NERSC: For compiling under gcc:
  || 
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ||   % ln -s config/lsu3shell/config-lsu3shell-gnu-nersc.mk config.mk
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  If you need to know more, e.g., to set up a new configuration file for your
  own compiler environment, please see the documentation in the comments at the
  start of the main configuration file:

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  config/lsu3shell/lsu3shell-customizations.mk
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

3) Build

  It is first necessary to carry out the module load or setenv commands needed
  to access the various libraries.  These are likely to be the same as defined
  in ndconfig, with possibly slight modifications.

  || @NDCRC:
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  || % source ~/code/spncci/config/ndconfig/env-gnu-ndcrc.csh
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  || @NERSC:
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  || % source ~/code/spncci/ndconfig/env-gnu-nersc.csh
  || % module load cray-libsci
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  Although a parallel make is much faster, the fortran libraries cannot be made
  in parallel.  So, it is fastest if you build them first:

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  % make libeigensolver_MFDn libsu3lib
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  Then, you can do a parallel make.  If you only need to build the handful of
  supporting programs for spncci, you can use the special target
  install-for-spncci:

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  % make install-for-spncci -j 8
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  Or, for a full installation of LSU3shell:
  
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  % make install -j 8
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  || @NERSC: We need to keep binaries for different architectures
  || separate.  The files will be installed to
  || install/$(CRAY_CPU_TARGET)/bin, e.g., install/haswell/bin or
  || install/mic-knl/bin.


