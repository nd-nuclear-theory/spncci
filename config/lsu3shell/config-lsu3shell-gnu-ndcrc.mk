include config/ndconfig/config-gnu-ndcrc.mk
include config/lsu3shell/lsu3shell-customizations.mk

# MPI configuration at ND CRC
CXX := mpicxx
FC := mpif77
fortran_libs += -lmpi_mpifh
