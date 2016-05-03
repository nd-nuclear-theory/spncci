#!/bin/csh

#$ -M amccoy@nd.edu # Email address for job notification
#$ -m abe			 # Send mail when job begins, ends and aborts
#$ -q long		  	 # Specify queue
#$ -N profile.ca.16.12     	 # Specify job name
#$ -pe smp 12

# Required modules
module load gcc/4.9.2	
module load gsl
module load boost/1.58
module load allinea
setenv OMP_NUM_THREADS 12
map --profile ./libraries/spncci/unit_tensor_test4 libraries/spncci/lgi_test/lgi0.4.dat 16
mv  *.map unit_tensor_profiles/.
