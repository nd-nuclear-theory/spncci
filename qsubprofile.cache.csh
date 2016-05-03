#!/bin/csh

#$ -M amccoy@nd.edu # Email address for job notification
#$ -m abe			 # Send mail when job begins, ends and aborts
#$ -q short		  	 # Specify queue
#$ -N profile.ca.10      # Specify job name
#$ -pe smp 8

# Required modules
module load gcc/4.9.2	
module load gsl
module load boost/1.58
module load allinea
set data_file=$HOME/projects/spncci/libraries/spncci/lgi_test/lgi0.4.dat
foreach c (20 40)
	foreach N (10)
		foreach t (1 8)
			setenv OMP_SCHEDULE "dynamic,$c"
			setenv OMP_NUM_THREADS $t
			map --profile ./libraries/spncci/unit_tensor_test $data_file $N
		end
	end
end
