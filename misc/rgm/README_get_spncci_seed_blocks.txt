First we need to calculate RMEs of symplectic generators and relative 2-body unit tensors in lsu3shell basis by running the script runsu3rme01. This script calls the lsu3shell programs SU3RME_MPI and ncsmSU3xSU2IrrepsTabular, which need to be copied to ~/Research/code/install/haswell/lsu3shell/bin, and spncci programs generate_lsu3shell_model_space, which needs to be copied to ~/Research/code/spncci/programs/unit_tensors/, and get_u3s_subspaces, which needs to be copied to ~/Research/code/spncci/programs/lgi/. Results, namely files lsu3shell_basis.dat, Brel.rme, Nrel.rme and relative_unit_<number>.rme, need to be moved to folder ~/Research/thesis/spncci/get_spncci_seed_blocks_dir/lsu3shell_rme/.

Then we need to calculate RMEs of one-body unit tensors in lsu3shell basis using SU3RME_Alexis and then SU3RME_Alexis_digest (see README.txt in ~/Research/thesis/lsu3shell/SU3RME_Alexis_dir/) and RMEs of two-body densities in lsu3shell basis using SU3RME_Alexis_2B and then SU3RME_Alexis_2B_digest (see README.txt in ~/Research/thesis/lsu3shell/SU3RME_Alexis_2B_dir/). Results, namely files a+N'taN_lm0_mu0_SS0_Tz.rme and N1_N2_N3_N4--lmf_muf_SSf--lmi_mui_SSi--irho0_lm0_mu0_SS0_Tz.rme need to be moved to folder ~/Research/thesis/spncci/get_spncci_seed_blocks_dir/lsu3shell_rme/.

Then batch script get_spncci_seed_blocks.sh should be edited accordingly. In this script, the last 3 numbers are Z (number of protons), N (number of neutrons) and Nmax. Then seeds can be calculated using get_spncci_seed_blocks by

sbatch get_spncci_seed_blocks.sh

Results appear in folder seeds, which needs to exist before.
