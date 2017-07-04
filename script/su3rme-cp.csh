#!/bin/csh

# beware SPNCCI_SU3RME_DIR may be colon delimited list, not single
# directory, in which case an attempt to copy to
# ${SPNCCI_SU3RME_DIR_TARGET}/${run} wll lead to very creative target
# directory names

set run = $1
mkdir --parents ${SPNCCI_SU3RME_DIR}/${run}
nohup cp --update --verbose ${MCSCRIPT_WORK_HOME}/${run}/results/* ${SPNCCI_SU3RME_DIR_TARGET}/${run}
