#!/bin/bash

# 1. Get the absolute path of THIS script's directory
# This ensures it works even if you call it from elsewhere
export FRAMEWORK_PATH=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)

# 2. Add the root and scripts folder to PYTHONPATH
# This lets you do "import my_module" from anywhere
# export PYTHONPATH="${FRAMEWORK_PATH}:${FRAMEWORK_PATH}/scripts:${PYTHONPATH}"
export PYTHONPATH="${FRAMEWORK_PATH}/lib:${FRAMEWORK_PATH}/scripts:${PYTHONPATH}"

# 3. Add 'scripts' to system PATH (only if not already there)
# This lets you run "myscript.py" directly as a command
if [[ ":$PATH:" != *":${FRAMEWORK_PATH}/scripts:"* ]]; then
    export PATH="${FRAMEWORK_PATH}/scripts:${PATH}"
fi

# 4. Optional: Custom prompt or alias
alias framework_cd="cd ${FRAMEWORK_PATH}"

echo "------------------------------------------------"
echo "Latinos framework loaded!"
echo "Location: ${FRAMEWORK_PATH}"
echo "------------------------------------------------"


source /cvmfs/sft.cern.ch/lcg/views/LCG_109a/x86_64-el9-gcc13-opt/setup.sh
if [[ "$(hostname -s)" == lxplus* ]]; then
    export MY_PATH_JOBS="/afs/cern.ch/user/a/abulla/private/LOG1/jobs/latinosRDF"
else
    export MY_PATH_JOBS="/gwpool/users/abulla/private/LOG1/jobs/latinosRDF" 
fi