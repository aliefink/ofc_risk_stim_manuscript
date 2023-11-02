#!/bin/sh

DATASET="${SGE_TASK}"

FUNCTION="pt_erp_caller"
SUBJECT="ST19"
METAID="ST19_B33"

# Analyze data set.
matlab_exec=matlab-2012b
func_call="${FUNCTION}('${SUBJECT}', '${METAID}','${DATASET}')"
echo ${func_call} > temp_${DATASET}.m
${matlab_exec} -nojvm -nodisplay -nosplash < temp_${DATASET}.m
rm temp_${DATASET}.m
