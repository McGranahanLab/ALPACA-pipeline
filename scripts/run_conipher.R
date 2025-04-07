library(CONIPHER)
library(tidyverse)

args = commandArgs(TRUE)

case_id = args[1]
prefix = args[2]
input_tsv_loc = args[3]
out_dir = args[4]
print('Running CONIPHER with the following parameters:')
print(paste('case_id:', case_id))
print(paste('prefix:', prefix))
print(paste('input_tsv_loc:', input_tsv_loc))
print(paste('out_dir:', out_dir))

CONIPHER::conipher_run(case_id = case_id,prefix = prefix,out_dir = out_dir,input_tsv_loc = input_tsv_loc)

print('Done running CONIPHER')