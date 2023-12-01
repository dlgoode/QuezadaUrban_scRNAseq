#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=<set here..> #180G
#SBATCH --time=10-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<mail id here>
#SBATCH --job-name="pyscenic"
#SBATCH --output="<log output.out>"  # lets call the file "<jobID>"
#SBATCH --error="<log output.err>"  # lets call the file "<jobID>"
#SBATCH --partition=prod
#if pySCENIC installed and configured as module...
module load pySCENIC 
module load python/3.8.1

res_dir="/pyscenic_results/results_dir"
input_dir="/pyscenic_results"

NUMBA_THREADING_LAYER='omp' arboreto_with_multiprocessing.py -m grnboost2 -o ${res_dir}/integratedAll_GRN_Step1.csv --num_workers 40 --seed 123123 ${input_dir}/integratedALL.loom ${input_dir}/allTFs_hg38.txt

#pyscenic grn --output ${res_dir}/integratedSubsampleV1_GRN_Step1.csv --num_workers 40 -m grnboost2 --seed 123123 ${input_dir}/integratedSubsampleV1.loom ${input_dir}/allTFs_hg38.txt


