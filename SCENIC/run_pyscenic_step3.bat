#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=<set here> #180G
#SBATCH --time=10-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<mail id here>
#SBATCH --job-name="pyscenic"
#SBATCH --output="<log output.out>"  # lets call the file "<jobID>"
#SBATCH --error="<log output.err>"  # lets call the file "<jobID>"
#SBATCH --partition=prod
#if pyscenic installed as module
module load pySCENIC
module load python/3.8.1

res_dir="/pyscenic_results/results_dir"
input_dir="/pyscenic_results"
#download the tfdb files at https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/
tfdb1="hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather"
tfdb2="hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather"

##Step1
#NUMBA_THREADING_LAYER='omp' arboreto_with_multiprocessing.py -m grnboost2 -o ${res_dir}/integratedSubsampleV1_GRN_Step1.csv --num_workers 40 --seed 123123 ${input_dir}/integratedSubsampleV1.loom ${input_dir}/allTFs_hg38.txt
###Step2
#NUMBA_THREADING_LAYER='omp' pyscenic ctx --output ${res_dir}/integratedSubsampleV1_CTX_regulons_Step2.csv --num_workers 40 --mode custom_multiprocessing --annotations_fname ${input_dir}/motifs-v9-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname ${input_dir}/integratedSubsampleV1.loom ${res_dir}/integratedSubsampleV1_GRN_Step1.csv ${tfdb1} ${tfdb2}
#Step3
NUMBA_THREADING_LAYER='omp' pyscenic aucell --output ${res_dir}/integratedALL_aucell_mtx_step3_withAR.csv --num_workers 40 --seed 123123 ${input_dir}/integratedALL.loom ${res_dir}/integratedALL_CTX_regulons_Step2_withAR.csv --nes_threshold 1.5
##step3 save file in loom as output
NUMBA_THREADING_LAYER='omp' pyscenic aucell --output ${res_dir}/integratedALL_aucell_mtx_step3_withAR.loom --num_workers 40 --seed 123123 ${input_dir}/integratedALL.loom ${res_dir}/integratedALL_CTX_regulons_Step2_withAR.csv --nes_threshold 1.5

