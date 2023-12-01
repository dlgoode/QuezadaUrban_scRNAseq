#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=120G #180G
#SBATCH --time=10-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<mail id here>
#SBATCH --job-name="pyscenic"
#SBATCH --output="<log output.out>"  # lets call the file "<jobID>"
#SBATCH --error="<log output.err>"  # lets call the file "<jobID>"
#SBATCH --partition=prod

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

###default code
#NUMBA_THREADING_LAYER='omp' pyscenic ctx --output ${res_dir}/integratedALL_CTX_regulons_Step2.csv --num_workers 40 --mode custom_multiprocessing --annotations_fname ${input_dir}/motifs-v9-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname ${input_dir}/integratedALL.loom ${res_dir}/integratedAll_GRN_Step1.csv ${tfdb1} ${tfdb2}

###To check if AR is getting PASS
#--thresholds
NUMBA_THREADING_LAYER='omp' pyscenic ctx --output ${res_dir}/integratedALL_CTX_regulons_Step2_withthresholds3.csv --num_workers 40 --mode custom_multiprocessing --thresholds 0.50 0.75 0.90 --rank_threshold 10000 --nes_threshold 0.0 --annotations_fname ${input_dir}/motifs-v9-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname ${input_dir}/integratedALL.loom ${res_dir}/integratedAll_GRN_Step1.csv ${tfdb1} ${tfdb2}


