#!/bin/bash
#SBATCH --partition=main

##SBATCH --partition=cmain

#SBATCH --job-name=directory2cleanAD_projects
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=3-00:00:00

##SBATCH --mem=30G                # Real memory (RAM) required (MB)
#SBATCH --mem=60G                # Real memory (RAM) required (MB)
#SBATCH --output=slurm_log/slurm.%N.%j.%x.out
#SBATCH --error=slurm_log/slurm.%N.%j.%x.err
#SBATCH --export=ALL
#SBATCH --constraint=oarc

pwd
source /home/as3190/anaconda3/etc/profile.d/conda.sh
conda activate dock3D_cud10

echo "FOR POSITIVE INTERACTIONS ONLY"

cur_output_file="output_data/"
curr_completed_file='completed.txt'
curr_prefix="clean_"

input_fasta_dir="input_fasta"
completed_fasta_dir="completed_fasta"

source_dir="/scratch/as3190/alphafold_multimer/output_data/"
new_dir="/projects/ccib/lamoureux/as3190/alphafold_multimer/output_data"

fasta_file_name_completed=${curr_prefix}'completed_fastas.txt'
completed_list=${curr_prefix}'completed_directory.txt'

filter_AF_results_path='/projects/ccib/lamoureux/as3190/alphafold_multimer/output_data/'
AF_output_process_path='/projects/ccib/lamoureux/as3190/alphafold_multimer/output_data/'

zip_msas_path="/projects/ccib/lamoureux/as3190/alphafold_multimer/output_data/"

run_afruns_process() {
  local cur_output_file=${1}
  local curr_completed_file=${2}
  local curr_prefix=${3}
  local input_fasta_dir=${4}
  local completed_fasta_dir=${5}
  local source_dir=${6}
  local new_dir=${7}
  local fasta_file_name_completed=${8}
  local completed_list=${9}
  local filter_AF_results_path=${10}
  local AF_output_process_path=${11}
  local zip_msas_path=${12}

  bash AFruns_process_automated.sh \
                                  "${cur_output_file}" \
                                  "${curr_completed_file}"\
                                  "${curr_prefix}"\
                                  "${input_fasta_dir}" \
                                  "${completed_fasta_dir}"\
                                  "${source_dir}" \
                                  "${new_dir}" \
                                  "${fasta_file_name_completed}"\
                                  "${completed_list}" \
                                  "${filter_AF_results_path}" \
                                  "${AF_output_process_path}" \
                                  "${zip_msas_path}"\ ;
}

run_afruns_process \
                  "${cur_output_file}" \
                  "${curr_completed_file}"\
                  "${curr_prefix}" \
                  "${input_fasta_dir}" \
                  "${completed_fasta_dir}" \
                  "${source_dir}" \
                  "${new_dir}" \
                  "${fasta_file_name_completed}"\
                  "${completed_list}" \
                  "${filter_AF_results_path}" \
                  "${AF_output_process_path}" \
                  "${zip_msas_path}"
