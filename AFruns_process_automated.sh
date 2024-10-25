# Get list of all completed AF runs from output_data directory in text file
#find "/scratch/as3190/alphafold_multimer/output_data/" -name ranked_0.pdb > completed.txt

echo "1 Output file directory : $1"
echo "2 Completed file: $2"
echo "3 Prefix: $3"
echo "4 input_fasta_dir: $4"
echo "5 completed_fasta_dir: $5"
echo "6 source_dir: $6"
echo "7 new_dir: $7"
echo "8 fasta_file_name_completed: $8"
echo "9 completed_directory_list: $9"
echo "10 filter_AF_results_path: ${10}"
echo "10 filter_AF_results_path: ${10}"
echo "11 AF_output_process_path: ${11}"
echo "12 zip_msas_path: ${12}"


find "${1}" -name ranked_9.pdb > "${2}"

cp ${1}/${2} "/scratch/as3190/alphafold_multimer/"

echo "!!!copied ${1}/${2} to /scratch/as3190/alphafold_multimer/ "

awk -F/ '{print $2}' "${2}" | sed 's/$/.fasta/' > "${3}completed_fastas.txt"
awk -F/ '{print $2}' "${2}" > "${3}completed_directory.txt"

echo "clean completed dir list " "${3}completed_directory.txt"

if bash copy_fastas_runs_combined.sh "${4}" "${5}" "${6}" "${7}" "${8}" "${9}"; then
    echo 'Copying completed runs from /scratch/ to /projects/ AND Copying completed fasta files from /input_fasta/ to /completed_fasta/'

    echo "copying /scratch/as3190/alphafold_multimer/clean_completed_directory.txt to /projects/ccib/lamoureux/as3190/alphafold_multimer/output_data/";

    echo 'copying from /scratch/as3190/alphafold_multimer/' "${9}" to ${7}

    cp "/scratch/as3190/alphafold_multimer/${9}"  "${7}"  #"/projects/ccib/lamoureux/as3190/alphafold_multimer/output_data/";

    echo 'Running filter_AF_results_pkl now! Filtering the PKL files'
    if srun -N1 -n1 python filter_AF_results_pkl.py cluster "${10}" "${9}"; then
        chmod +x class_AF_output_processing.py
        srun -N1 -n1 python class_AF_output_processing.py cluster "${11}" "${9}"
    else
        echo 'Error: filter_AF_results_pkl.sh failed. Aborting further operations.'
    fi
else
    echo 'Error: copy_fastas_runs_combined.py failed. Aborting further operations.'
fi

#echo  'Done with all previous steps:: Zipping all MSA dirs'
#parent_dir2="${12}"
##parent_dir2="/projects/ccib/lamoureux/as3190/alphafold_multimer/output_data/"
#find "$parent_dir2" -mindepth 1 -type d -name 'msas' -exec tar czvf {}.tar.gz {} \; -exec rm -r {} \;
#



