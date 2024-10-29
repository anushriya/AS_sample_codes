README

**Python script samples:**

`class_AF_output_processing.py` - A script for processing AlphaFold-Multimer outputs. Along with regular pLDDT and PAE parsing, it includes custom analysis for sub-domains using pLDDT. It then parses residues within sub-domains for low PAE score regions. 
Any residues not meeting the confidence threshold are trimmed in the final structure, which are processed with the `postprocess_pdbs.py` script.

`metrics_analysis.py` - Parses AlphaFold outputs from the above script to generate a pairwise relationship plot for all metrics considered (PAE, pLDDT, average PAE, inter- and intra-chain PAE averages, among others).

**Bash scripts Samples:**

`AFruns_process_automated.sh` and `run_AFruns_process_automated.sh` automate handling post AlphaFold-Multimer run files. 
This includes automated trimming of output files to retain only necessary metrics and transferring files to permanent directories.

**Others:**
`rcsb_api_search.py` - This script uses the RCSB API to query and process various ligands and proteins from the database directly.
