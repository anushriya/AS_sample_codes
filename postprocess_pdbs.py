import os
from pymol import cmd, util


def write_filter_pdb(pdb_path, complex_name, chainA, chainB, save_png=False):
    """
    Args:
        pdb_path: Path to pdb files
        complex_name: name of current complex
        chainA: List of residues for chain A to select for trimming
        chainB: List of residues for chain B to select for trimming
        save_png: if True,s ave png of plot
    Returns:

    """

    # starting_name = complex_name + "_ranked_0"
    starting_name = "ranked_0"

    for filename in os.listdir(pdb_path):
        if filename.startswith(starting_name) and filename.endswith(".pdb"):
            complete_path = os.path.join(pdb_path, filename)
            structure_name = 'my_structure'
            cmd.load(complete_path, structure_name)

            chainA_res = chainA
            chainB_res = chainB

            cmd.select('chainA_residue_selection',
                       f'{structure_name} and chain A and resi {"+".join(map(str, chainA_res))}')
            cmd.create('chainA_selected_structure', f'{structure_name} and chainA_residue_selection')

            cmd.select('chainB_residue_selection',
                       f'{structure_name} and chain B and resi {"+".join(map(str, chainB_res))}')
            cmd.create('chainB_selected_structure', f'{structure_name} and chainB_residue_selection')

            cmd.create('combined_structure', 'chainA_selected_structure or chainB_selected_structure')

            # Save trimmed pdbs
            output_pdb_file = "CMP_" + complex_name + "_filtered.pdb"
            output_complete_path = os.path.join(pdb_path, output_pdb_file)

            print("saving new pdb", output_complete_path)
            cmd.save(output_complete_path, "combined_structure")

            cmd.delete('my_structure')

            if save_png:
                cmd.show_as('cartoon', 'combined_structure')
                cmd.reset()
                cmd.center('combined_structure')
                cmd.set('depth_cue', 0)
                util.cbc(selection='combined_structure', first_color=7, quiet=1, _self=cmd)
                image_name = complex_name + "_filtered.png"
                image_complete_path = os.path.join(pdb_path, image_name)
                cmd.png(image_complete_path, width=230, height=230, dpi=-1, quiet=1, ray=0)

            cmd.delete("all")
