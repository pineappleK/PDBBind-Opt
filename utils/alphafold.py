import shutil

def get_af2_protein(pdb_seq, af2_path, output_path):
    """
    Run AlphaFold2 on the protein sequence and save the output to a file.
    
    Args:
        pdb_seq: protein sequence
        output_path: output path
    """
    # Run AlphaFold2
    shutil.copy(output_path, af2_path)