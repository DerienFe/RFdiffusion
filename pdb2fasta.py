import os
import argparse
from Bio import SeqIO
from Bio.PDB import PDBParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def extract_fasta_from_pdb(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('', pdb_file)
    sequences = {}
    
    for model in structure:
        for chain in model:
            seq = Seq("".join([residue.resname for residue in chain if residue.id[0] == " "]))
            chain_id = chain.id
            sequences[chain_id] = seq

    return sequences

def pdb_to_fasta(input_dir, output_dir):
    for filename in os.listdir(input_dir):
        if filename.endswith('.pdb'):
            pdb_file = os.path.join(input_dir, filename)
            sequences = extract_fasta_from_pdb(pdb_file)
            
            # Prepare the output file path with the same name but .fasta extension
            fasta_filename = f"{os.path.splitext(filename)[0]}.fasta"
            output_file = os.path.join(output_dir, fasta_filename)
            
            fasta_sequences = [
                SeqRecord(sequence, id=f"{filename}_{chain_id}", description="")
                for chain_id, sequence in sequences.items()
            ]
            
            # Write each PDB file's sequence to its own .fasta file
            with open(output_file, "w") as fasta_output:
                SeqIO.write(fasta_sequences, fasta_output, "fasta")
            
            print(f"FASTA sequence written to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert each PDB file in a directory to a separate FASTA file")
    parser.add_argument("-input_dir", required=True, help="Directory containing PDB files")
    parser.add_argument("-output_dir", required=True, help="Directory to save the output FASTA files")
    
    args = parser.parse_args()
    
    # Ensure the output directory exists
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Run the function with the provided arguments
    pdb_to_fasta(args.input_dir, args.output_dir)
