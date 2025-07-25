import os
import subprocess
import argparse
import tempfile
from Bio import PDB

def run_anarci(chain_fasta_path, output_prefix, scheme):
    """Runs ANARCI and returns the path to the main output file."""
    anarci_output_path = f"{output_prefix}" 
    anarci_cmd = [
        "ANARCI", "-i", chain_fasta_path, "-o", output_prefix, "--scheme", scheme
    ]
    try:
        subprocess.run(anarci_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if os.path.exists(anarci_output_path):
            return anarci_output_path
    except subprocess.CalledProcessError as e:
        # ANARCI often prints warnings to stderr for non-IG C-terminal tails, but still produces a valid file.
        if "Warning: Non IG chains cannot be numbered" in e.stderr and os.path.exists(anarci_output_path):
            return anarci_output_path
        else:
            print(f"ANARCI Error: {e.stderr}")
    return None

def parse_anarci_output(anarci_file):
    """Parses ANARCI output to get a mapping of residue index to ANARCI numbering, sequence, and chain type (H/L)."""
    anarci_map = {}
    anarci_seq = []
    chain_type = None
    if not anarci_file:
        return anarci_map, "".join(anarci_seq), chain_type
        
    with open(anarci_file, 'r') as f:
        for line in f:
            if line.startswith('#') or line.strip() == "":
                continue
            
            parts = line.split()
            if len(parts) >= 3 and parts[-1] != '-': # Ignore deletions
                # Set chain type from the first character of the first valid line
                if chain_type is None:
                    chain_type = parts[0]

                res_num_str = parts[1]
                icode = ' '
                if len(parts) == 4:
                    icode = parts[2]
                
                anarci_map[len(anarci_seq)] = (int(res_num_str), icode)
                anarci_seq.append(parts[-1])
                
    return anarci_map, "".join(anarci_seq), chain_type

def renumber_and_merge_chains(pdb_path, scheme):
    """Identifies, renames (H/L), renumbers, and merges antibody chains from a PDB file."""
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('input_pdb', pdb_path)
    
    pdb_stem = os.path.splitext(os.path.basename(pdb_path))[0]
    output_filename = f"{pdb_stem}_{scheme}_renumbered.pdb"
    
    renumbered_chains = []

    with tempfile.TemporaryDirectory() as temp_dir:
        for model in structure:
            for chain in model:
                original_chain_id = chain.get_id()
                print(f"Processing Chain {original_chain_id}...")
                
                pdb_seq = "".join(PDB.Polypeptide.three_to_one(res.get_resname()) for res in chain if PDB.is_aa(res))
                if not pdb_seq:
                    continue

                temp_fasta_path = os.path.join(temp_dir, f"chain_{original_chain_id}.fasta")
                with open(temp_fasta_path, 'w') as f_fasta:
                    f_fasta.write(f">{original_chain_id}\n{pdb_seq}\n")

                anarci_output_prefix = os.path.join(temp_dir, f"anarci_{original_chain_id}")
                anarci_file = run_anarci(temp_fasta_path, anarci_output_prefix, scheme)

                if anarci_file:
                    anarci_map, anarci_seq, chain_type = parse_anarci_output(anarci_file)
                    print(f"  Chain {original_chain_id} identified as a {chain_type}-type antibody chain.")
                    
                    if pdb_seq.startswith(anarci_seq):
                        print(f"  Sequence is a partial match. Renumbering and renaming to Chain {chain_type}...")
                        
                        # Build a completely new chain to avoid residue ID conflicts
                        new_renumbered_chain = PDB.Chain.Chain(chain_type)
                        aa_residues = [res for res in chain if PDB.is_aa(res)]
                        
                        num_anarci_residues = len(anarci_seq)
                        
                        # 1. Renumber the part matched by ANARCI
                        for i in range(num_anarci_residues):
                            new_residue = aa_residues[i].copy()
                            if i in anarci_map:
                                new_res_num, new_icode = anarci_map[i]
                                new_residue.id = (new_residue.id[0], new_res_num, new_icode)
                            new_renumbered_chain.add(new_residue)
                        
                        # 2. Sequentially renumber the C-terminal tail
                        if num_anarci_residues > 0:
                            last_anarci_res_num, _ = anarci_map[num_anarci_residues - 1]
                            next_res_num = last_anarci_res_num + 1

                            for i in range(num_anarci_residues, len(aa_residues)):
                                new_residue = aa_residues[i].copy()
                                new_residue.id = (new_residue.id[0], next_res_num, ' ')
                                new_renumbered_chain.add(new_residue)
                                next_res_num += 1
                        
                        renumbered_chains.append(new_renumbered_chain)
                    else:
                        print(f"  ERROR: Residue sequence is a complete mismatch for Chain {original_chain_id}. Skipping.")
                else:
                    print(f"  Chain {original_chain_id} is not an antibody chain or ANARCI failed. Skipping.")
    
    if renumbered_chains:
        print(f"\nMerging {len(renumbered_chains)} renumbered chains into {output_filename}...")
        final_structure = PDB.Structure.Structure("renumbered")
        final_model = PDB.Model.Model(0)
        for chain_obj in renumbered_chains:
            final_model.add(chain_obj)
        final_structure.add(final_model)

        io_final = PDB.PDBIO()
        io_final.set_structure(final_structure)
        io_final.save(output_filename)
        print("Done!")
    else:
        print("\nNo antibody chains were successfully renumbered. No output file created.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Renumber and rename antibody chains (H/L) in a PDB file using ANARCI and merge them.")
    parser.add_argument("pdb_file", help="Path to the input PDB file.")
    parser.add_argument("--scheme", "-s", default="kabat", help="Numbering scheme to use (e.g., kabat, chothia). Default is 'kabat'.")
    args = parser.parse_args()
    
    renumber_and_merge_chains(args.pdb_file, args.scheme)
