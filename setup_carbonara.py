#!/usr/bin/env python3

import argparse
import os
import shutil
import sys
import wiggle.CarbonaraDataTools as cdt

def main():
    parser = argparse.ArgumentParser(description='Setup Carbonara processing pipeline')
    parser.add_argument('-p', '--pdb', required=True, help='Path to input PDB file')
    parser.add_argument('-s', '--saxs', required=True, help='Path to input SAXS data file')
    parser.add_argument('-n', '--name', required=True, help='Name for this protein/refinement')
    parser.add_argument('-d', '--dir', default=os.getcwd(), help='Base directory (default: current directory)')
    args = parser.parse_args()
    
    try:
        # Setup master directory
        fit_master_dir = cdt.setup_fit_master_dir(root_dir=args.dir, fit_master_name='carbonara_runs')
        
        # Setup refinement directory
        refine_dir = cdt.setup_refinement_dir(args.name, fit_master_dir)
        print(f"Created directory structure in: {refine_dir}")
        
        # Process PDB and extract structure information - the secondary structure identification is not amazing atm
        coords_chains, sequence_chains, secondary_structure_chains, missing_residues_chains = cdt.pull_structure_from_pdb(args.pdb)
        
        # Give warning if missing residues are found
        for coords in coords_chains:
            breaking_indices = cdt.missing_ca_check(coords, threshold_dist_Å=10)
            if len(breaking_indices) > 0: print("Warning: Missing segments of chain found: ", len(breaking_indices), breaking_indices)
        
        # write coordinates files ( coordinates1.dat, coordinates2.dat, etc [new files for each chain])
        coords_files = []
        for i, coords in enumerate( coords_chains ):
            coords_files.append( cdt.write_coordinates_file(coords, working_path=refine_dir, carb_index=i+1) )
        
        # Write fingerprint file
        number_of_chains = len(coords_chains)
        fingerprint_file = cdt.write_fingerprint_file( number_chains=number_of_chains, sequence=sequence_chains,
                                                       secondary_structure=secondary_structure_chains, working_path=refine_dir )

        # write mixture file - used for ensemble refinement, currently not used - writes 1 to mixture file
        mixture_file = cdt.write_mixture_file(working_path=refine_dir)
        
        # Copy SAXS file to Saxs.dat (this is the file that Carbonara will use)
        shutil.copy2(args.saxs, os.path.join(refine_dir, 'Saxs.dat'))

        # auto select flexible linker chains that dont break inter-beta sheets
        varying_linker_chains = []
        for coord_file in coords_files: 
            varying_linker_chains.append( cdt.auto_select_varying_linker(coord_file, fingerprint_file) )

        # write flexible linkers to files (varysections1.dat, varysections2.dat, etc [each file is for a different chain])
        varying_section_files = []
        for varying_linkers in varying_linker_chains:
            varying_section_files.append( cdt.write_varysections_file(varying_linkers, refine_dir) )


        print("\nSetup completed successfully!")
        print(f"All files have been created in: {refine_dir}")
        
    except Exception as e:
        print(f"Error during setup: {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main() 