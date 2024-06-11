import CarbonaraDataTools as CDT
import numpy as np
import os, subprocess, shutil, sys
from datetime import datetime

### sys.argv[1] = PDB File
### sys.argv[2] = SAXS File

def defaultRun(pdb_file,saxs_file):
    ### File paths
    mol_name = os.path.basename(pdb_file)[:-4]
    now = datetime.now()
    run_name = now.strftime("%Y%m%d_%H%M%S_%f")
    working_path = 'newFitData/'+mol_name
    if not os.path.isdir(working_path):
        os.makedirs(working_path)
    ### Organise chain coords
    coords_chains_in, sequence_chains_in, secondary_structure_chains_in, missing_residues_chains_in = CDT.pull_structure_from_pdb(pdb_file)
    for j in range(len(coords_chains_in)):
        breaking_indices = CDT.missing_ca_check(coords_chains_in[j])
        if len(breaking_indices) > 0:
            chains = [i+1 for i in range(len(breaking_indices)+1)]
            print("We think your PDB file has", len(breaking_indices)+1, 'chains')
            coords_chains_in[j], sequence_chains_in[j] = CDT.break_into_chains(coords_chains_in[j],sequence_chains_in[j],breaking_indices)

    coords_chains,sequence_chains = [],[]
    for i in range(len(coords_chains_in)):
        if isinstance(coords_chains_in[i],list):
            for j in range(len(coords_chains_in[i])):
                coords_chains.append(coords_chains_in[i][j])
                sequence_chains.append(sequence_chains_in[i][j])
        else:
            coords_chains.append(coords_chains_in[i])
            sequence_chains.append(sequence_chains_in[i])

    coords=[]
    for i in range(len(coords_chains)):
        for j in range(len(coords_chains[i])):
            coords.append(coords_chains[i][j])
    coords = np.array(coords)
    ### Use kappa tau method to identify secondary structure
    secondarysplit= [CDT.findSimilar(CDT.getKapTauList(coords_chains[i]),0.2) for i in range(len(coords_chains))]
    ### Writes all the coords, fingerprint, mixture & SAXS in Carbonara's expected input / file struct
    CDT.write_fingerprint_file(len(coords_chains), sequence_chains, secondarysplit, working_path)    
    CDT.write_coordinates_file(coords,working_path)
    CDT.write_mixture_file(working_path)
    CDT.write_saxs(SAXS_file, working_path)
    chains = [i+1 for i in range(len(coords_chains))]
    varyingSection_tensor = [[] for i in range(len(breaking_indices)+1)]
    ### SAXS range
    minq = 0.01
    startq = 0.15
    maxq = round(np.genfromtxt(working_path+'/Saxs.dat')[:,0].max(),2)
    ### Varying Sections
    allowed_linker, linker_indices = CDT.find_non_varying_linkers(working_path+'/coordinates1.dat',working_path+'/fingerPrint1.dat')
    CDT.write_varysections_file(allowed_linker,working_path)
    ### Run length and repeats
    no_runs = '10'
    fit_steps = '10000'
    ### Rotations for multimers, if monomer will do nothing
    rotation = True
    ### Writhe the executable
    CDT.write_run_sh_file(working_path,
                                    mol_name,
                                    run_name,
                                    1,
                                    minq,
                                    maxq,
                                    startq,
                                    int(fit_steps),
                                    int(no_runs),
                                    pairedQ=False,
                                    rotation=rotation)
    
