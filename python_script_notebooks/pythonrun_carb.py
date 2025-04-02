from molecule_wrapper import init_molecule, try_update, accept_update, get_current_coords

# Initialize
coords = '/Users/josh/Documents/PhD/DevDungeon/carbonara/newFitData/SASDCP8/SAS_8_3segs/mol1_sub_0_end_xyz.dat'
fp =  '/Users/josh/Documents/PhD/DevDungeon/carbonara/newFitData/SASDCP8/fingerPrint1.dat'
init_molecule(fp, coords)

# Try update
coords = try_update(26)
if coords is not None:

    
    print("Update successful")
    print(coords)
