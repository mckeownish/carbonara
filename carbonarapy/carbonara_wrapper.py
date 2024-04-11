import subprocess
import glob
import os

def run_predict_structure(scatter_file, file_locs, initial_coords_file, paired_predictions, fixed_sections, no_structures,
                           within_monomer_hydro_cover, between_monomer_hydro_cover, kmin, kmax, max_no_fit_steps, prediction_file,
                             scatter_out, mixture_file, prev_fit_str, log_loc, end_line_prev_log, affine_trans, number_runs=3, verbose=False):
    
    existing_files = glob.glob(log_loc+'/*')

    # if len(existing_files) > 0:

    for file_path in existing_files:
        if os.path.isfile(file_path):
            os.remove(file_path)


    result_lst = []

    executable = "./predictStructureQvary"
    for i in range(1, number_runs+1):

        # Constructing the command with all args!
        cmd = [
            executable,
            scatter_file,
            file_locs,
            initial_coords_file,
            paired_predictions,
            fixed_sections,
            str(no_structures),
            within_monomer_hydro_cover,
            between_monomer_hydro_cover,
            str(kmin),
            str(kmax),
            str(max_no_fit_steps),
            f"{prediction_file}/mol{i}",
            f"{scatter_out}/scatter{i}.dat",
            mixture_file,
            prev_fit_str,
            f"{log_loc}/fitLog{i}.dat",
            end_line_prev_log,
            str(affine_trans)
        ]

        if verbose:
            print(f"Executing: {' '.join(cmd)}")

        # Execute the command
        result = subprocess.run(cmd, capture_output=True, text=True)
        result_lst.append({'cmd': cmd, 'stdout': result.stdout, 'stderr': result.stderr, 'returncode': result.returncode})

        # Check for errors
        if result.returncode != 0:
            # Log or print the error and exit
            print(f"Error executing command: {' '.join(cmd)}")
            print(f"Stderr: {result.stderr}")
            raise Exception("C++ execution failed")

        # check a logfile has been written
        output_file = f"{log_loc}/fitLog{i}.dat"
        if not os.path.exists(output_file):
            raise Exception(f"Expected output file {output_file} not found after execution.")

        if verbose:
            print(f"Execution completed for run number {i}")
    
    return result_lst

# Example usage
# run_predict_structure("/path/to/Saxs.dat", "/path/to/fileLocs", "frompdb", "False", "/path/to/fixedsections.dat", 1, "none", "none", 0.01, 0.25, 5, "/path/to/predictionFile", "/path/to/scatterOut", "/path/to/mixtureFile.dat", "/path/to/prevFitStr", "/path/to/logLoc", "null", "False")
