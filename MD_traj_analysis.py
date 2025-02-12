"""
This Python script is designed for processing and analyzing molecular dynamics (MD) trajectories of protein-nucleic acid complexes.
It automates the workflow of trajectory preparation, essential analyses, and visualization, focusing on extracting insights
into the dynamics and energetics of these biomolecular systems.

The script performs the following main tasks:

1. Trajectory Processing (using the `process_trajectory` function):
   - Converts PDB structure to GRO format if a GRO file is not initially available.
   - Reads molecular dynamics trajectory files in GRO and XTC formats.
   - Selects atoms belonging to protein and nucleic acid components from the system.
   - Applies a series of transformations to the trajectory to prepare it for analysis:
     - Unwraps the complex to correct for periodic boundary conditions.
     - Centers the nucleic acid within the simulation box based on its mass center.
     - Wraps the complex back into the box, ensuring molecular integrity.
     - Fits and rotates the protein based on a reference protein structure to remove overall translational and rotational motion,
       facilitating the analysis of internal dynamics.
   - Saves the processed trajectory as new GRO and XTC files, ready for subsequent analysis.

2. Trajectory Analysis (using the `analyze_trajectory` function):
   - Calculates key metrics to characterize the dynamics of the nucleic acid component:
     - Root Mean Square Deviation (RMSD) of the nucleic acid over time, indicating conformational changes.
     - Root Mean Square Fluctuation (RMSF) per residue of the nucleic acid, highlighting flexible regions.
     - Energy landscape projections based on RMSD and Radius of Gyration (Rg) of the nucleic acid. These landscapes provide insights
       into the free energy distribution and conformational preferences of the nucleic acid during the simulation.
   - Generates and saves plots for RMSD over time, RMSF per residue, and the energy landscapes.
   - Outputs the energy landscape data to text files for further quantitative analysis.
   - Consolidates all generated plots into a single image file for easy review and saves a copy to a designated image output directory.

3. Batch Processing via Main Function (`main` function):
   - Sets user-defined parameters such as simulation temperature, output filename prefixes, and image output directory.
   - Allows for flexible specification of molecular dynamics simulation directories to be processed, either through automatic directory
     discovery based on naming conventions or by manually listing directory paths.
   - Iterates through each specified simulation directory to:
     - Determine if trajectory processing and/or analysis is needed based on the presence of processed files or simulation state files.
     - Executes the `process_trajectory` and `analyze_trajectory` functions accordingly.
   - Skips directories that do not contain relevant simulation files or are not valid directories, providing informative messages.
   - Completes the analysis for all specified directories in a batch manner.

This script is intended to streamline the analysis of molecular dynamics simulations, providing researchers with automated tools
to assess the conformational dynamics and energetic landscapes of protein-nucleic acid complexes, which are crucial for understanding
their biological functions and interactions.
"""
import os
import MDAnalysis as mda
from MDAnalysis.transformations import unwrap, center_in_box, wrap, fit_rot_trans
from MDAnalysis.analysis.rms import RMSD, RMSF
from MDAnalysis.analysis import align
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants  # Import physical constants library

# --- User-modifiable parameters ---
TEMPERATURE_K = 300  # Simulation temperature, Kelvin (K)
OUTPUT_FILENAME_PREFIX = 'nucleic_analysis'  # Output filename prefix
IMAGE_OUTPUT_DIR_NAME = 'outpng' # Image output directory name
# --- End of parameter modifications ---


def process_trajectory(md_dir):
    """
    Process molecular dynamics trajectory, extract protein and nucleic acid, and perform fit and rotate.

    Parameters:
    md_dir (str): Path to the molecular dynamics simulation directory.
    """
    print(f"Start processing trajectory files, directory: {md_dir}")

    input_PDB = os.path.join(md_dir, "input.pdb")
    input_gro = os.path.join(md_dir, "input.gro")
    input_xtc = os.path.join(md_dir, "trajectory.xtc")
    pro_nucle_gro = os.path.join(md_dir, "protein_nucleic.gro")
    pro_nucle_xtc = os.path.join(md_dir, "protein_nucleic.xtc")
    pro_nucle_gro_fit = os.path.join(md_dir, "protein_nucleic_fit.gro")
    pro_nucle_xtc_fit = os.path.join(md_dir, "protein_nucleic_fit.xtc")

    # 1. PDB to GRO conversion (if gro file does not exist)
    if not os.path.exists(input_gro):
        print(f"Converting PDB file to GRO file: {input_PDB} -> {input_gro}")
        try:
            u_pdb = mda.Universe(input_PDB, format="pdb")
            u_pdb.atoms.write(input_gro)
        except Exception as e:
            print(f"PDB to GRO conversion failed: {e}")
            return

    # 2. Read gro and xtc files
    try:
        u = mda.Universe(input_gro, input_xtc)
    except Exception as e:
        print(f"Failed to read GRO and XTC files: {e}")
        return

    # 3. Select protein and nucleic acid and write new gro and xtc (no frame skipping)
    protein_and_nucleic = u.select_atoms('protein or nucleic')
    protein_and_nucleic.dimensions = u.dimensions # Keep box dimension information
    try:
        protein_and_nucleic.write(pro_nucle_gro) # Save as new gro
        with mda.Writer(pro_nucle_xtc, protein_and_nucleic.atoms.n_atoms) as W: # Save as new xtc (skipping frames every 100 frames)
            for ts in u.trajectory[::100]: # Here is still frame skipping, according to your notebook code, keep frame skipping
                W.write(protein_and_nucleic.atoms)
        print(f"Saved GRO and XTC files for protein and nucleic acid: {pro_nucle_gro}, {pro_nucle_xtc}")
    except Exception as e:
        print(f"Failed to save GRO/XTC files for protein and nucleic acid: {e}")
        return


    # 4. Fit and Rotate trajectory processing
    try:
        u_fit = mda.Universe(pro_nucle_gro, pro_nucle_xtc) # Use new gro and xtc
        protein = u_fit.select_atoms('protein')
        dna = u_fit.select_atoms('nucleic')
        complex_fit = u_fit.select_atoms('protein or nucleic')
        u_fit.atoms.guess_bonds()
        ref_u = u_fit.copy()
        reference = ref_u.select_atoms("protein")

        workflow = (
            unwrap(complex_fit),
            center_in_box(dna, center='mass'),
            wrap(complex_fit, compound='fragments'),
            fit_rot_trans(protein, reference)
        )
        u_fit.trajectory.add_transformations(*workflow)

        u_fit.atoms.write(pro_nucle_gro_fit) # Save fit gro
        with mda.Writer(pro_nucle_xtc_fit, u_fit.atoms.n_atoms) as W: # Save fit xtc (all frames)
            for ts in u_fit.trajectory:
                W.write(u_fit.atoms)
        print(f"Saved fit & rotated GRO and XTC files: {pro_nucle_gro_fit}, {pro_nucle_xtc_fit}")

    except Exception as e:
        print(f"Trajectory Fit & Rotate processing failed: {e}")
        return

    print(f"Trajectory file processing completed, directory: {md_dir}")



def analyze_trajectory(pro_nucle_gro_fit, pro_nucle_xtc_fit, output_dir, output_filename_prefix, image_output_dir):
    """
    Analyze molecular dynamics trajectory, calculate RMSD, RMSF and Energy Landscape.

    Parameters:
    pro_nucle_gro_fit (str): Path to the fitted GRO file.
    pro_nucle_xtc_fit (str): Path to the fitted XTC file.
    output_dir (str): Directory to save output files.
    output_filename_prefix (str): Output filename prefix.
    image_output_dir (str): Image output directory.
    """
    print(f"Start analyzing trajectory files: GRO={pro_nucle_gro_fit}, XTC={pro_nucle_xtc_fit}")

    # Load trajectory and reference structure
    try:
        u = mda.Universe(pro_nucle_gro_fit, pro_nucle_xtc_fit)
        ref = mda.Universe(pro_nucle_gro_fit)  # Use initial structure as reference
    except Exception as e:
        print(f"Failed to load trajectory files: {e}")
        return

    nucleic = u.select_atoms('nucleic')
    ref_nucleic = ref.select_atoms('nucleic')

    # Create a figure and axes object with 1 row and 4 columns of subplots
    fig, axes = plt.subplots(1, 4, figsize=(24, 6))  # Adjust figsize to make image wider

    # 1. Calculate RMSD over time and plot
    print("Start calculating RMSD...")
    mobile_atoms_rmsd = u.select_atoms('nucleic')
    reference_atoms_rmsd = ref.select_atoms('nucleic')

    rmsd_analysis = RMSD(mobile_atoms_rmsd, reference_atoms_rmsd,
                          center=True,
                          superposition=True)
    try:
        rmsd_analysis.run()
        rmsd_values = rmsd_analysis.rmsd[:, 2]
        time_steps = range(len(rmsd_values))

        axes[0].plot(time_steps, rmsd_values)  # Use axes[0] to plot
        axes[0].set_xlabel('Time Step')
        axes[0].set_ylabel('RMSD (Å)')
        axes[0].set_title('RMSD of Nucleic Over Time')
        print(f"RMSD over time plot is ready")
    except Exception as e:
        print(f"RMSD calculation failed: {e}")


    # 2. Calculate RMSF and plot
    print("Start calculating RMSF...")
    try:
        aligner_rmsf = align.AlignTraj(u, ref, select='nucleic', in_memory=True).run()  # Align for RMSF calculation
        rmsf_analyzer = RMSF(nucleic).run()

        residues = nucleic.residues.resids
        res_rmsf = []
        for res in nucleic.residues:
            atom_indices = res.atoms.ix - nucleic[0].ix
            res_rmsf.append(rmsf_analyzer.results.rmsf[atom_indices].mean())

        axes[1].plot(residues, res_rmsf, 'o-')  # Use axes[1] to plot
        axes[1].set_xlabel('Residue ID')
        axes[1].set_ylabel('RMSF (Å)')
        axes[1].set_title('RMSF per Residue (Nucleic)')
        print(f"RMSF per Residue plot is ready")
    except Exception as e:
        print(f"RMSF calculation failed: {e}")


    # 3. Calculate RMSD-based Energy Landscape and plot
    print("Start calculating RMSD-based Energy Landscape...")
    try:
        rmsd_analyzer_el = RMSD(u, ref, select='nucleic', in_memory=True)  # Recalculate RMSD, or use the previously calculated rmsd_analyzer
        rmsd_analyzer_el.run()
        rmsd_values_el = rmsd_analyzer_el.rmsd[:, 2]

        kB = scipy.constants.Boltzmann
        kT = kB * TEMPERATURE_K

        hist_rmsd_el, bin_edges_rmsd_el = np.histogram(rmsd_values_el, bins=50, density=True)
        bin_centers_rmsd_el = (bin_edges_rmsd_el[:-1] + bin_edges_rmsd_el[1:]) / 2

        free_energy_rmsd_el = -kT * np.log(hist_rmsd_el)
        free_energy_kJ_mol_rmsd_el = free_energy_rmsd_el / 1000 / 4.184
        min_free_energy_rmsd_el = np.min(free_energy_kJ_mol_rmsd_el[np.isfinite(free_energy_kJ_mol_rmsd_el)])
        free_energy_land_rmsd_el = free_energy_kJ_mol_rmsd_el - min_free_energy_rmsd_el

        axes[2].plot(bin_centers_rmsd_el, free_energy_land_rmsd_el, '-')  # Use axes[2] to plot
        axes[2].set_xlabel('RMSD (Å)')
        axes[2].set_ylabel('Free Energy (kJ/mol)')
        axes[2].set_title('Energy Landscape (RMSD)')

        output_data_filename_rmsd_el = os.path.join(output_dir, f"{output_filename_prefix}_rmsd_energy_landscape_data.txt")
        np.savetxt(output_data_filename_rmsd_el, np.column_stack([bin_centers_rmsd_el, free_energy_land_rmsd_el]),
                    header='RMSD (Å)  Free Energy (kJ/mol)', fmt='%10.5f')
        print(f"RMSD-based energy landscape data saved to file: {output_data_filename_rmsd_el}")
        print(f"RMSD-based energy landscape plot is ready")

    except Exception as e:
        print(f"RMSD-based Energy Landscape calculation failed: {e}")


    # 4. Calculate Rg-based Energy Landscape and plot
    print("Start calculating Rg-based Energy Landscape...")
    try:
        rg_values = []
        for ts in u.trajectory:
            rg = nucleic.radius_of_gyration()
            rg_values.append(rg)
        rg_values = np.array(rg_values)

        hist_rg_el, bin_edges_rg_el = np.histogram(rg_values, bins=50, density=True)
        bin_centers_rg_el = (bin_edges_rg_el[:-1] + bin_edges_rg_el[1:]) / 2
        free_energy_rg_el = -kT * np.log(hist_rg_el)
        free_energy_kJ_mol_rg_el = free_energy_rg_el / 1000 / 4.184
        min_free_energy_rg_el = np.min(free_energy_kJ_mol_rg_el[np.isfinite(free_energy_kJ_mol_rg_el)])
        free_energy_land_rg_el = free_energy_kJ_mol_rg_el - min_free_energy_rg_el

        axes[3].plot(bin_centers_rg_el, free_energy_land_rg_el, '-')  # Use axes[3] to plot
        axes[3].set_xlabel('Radius of Gyration (Rg, Å)')
        axes[3].set_ylabel('Free Energy (kJ/mol)')
        axes[3].set_title('Energy Landscape (Rg)')

        output_data_filename_rg_el = os.path.join(output_dir, f"{output_filename_prefix}_rg_energy_landscape_data.txt")
        np.savetxt(output_data_filename_rg_el, np.column_stack([bin_centers_rg_el, free_energy_land_rg_el]),
                    header='Rg (Å)  Free Energy (kJ/mol)', fmt='%10.5f')
        print(f"Rg-based energy landscape data saved to file: {output_data_filename_rg_el}")
        print(f"Rg-based energy landscape plot is ready")

    except Exception as e:
        print(f"Rg-based Energy Landscape calculation failed: {e}")


    plt.tight_layout()  # Automatically adjust subplot layout to prevent overlap
    output_plot_filename_all = os.path.join(output_dir, f"{output_filename_prefix}_all_plots.png")
    fig.savefig(output_plot_filename_all)  # Save figure containing all subplots
    print(f"All analysis plots saved as: {output_plot_filename_all}")


    # Copy the generated all_plots image to the specified image output folder
    try:
        md_dirname = os.path.basename(output_dir) # Get md_dir directory name
        image_output_path_all_plots = os.path.join(image_output_dir, f"{md_dirname}_all_plots.png") # Use md_dirname as new filename
        fig.savefig(image_output_path_all_plots) # Save all plots image
        print(f"Analysis images copied to image output folder: {image_output_dir}")
    except Exception as e:
        print(f"Failed to copy analysis images to output folder: {e}")


    # plt.show()  #  Display all plots together # Comment out, no need to display plots

    print(f"Trajectory analysis completed. Result files saved to: {output_dir}, images saved to: {image_output_dir}")



def main():
    # --- List of md_dir to be processed ---
    # Here you need to modify the md_dir list according to your actual situation.
    # It can be a manually specified list, or automatically obtain all subdirectories that meet the conditions under a certain directory.
    image_output_dir = os.path.join("path/to/image/output/parent/directory", IMAGE_OUTPUT_DIR_NAME)# Image output folder name
    md_parent_dir = "path/to/parent/directory/containing/md_directories" # Please replace with the parent directory containing all md_dir
    # md_dirs = [os.path.join(md_parent_dir, d) for d in os.listdir(md_parent_dir) if os.path.isdir(os.path.join(md_parent_dir, d)) and d.startswith('iter2_')] # Automatically find directories starting with 'hdockoutpdb_relax_from_index'
    md_dirs = []
    parent_dir_contents = os.listdir(md_parent_dir)

    for item_name in parent_dir_contents:
        item_path = os.path.join(md_parent_dir, item_name)
        if os.path.isdir(item_path):
            if item_name.startswith('iter2_'):
                md_dirs.append(item_path)


    # If you want to manually specify the md_dir list, you can use the following way, and comment out the above automatic search way
    # md_dirs = [
    #     "/path/to/md_dir1",
    #     "/path/to/md_dir2",
    #     # ... more md_dir paths
    # ]
    # --- End of md_dir list configuration ---


    if not os.path.exists(image_output_dir):
        os.makedirs(image_output_dir)  # If image output folder does not exist, create it

    for md_dir in md_dirs:
        if not os.path.isdir(md_dir): # Ensure it is directory
            print(f"Skipping non-directory path: {md_dir}")
            continue

        print(f"Start processing directory: {md_dir}")
        simulation_state_pkl = os.path.join(md_dir, "simulation_state.pkl") # Check if simulation_state.pkl exists
        pro_nucle_gro_fit = os.path.join(md_dir, "protein_nucleic_fit.gro") # Check if PD_fit.gro exists
        pro_nucle_xtc_fit = os.path.join(md_dir, "protein_nucleic_fit.xtc") # Check if PD_fit.xtc exists
        output_dir = md_dir # Output file path is md_dir itself

        if os.path.exists(pro_nucle_gro_fit) and os.path.exists(pro_nucle_xtc_fit): # If fit gro and xtc files already exist, directly perform analysis
            print(f"Found existing fitted files: {pro_nucle_gro_fit}, {pro_nucle_xtc_fit}. Skipping trajectory processing, proceeding directly to analysis.")
            analyze_trajectory(pro_nucle_gro_fit, pro_nucle_xtc_fit, output_dir, OUTPUT_FILENAME_PREFIX, image_output_dir)
        elif os.path.exists(simulation_state_pkl): # If simulation_state.pkl exists, but fit files do not exist, then perform trajectory processing and analysis
            print(f"Found simulation_state.pkl file: {simulation_state_pkl}. Performing trajectory processing and analysis.")
            process_trajectory(md_dir) # Perform trajectory processing
            analyze_trajectory(pro_nucle_gro_fit, pro_nucle_xtc_fit, output_dir, OUTPUT_FILENAME_PREFIX, image_output_dir) # Perform analysis
        else: # If simulation_state.pkl and fit files do not exist, then skip this directory
            print(f"Did not find simulation_state.pkl or processed files in directory {md_dir}, skipping this directory.")
        print(f"Directory {md_dir} processing completed.\n")

    print("All directories processing completed.")


if __name__ == "__main__":
    main()
