import matplotlib.pyplot as plt
import MDAnalysis as mda
import matplotlib
matplotlib.use('Agg')  
import os
import sys
import numpy as np
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from defect import radius2
import numpy as np
import os
os.chdir(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

def plot_histogram(defects, label, color=None):
    h, _ = np.histogram(defects, bins=np.linspace(0, 150, 600))
    h[0] = 0  # Set the first bin to 0 for better visualization
    binp = 0.5 * (_[1:] + _[:-1])  # Midpoints of bins
    plt.scatter(binp, h / np.sum(h), label=label, color=color)
    plt.yscale('log')

if __name__ == "__main__":
    frame_start = 0
    frame_end = 10
    protein_atom_count = 626
    base_directory = 'output/rep3_skip100/'
    lipid_types = ['PLacyl', 'TGacyl', 'TGglyc']
    output_base_dir = 'filtered_gro_paper'
    apply_protein_cutoff = False

    for lipid_type in lipid_types:
        directory_prefix = os.path.join(base_directory, f"{lipid_type}")
        output_dir = os.path.join(output_base_dir, lipid_type)

        print(f"Processing {lipid_type}...")
        output_files = radius2.process_frames(frame_start, frame_end, protein_atom_count, directory_prefix, lipid_type, output_dir, use_cutoff=apply_protein_cutoff)
        print(f"Renumbering {lipid_type} files...")
        renumbered_files = radius2.renumber_all_gro_files(output_files)
        print(f"Completed processing for {lipid_type}.")

    # Store all defect data
    processed_defects_up = {}
    processed_defects_down = {}
    processed_defects_combined = {}

    for lipid_type in lipid_types:
        processed_directory_prefix = os.path.join(output_base_dir, lipid_type)
        defects_up_all = []
        defects_down_all = []
        defects_combined_all = []

        for frame_idx in range(frame_start, frame_end + 1):
            processed_file_name = f"renumbered_{lipid_type}_corrected_frame_{frame_idx}.gro"
            processed_file_path = os.path.join(processed_directory_prefix, processed_file_name)

            if os.path.exists(processed_file_path):
                u = mda.Universe(processed_file_path)
                defects_up, defects_down = radius2.calculate_defects(u, combine=False)
                defects_up_all.extend(defects_up)
                defects_down_all.extend(defects_down)
                defects_combined_all.extend(defects_up + defects_down)  # Combine here

        processed_defects_up[lipid_type] = defects_up_all
        processed_defects_down[lipid_type] = defects_down_all
        processed_defects_combined[lipid_type] = defects_combined_all

    # Plot top and bottom defects together
    plt.figure(figsize=(8, 6))
    for lipid_type, color in zip(lipid_types, ['red', 'blue', 'green']):
        plot_histogram(processed_defects_up[lipid_type], label=f'{lipid_type} Up', color=color)
        plot_histogram(processed_defects_down[lipid_type], label=f'{lipid_type} protein', color=f'dark{color}')
    plt.legend()
    plt.title('Top and Bottom Defects')
    plt.xlabel('Defect Size')
    plt.ylabel('Frequency (log scale)')
    plt.savefig('top_bottom_defect_histogram2.png', dpi=300, bbox_inches='tight')
    plt.close()

    # Plot combined defects together
    plt.figure(figsize=(8, 6))
    for lipid_type, color in zip(lipid_types, ['red', 'blue', 'green']):
        plot_histogram(processed_defects_combined[lipid_type], label=f'{lipid_type} Combined', color=color)
    plt.legend()
    plt.title('Combined Defects')
    plt.xlabel('Defect Size')
    plt.ylabel('Frequency (log scale)')
    plt.savefig('combined_defect_histogram2.png', dpi=300, bbox_inches='tight')
    plt.close()
