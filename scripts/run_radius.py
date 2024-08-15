from ../defect/radius import process_frames, renumber_all_gro_files, calculate_defects
import matplotlib.pyplot as plt
import MDAnalysis as mda
import os
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from defect.radius import PackingDefect2, PackingDefect2PMDA

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
    base_directory = 'GRO_paper/rep3_skip100/'
    lipid_types = ['PLacyl', 'TGacyl', 'TGglyc']
    output_base_dir = 'filtered_gro_paper'

    for lipid_type in lipid_types:
        directory_prefix = os.path.join(base_directory, f"{lipid_type}")
        output_dir = os.path.join(output_base_dir, lipid_type)

        print(f"Processing {lipid_type}...")
        output_files = process_frames(frame_start, frame_end, protein_atom_count, directory_prefix, lipid_type, output_dir)
        print(f"Renumbering {lipid_type} files...")
        renumbered_files = renumber_all_gro_files(output_files)
        print(f"Completed processing for {lipid_type}.")

    processed_defects_up = {}
    processed_defects_down = {}

    for lipid_type in lipid_types:
        processed_directory_prefix = os.path.join(output_base_dir, lipid_type)
        defects_up_all = []
        defects_down_all = []

        for frame_idx in range(frame_start, frame_end + 1):
            processed_file_name = f"renumbered_{lipid_type}_corrected_frame_{frame_idx}.gro"
            processed_file_path = os.path.join(processed_directory_prefix, processed_file_name)

            if os.path.exists(processed_file_path):
                u = mda.Universe(processed_file_path)
                defects_up, defects_down = calculate_defects(u)
                defects_up_all.extend(defects_up)
                defects_down_all.extend(defects_down)

        processed_defects_up[lipid_type] = defects_up_all
        processed_defects_down[lipid_type] = defects_down_all


    plot_histogram(processed_defects_up['PLacyl'], label='PLacyl', color='red')
    plot_histogram(processed_defects_down['PLacyl'], label='PLacyl with protein', color='darkred')
    plot_histogram(processed_defects_up['TGacyl'], label='TGacyl', color='blue')
    plot_histogram(processed_defects_down['TGacyl'], label='TGacyl with protein', color='darkblue')
    plot_histogram(processed_defects_up['TGglyc'], label='TGglyc', color='green')
    plot_histogram(processed_defects_down['TGglyc'], label='TGglyc with protein', color='darkgreen')

    plt.legend()
    plt.show()