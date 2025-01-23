import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt
import os
from MDAnalysis.lib import distances
import subprocess
from .utils import _make_graph, _dfs



def calculate_defects(u):
    defects = []
    defects_up = []
    defects_down = []
    protein_atoms = u.select_atoms("protein")

    for ts in u.trajectory:
        ag = u.select_atoms('prop z > 0')
        hz = np.average(ag.positions[:,2])
        agup = u.select_atoms('prop z > %f' %hz)
        agdw = u.select_atoms('prop z < %f' %hz)

        xarray = np.arange(0, u.dimensions[0], 1)
        yarray = np.arange(0, u.dimensions[1], 1)
        xx, yy = np.meshgrid(xarray, yarray)
        Mup = np.zeros_like(xx)
        Mdw = np.zeros_like(xx)
        
        ### UP
        xind = np.minimum(agup.positions[:,0].astype(np.int64), Mup.shape[0]-1)
        yind = np.minimum(agup.positions[:,1].astype(np.int64), Mup.shape[1]-1)
        Mup[xind, yind] = 1

                    
        graph = _make_graph(Mup)
        visited = set([])
        for n in graph:
            if n not in visited:
                defect_loc = _dfs(graph, n)
                visited = visited.union(defect_loc)
                defects_up.append(len(defect_loc))
        
        ### DW  
        xind = np.minimum(agdw.positions[:,0].astype(np.int64), Mdw.shape[0]-1)
        yind = np.minimum(agdw.positions[:,1].astype(np.int64), Mdw.shape[1]-1)
        Mdw[xind, yind] = 1

        graph = _make_graph(Mdw)
        visited = set([])
        for n in graph:
            if n not in visited:
                defect_loc = _dfs(graph, n)
                visited = visited.union(defect_loc)
                defects_down.append(len(defect_loc))

    return defects_up, defects_down

def calculate_defects_from_gro(directory_prefix, file_prefix, start_frame=0):
    defects_up = []
    defects_down = []
    print(f'Working in directory: {directory_prefix}')

    frame_idx = start_frame
    while True:
        gro_file_path = os.path.join(directory_prefix, f"{file_prefix}_frame_{frame_idx}.gro")
        print(f'Checking for file: {gro_file_path}')  # Debugging print

        if not os.path.exists(gro_file_path):
            print(f'File not found: {gro_file_path}')
            break

        u = mda.Universe(gro_file_path)
        up, down = calculate_defects(u)
        defects_up.extend(up)
        defects_down.extend(down)

        frame_idx += 1

    return defects_up, defects_down

def process_directories(base_directory):
    suffix_pairs = {
        'resultsPLacyl': 'PLacyl',
        'resultsTGacyl': 'TGacyl',
        'resultsTGglyc': 'TGglyc'
    }
    
    all_defects_up = {}
    all_defects_down = {}
    
    for dir_suffix, file_suffix in suffix_pairs.items():
        directory_prefix = os.path.join(base_directory, dir_suffix)
        defects_up, defects_down = calculate_defects_from_gro(directory_prefix, file_suffix)
        all_defects_up[dir_suffix] = defects_up
        all_defects_down[dir_suffix] = defects_down

    return all_defects_up, all_defects_down

def compute_distances(positions1, positions2):
    diff = positions1[:, np.newaxis, :] - positions2
    return np.sqrt(np.sum(diff**2, axis=2))

def write_filtered_gro_by_atom_count(input_file, output_file, cutoff_distance=1.5, protein_atom_count=627):
    with open(input_file, 'r') as f:
        lines = f.readlines()

    header = lines[0].strip()
    footer = lines[-1].strip()

    protein_atoms = lines[2:2+protein_atom_count]
    defect_atoms = lines[2+protein_atom_count:-1]

    # Extract positions of protein and defect atoms
    protein_positions = np.array([[float(line[20:28].strip()), float(line[28:36].strip()), float(line[36:44].strip())] for line in protein_atoms])
    defect_positions = np.array([[float(line[20:28].strip()), float(line[28:36].strip()), float(line[36:44].strip())] for line in defect_atoms])
    min_distances = np.min(compute_distances(defect_positions, protein_positions), axis=1)
    
    # Filter defect atoms based on the cutoff distance
    filtered_defect_atoms = [atom for i, atom in enumerate(defect_atoms) if min_distances[i] > cutoff_distance]

    with open(output_file, 'w') as f:
        f.write(header + '\n')
        f.write(f"{len(protein_atoms) + len(filtered_defect_atoms)}\n")
        f.writelines(protein_atoms)
        f.writelines(filtered_defect_atoms)
        f.write(footer)
    

def process_frames(frame_start, frame_end, protein_atom_count, directory_prefix, lipid_type, output_dir, min_cutoff_distance=1.0):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    output_files = []
    for frame_idx in range(frame_start, frame_end + 1):
        input_file_path = f"{directory_prefix}/{lipid_type}_frame_{frame_idx}.gro"
        output_file_path = f"{output_dir}/{lipid_type}_corrected_frame_{frame_idx}.gro"
        write_filtered_gro_by_atom_count(input_file_path, output_file_path, min_cutoff_distance, protein_atom_count)
        output_files.append(output_file_path)
    return output_files

def renumber_gro(input_file, output_file):
    try:
        subprocess.run(['gmx', 'genconf','-f', input_file, '-o', output_file, '-renumber'], check=True)
        print(f'Renumbered gro saved to {output_file}')
    except subprocess.CalledProcessError as e:
        print(f'Error in renumbering: {e}')

def renumber_all_gro_files(input_files):
    renumbered_files = []
    for input_file in input_files:
        output_file = os.path.join(os.path.dirname(input_file), "renumbered_" + os.path.basename(input_file))
        renumber_gro(input_file, output_file)
        renumbered_files.append(output_file)
        try:
            os.remove(input_file)
            print(f"Deleted file: {input_file}")
        except OSError as e:
            print(f"Error deleting file {input_file}: {e}")
    return renumbered_files