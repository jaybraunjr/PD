import numpy as np
import MDAnalysis as mda
from MDAnalysis import Universe
import numpy as np

def apply_pbc(positions, box):
    """Apply periodic boundary conditions to positions."""
    box_xy = np.array([box[0], box[1], 0])
    box_xyz = np.array([box[0], box[1], box[2]])
    return positions - box_xy * np.floor(positions / box_xyz)


def calculate_bounding_box(atomgroup, padding=10):
    """Calculate the bounding box of an atom group with padding."""
    bbox = atomgroup.bbox()
    x_min, y_min, _ = bbox[0]
    x_max, y_max, _ = bbox[1]
    return {
        "x_min": x_min - padding,
        "x_max": x_max + padding,
        "y_min": y_min - padding,
        "y_max": y_max + padding,
    }


def initialize_grid(box, dx=1, dy=1, hz=0):
    """Initialize a 2D grid based on box dimensions."""
    xarray = np.arange(0, box[0], dx)
    yarray = np.arange(0, box[1], dy)
    xx, yy = np.meshgrid(xarray, yarray)
    return {
        'xx': xx,
        'yy': yy,
        'up': np.zeros_like(xx),
        'dw': np.zeros_like(xx),
        'z_up': np.zeros_like(xx) + hz,
        'z_dw': np.zeros_like(xx) + hz,
    }


def update_defect_matrix(grid, xatom, yatom, zatom, radius, acyl, leaflet, dx, dy):
    """Update defect matrices for a single atom."""
    dxx = grid['xx'] - xatom
    dxx -= grid['xx'].shape[0] * np.round(dxx / grid['xx'].shape[0])
    dyy = grid['yy'] - yatom
    dyy -= grid['yy'].shape[1] * np.round(dyy / grid['yy'].shape[1])

    dist_meet = (np.sqrt(dx**2 + dy**2) / 2 + radius)**2
    bAr = dxx**2 + dyy**2 < dist_meet

    if leaflet == 'up':
        bAnP = grid['up'] >= 0
        baZ = zatom > grid['z_up']
        bA = bAr & bAnP & baZ
        grid['up'][bA] = acyl
        grid['z_up'][bA] = zatom
    elif leaflet == 'dw':
        bAnP = grid['dw'] >= 0
        baZ = zatom < grid['z_dw']
        bA = bAr & bAnP & baZ
        grid['dw'][bA] = acyl
        grid['z_dw'][bA] = zatom


def compute_pairwise_distances(positions1, positions2):
    """Compute pairwise distances between two sets of positions."""
    diff = positions1[:, np.newaxis, :] - positions2
    return np.sqrt(np.sum(diff**2, axis=2))


def _dfs(graph, start):
    """Depth First Search for connected components in a graph."""
    visited, stack = set(), [start]
    while stack:
        vertex = stack.pop()
        if vertex not in visited:
            visited.add(vertex)
            stack.extend(graph[vertex] - visited)
    return visited


def _make_graph(matrix):
    """Convert a binary matrix to a graph representation."""
    graph = {}
    xis, yis = matrix.shape  # Get the dimensions of the matrix

    for (xi, yi), value in np.ndenumerate(matrix):
        if value == 0:
            continue  

        node_index = xi * yis + yi
        neighbor_list = []  
        for dx in [-1, 0, 1]:
            for dy in [-1, 0, 1]:

                x, y = divmod(xi + dx, xis)[1], divmod(yi + dy, yis)[1]

                if matrix[x, y] == 1 and (x, y) != (xi, yi):
                    neighbor_node_index = x * yis + y  
                    neighbor_list.append(neighbor_node_index)  

        graph[node_index] = set(neighbor_list)  

    return graph


def filter_defects_by_distance(defects, positions, threshold):
    """Filter defects based on distance to a set of positions."""
    distances = compute_pairwise_distances(defects, positions)
    return [defect for i, defect in enumerate(defects) if np.all(distances[i] > threshold)]


def write_to_file(output_file, data):
    """Write data to a text file."""
    with open(output_file, 'w') as f:
        for line in data:
            f.write(f"{line}\n")

def validate_defect_thresholds(defect_types, defect_thresholds):
    for dt in defect_types:
        if dt not in defect_thresholds:
            raise ValueError(f"Missing threshold for defect type: {dt}")


def populate_grid_with_atoms(grid, atom_group, radii_lookup, leaflet, dx, dy):
    for atom in atom_group:
        x, y, z = atom.position
        try:
            radius, acyl = radii_lookup[atom.resname][atom.name]
        except KeyError:
            continue  # skip atoms not present in radii lookup
        update_defect_matrix(grid, x, y, z, radius, acyl, leaflet, dx, dy)

def get_defect_coordinates(grid, threshold):
    mask = (grid == threshold)
    ys, xs = np.where(mask)
    return xs, ys

def write_combined_gro(protein_atoms, defect_atoms, dimensions, filepath):
    combined = mda.Merge(protein_atoms, defect_atoms)
    combined.atoms.positions[:len(protein_atoms)] = protein_atoms.positions
    combined.atoms.positions[len(protein_atoms):] = defect_atoms.positions
    combined.trajectory.ts.dimensions = dimensions
    combined.atoms.write(filepath)



def initialize_empty_defect_universe(n_atoms, nframes, dims, dt):
    fac = np.zeros((nframes, n_atoms, 3))

    df = Universe.empty(
        n_atoms=n_atoms,
        n_residues=n_atoms,
        atom_resindex=np.arange(n_atoms),
        residue_segindex=[0] * n_atoms,
        trajectory=True,
    )
    df.add_TopologyAttr('resname', ['O'] * n_atoms)
    df.add_TopologyAttr('name', ['O'] * n_atoms)
    df.add_TopologyAttr('resid', np.arange(n_atoms) + 1)
    df.load_new(fac, order='fac')
    df.trajectory[0].dt = dt

    for i, ts in enumerate(df.trajectory):
        df.trajectory[i].dimensions = dims[i]

    return df
