import numpy as np
from MDAnalysis import Universe
import MDAnalysis as mda
import warnings
import os
import json
from .utils import (
    apply_pbc,
    calculate_bounding_box,
    initialize_grid,
    update_defect_matrix,
    compute_pairwise_distances
)

warnings.filterwarnings("ignore")


class PackingDefect2:

    def __init__(self, radii_file):

        with open(radii_file, 'r') as file:
            self.types_radii = json.load(file)

    def read_top(self, resname, topology_file):

        tails = ['C2{}'.format(i) for i in range(2, 23)] + \
                ['C3{}'.format(i) for i in range(2, 23)] + \
                ['H{}{}'.format(i, suffix) for i in range(2, 23) for suffix in ['R', 'S', 'X', 'Y']] + \
                ['H16Z', 'H18T', 'H91', 'H101', 'H18Z', 'H20T']

        TGglyc = ['O11', 'O21', 'O31', 'O12', 'O22', 'O32', 'C1', 'C2', 'C3', 'C11', 'C21', 'C31', 'HA', 'HB', 'HS', 'HX', 'HY']

        COside = ['C{}'.format(i) for i in range(1, 28)] + \
                 ['H{}'.format(i) for i in range(1, 28)] + \
                 ['H{}{}'.format(i, suffix) for i in range(1, 28) for suffix in ['A', 'B']] + \
                 ['H3', 'H19C', 'H18C', 'H6', 'H14', 'H17', 'H25', 'H27C', 'H26C']

        with open(topology_file) as file:
            startread = False
            output = {}
            for line in file:
                if line.startswith('!'):
                    continue
                if line.startswith(f'RESI {resname}'):
                    startread = True
                elif startread and line.startswith('BOND'):
                    break
                elif startread and line.startswith('ATOM'):
                    atom_name, atom_type = line.split()[1:3]
                    acyl = -1
                    if resname != 'TRIO':
                        if atom_name in tails:
                            acyl = 1
                    elif atom_name in TGglyc:
                        acyl = 2
                    else:
                        acyl = 3

                    output[atom_name] = [self.types_radii[atom_type], acyl]

        return output


class PackingDefect2Sequential:

    def __init__(self, atomgroups, radii, prefix='./', leaflet='both', defect_types=None, defect_thresholds=None):

        self.N = 10000 
        self.universe = atomgroups[0].universe
        u = atomgroups[0].universe
        self.dt = u.trajectory[0].dt 
        self.atomgroups = atomgroups
        self.dx = 1 
        self.dy = 1  
        self.radii = radii
        self.prefix = prefix
        self.leaflet = leaflet
        self.protein_atoms = self.universe.select_atoms("protein", updating=True)
        self.bbox_data = {"x_min": 0, "x_max": 0, "y_min": 0, "y_max": 0} 
        self._results = []
        # Allow customizable defect types and thresholds
        self.defect_types = defect_types or ['PLacyl', 'TGglyc', 'TGacyl']
        self.defect_thresholds = defect_thresholds or {t: i + 1 for i, t in enumerate(self.defect_types)}

    def process(self):

        self._results = []  
        for ts in self.universe.trajectory:
            print(f"Processing frame {ts.frame}, time: {ts.time:.3f}, pbc: {ts.dimensions[:3]}")
            frame_result = self._single_frame(ts)  
            if frame_result:  
                self._results.append([frame_result])  # Append as a list containing the tuple
        if not self._results:
            print("No frames were processed. Exiting...")
            return
        self._conclude()

    def _single_frame(self, ts):

        ag = self.atomgroups[0]
        protein_atoms = ag.universe.select_atoms("protein", updating=True)
        dim = ts.dimensions.copy()
        pbc = dim[:3]
        print(f"time: {ts.time / 1000:.3f}    pbc: {pbc[0]:.3f} {pbc[1]:.3f} {pbc[2]:.3f}")

        # Apply PBC using the generalized utility
        ag.universe.atoms.positions = apply_pbc(ag.universe.atoms.positions, pbc)

        hz = np.average(ag.select_atoms('name P').positions[:, 2])
        self.bbox_data = calculate_bounding_box(protein_atoms)
        grid = initialize_grid(pbc, self.dx, self.dy, hz)
        zlim, PL = self._analyze_leaflets(ag, hz, grid)

        return grid['up'], grid['dw'], PL['up'] + 5, PL['dw'] - 5, dim

    def _analyze_leaflets(self, ag, hz, grid):

        zlim = {'up': np.max(ag.positions[:, 2]), 'dw': np.min(ag.positions[:, 2])}
        PL = {
            'up': ag.select_atoms(f'name P and prop z > {hz}').center_of_mass()[2],
            'dw': ag.select_atoms(f'name P and prop z < {hz}').center_of_mass()[2],
        }
        atoms = {}
        if self.leaflet in ['both', 'up']:
            atoms['up'] = ag.select_atoms(f'prop z > {PL["up"] - 20}')
        if self.leaflet in ['both', 'dw']:
            atoms['dw'] = ag.select_atoms(f'prop z < {PL["dw"] + 20}')
        for leaflet, atom_group in atoms.items():
            for atom in atom_group:
                xatom, yatom, zatom = atom.position
                radius, acyl = self.radii[atom.resname][atom.name]
                update_defect_matrix(
                    grid, xatom, yatom, zatom, radius, acyl, leaflet, self.dx, self.dy
                )

        return zlim, PL

    def _conclude(self):

        print("Concluding...")
        Mup, Mdw, zlimup, zlimdw, dim = self._aggregate_results()
        df = self._initialize_defect_universe(dim)
        self._process_defects(df, Mup, Mdw, zlimup, zlimdw, dim)

    def _aggregate_results(self):

        Mup, Mdw, zlimup, zlimdw, dim = [], [], [], [], []
        for r in self._results:
            for rr in r:
                if rr[0] is None:
                    continue
                Mup.append(rr[0])
                Mdw.append(rr[1])
                zlimup.append(rr[2])
                zlimdw.append(rr[3])
                dim.append(rr[4])

        return Mup, Mdw, zlimup, zlimdw, dim

    def _initialize_defect_universe(self, dim):
        nframes = len(dim)
        fac = np.zeros((nframes, self.N, 3))

        print(f"Number of frames: {nframes}")
        print(f"Number of atoms: {self.N}")
        print(f"Shape of fac array: {fac.shape}")

        df = Universe.empty(
            n_atoms=self.N,
            n_residues=self.N,
            atom_resindex=np.arange(self.N),
            residue_segindex=[0] * self.N,
            trajectory=True,
        )

        df.add_TopologyAttr('resname', ['O'] * self.N)
        df.add_TopologyAttr('name', ['O'] * self.N)
        df.add_TopologyAttr('resid', np.arange(self.N) + 1)

        df.load_new(fac, order='fac')
        df.trajectory[0].dt = self.dt

        for i, ts in enumerate(df.trajectory):
            df.trajectory[i].dimensions = dim[i]

        return df

    def _process_defects(self, df, Mup, Mdw, zlimup, zlimdw, dim):
        defect_uni = {d: df.copy() for d in self.defect_types}

        for d in self.defect_types:
            threshold = self.defect_thresholds[d]
            for i, ts in enumerate(defect_uni[d].trajectory):
                num = 0
                num = self._process_leaflet(threshold, Mup[i], zlimup[i], defect_uni[d], num)
                self._process_leaflet(threshold, Mdw[i], zlimdw[i], defect_uni[d], num)
        self._save_outputs(defect_uni, self.defect_types)

    def _process_leaflet(self, threshold, M, zlim, defect_universe, num):
        bA = (M == threshold)
        ind = np.where(bA)
        xs, ys = ind[1], ind[0]

        for x1, y1 in zip(xs, ys):
            if num >= self.N:
                break
            pos = np.array([x1, y1, zlim])
            defect_universe.atoms[num].position = pos
            num += 1

        return num

    def _save_outputs(self, defect_uni, defects):
        output_base_dir = 'GRO_paper'

        for d in defects:
            output_dir = os.path.join(self.prefix, d)
            os.makedirs(output_dir, exist_ok=True)

            u = defect_uni[d]

            for i, ts in enumerate(u.trajectory):
                self.protein_atoms.universe.trajectory[i]
                combined_universe = mda.Merge(self.protein_atoms, u.atoms)

                combined_universe.atoms.positions[len(self.protein_atoms):] = ts.positions
                combined_universe.atoms.positions[:len(self.protein_atoms)] = self.protein_atoms.positions
                combined_universe.trajectory.ts.dimensions = ts.dimensions

                output_filepath = os.path.join(output_dir, f"{d}_frame_{i}.gro")
                combined_universe.atoms.write(output_filepath)
