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
    compute_pairwise_distances,
    validate_defect_thresholds,
    populate_grid_with_atoms,
    initialize_empty_defect_universe,
    get_defect_coordinates,
    write_combined_gro
)

warnings.filterwarnings("ignore")


class PackingDefect2:
    def __init__(self, radii_file, classification_rules=None):
        with open(radii_file, 'r') as file:
            self.types_radii = json.load(file)
        self.classification_rules = classification_rules

    def default_classify(self, resname, atom_name):
        tails = [f'C2{i}' for i in range(2, 23)] + \
                [f'C3{i}' for i in range(2, 23)] + \
                [f'H{i}{s}' for i in range(2, 23) for s in ['R', 'S', 'X', 'Y']] + \
                ['H16Z', 'H18T', 'H91', 'H101', 'H18Z', 'H20T']
        TGglyc = ['O11', 'O21', 'O31', 'O12', 'O22', 'O32', 'C1', 'C2', 'C3',
                  'C11', 'C21', 'C31', 'HA', 'HB', 'HS', 'HX', 'HY']

        if resname != 'TRIO':
            return 1 if atom_name in tails else -1
        else:
            return 2 if atom_name in TGglyc else 3

    def read_top(self, resname, topology_file):
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
                    if callable(self.classification_rules):
                        acyl = self.classification_rules(resname, atom_name)
                    elif isinstance(self.classification_rules, dict):
                        acyl = self.classification_rules.get(atom_name, -1)
                    else:
                        acyl = self.default_classify(resname, atom_name)
                    output[atom_name] = [self.types_radii[atom_type], acyl]

        return output


class PackingDefect2Sequential:
    def __init__(self, atomgroups, radii, prefix='./', leaflet='both', defect_types=None, defect_thresholds=None):
        self.N = 10000 
        self.universe = atomgroups[0].universe
        self.dt = self.universe.trajectory[0].dt 
        self.atomgroups = atomgroups
        self.dx = 1 
        self.dy = 1  
        self.radii = radii
        self.prefix = prefix
        self.leaflet = leaflet
        self.protein_atoms = self.universe.select_atoms("protein", updating=True)
        self.bbox_data = {"x_min": 0, "x_max": 0, "y_min": 0, "y_max": 0} 
        self._results = []

        self.defect_types = defect_types or ['PLacyl', 'TGglyc', 'TGacyl']
        self.defect_thresholds = defect_thresholds or {t: i + 1 for i, t in enumerate(self.defect_types)}
        validate_defect_thresholds(self.defect_types, self.defect_thresholds)

    def process(self):
        self._results = []  
        for ts in self.universe.trajectory:
            print(f"Processing frame {ts.frame}, time: {ts.time:.3f}, pbc: {ts.dimensions[:3]}")
            frame_result = self._single_frame(ts)  
            if frame_result:  
                self._results.append([frame_result])
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
            populate_grid_with_atoms(grid, atom_group, self.radii, leaflet, self.dx, self.dy)

        return zlim, PL

    def _conclude(self):
        print("Concluding...")
        Mup, Mdw, zlimup, zlimdw, dim = self._aggregate_results()
        df = initialize_empty_defect_universe(self.N, len(dim), dim, self.dt)
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
        xs, ys = get_defect_coordinates(M, threshold)

        for x1, y1 in zip(xs, ys):
            if num >= self.N:
                break
            pos = np.array([x1, y1, zlim])
            defect_universe.atoms[num].position = pos
            num += 1

        return num

    def _save_outputs(self, defect_uni, defects):
        for d in defects:
            output_dir = os.path.join(self.prefix, d)
            os.makedirs(output_dir, exist_ok=True)
            u = defect_uni[d]

            for i, ts in enumerate(u.trajectory):
                self.protein_atoms.universe.trajectory[i]
                output_filepath = os.path.join(output_dir, f"{d}_frame_{i}.gro")
                write_combined_gro(self.protein_atoms, u.atoms, ts.dimensions, output_filepath)

