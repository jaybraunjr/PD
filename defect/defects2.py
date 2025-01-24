import numpy as np
from MDAnalysis import Universe
import MDAnalysis as mda
import warnings
import os
import json

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
        self._apply_pbc(ag, pbc)
        hz = np.average(ag.select_atoms('name P').positions[:, 2])
        self._calculate_bounding_box(protein_atoms)
        xx, yy, M, Z = self._initialize_grid(pbc, hz)
        zlim, PL = self._analyze_leaflets(ag, hz, xx, yy, M, Z)

        return M['up'], M['dw'], PL['up'] + 5, PL['dw'] - 5, dim

    # --- Helper Functions ---
    def _apply_pbc(self, ag, pbc):

        pbc_xy0 = np.array([pbc[0], pbc[1], 0])
        pbc_xyz = np.array([pbc[0], pbc[1], pbc[2]])
        ag.universe.atoms.positions -= pbc_xy0 * np.floor(ag.universe.atoms.positions / pbc_xyz)

    def _calculate_bounding_box(self, protein_atoms):

        protein_bbox = protein_atoms.bbox()
        x_min, y_min, _ = protein_bbox[0]
        x_max, y_max, _ = protein_bbox[1]
        x_min -= 10
        x_max += 10
        y_min -= 10
        y_max += 10

        self.bbox_data = {
            "x_min": x_min,
            "x_max": x_max,
            "y_min": y_min,
            "y_max": y_max,
        }

    def _initialize_grid(self, pbc, hz):

        xarray = np.arange(0, pbc[0], self.dx)
        yarray = np.arange(0, pbc[1], self.dy)
        xx, yy = np.meshgrid(xarray, yarray)
        M = {'up': np.zeros_like(xx), 'dw': np.zeros_like(xx)}
        Z = {'up': np.zeros_like(xx) + hz, 'dw': np.zeros_like(xx) + hz}

        return xx, yy, M, Z

    def _analyze_leaflets(self, ag, hz, xx, yy, M, Z):

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
                self._update_defect_matrix(
                    leaflet, xatom, yatom, zatom, radius, acyl, xx, yy, M, Z, zlim
                )

        return zlim, PL

    def _update_defect_matrix(self, leaflet, xatom, yatom, zatom, radius, acyl, xx, yy, M, Z, zlim):

        dxx = xx - xatom
        dxx -= xx.shape[0] * np.round(dxx / xx.shape[0])
        dyy = yy - yatom
        dyy -= yy.shape[1] * np.round(dyy / yy.shape[1])

        dist_meet = (np.sqrt(self.dx**2 + self.dy**2) / 2 + radius)**2
        bAr = dxx**2 + dyy**2 < dist_meet

        if acyl == -1:
            M[leaflet][bAr] = acyl
            return

        bAnP = M[leaflet] >= 0
        baZ = zatom > Z[leaflet] if leaflet == 'up' else zatom < Z[leaflet]
        bA = bAr & bAnP & baZ
        M[leaflet][bA] = acyl
        Z[leaflet][bA] = zatom


    def _conclude(self):

        print("Concluding...")
        Mup, Mdw, zlimup, zlimdw, dim = self._aggregate_results()
        # self.N = self._calculate_max_defects(Mup, Mdw)
        df = self._initialize_defect_universe(dim)
        defects = ['PLacyl', 'TGglyc', 'TGacyl']
        defect_thr = {'PLacyl': 1, 'TGglyc': 2, 'TGacyl': 3}
        self._process_defects(df, Mup, Mdw, zlimup, zlimdw, dim, defects, defect_thr)

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


    def _calculate_max_defects(self, Mup, Mdw):
        max_up = max(np.sum(m > 0) for m in Mup) if Mup else 0
        max_dw = max(np.sum(m > 0) for m in Mdw) if Mdw else 0
        return max(max_up, max_dw)

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


    def _process_defects(self, df, Mup, Mdw, zlimup, zlimdw, dim, defects, defect_thr):
        defect_uni = {d: df.copy() for d in self.defect_types}
        defect_clu = {d: [] for d in self.defect_types}

        for d in self.defect_types:
            threshold = self.defect_thresholds[d]
            for i, ts in enumerate(defect_uni[d].trajectory):
                num = 0
                num = self._process_leaflet(
                    defect_thr[d], Mup[i], zlimup[i], defect_uni[d], num
                )
                self._process_leaflet(
                    defect_thr[d], Mdw[i], zlimdw[i], defect_uni[d], num
                )
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

                # Update positions
                combined_universe.atoms.positions[len(self.protein_atoms):] = ts.positions
                combined_universe.atoms.positions[:len(self.protein_atoms)] = self.protein_atoms.positions
                combined_universe.trajectory.ts.dimensions = ts.dimensions

                # Save to file
                output_filepath = os.path.join(output_dir, f"{d}_frame_{i}.gro")
                combined_universe.atoms.write(output_filepath)




