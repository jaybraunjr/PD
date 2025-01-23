from pmda.parallel import ParallelAnalysisBase
import numpy as np
from MDAnalysis import Universe
import MDAnalysis as mda
import warnings
import os
from .utils import _make_graph, _dfs
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
                    break  # Exit the loop once we reach 'BOND'
                elif startread and line.startswith('ATOM'):
                    atom_name, atom_type = line.split()[1:3]

                    acyl = -1  # Default value
                    if resname != 'TRIO':
                        if atom_name in tails:
                            acyl = 1
                    elif atom_name in TGglyc:
                        acyl = 2
                    else:
                        acyl = 3

                    output[atom_name] = [self.types_radii[atom_type], acyl]

        return output
    

    def defect_size(self, matrices, nbins, bin_max, prob=True):

        bins = np.linspace(0, bin_max, nbins)

        defects = []
        for matrix in matrices:
            graph = _make_graph(matrix)
            visited = set()
            for n in graph:
                if n not in visited:
                    defect_loc = _dfs(graph, n)
                    visited.update(defect_loc)
                    defects.append(len(defect_loc))

        hist, bin_edges = np.histogram(defects, bins)
        hist = hist.astype(np.float64)
        binp = 0.5 * (bin_edges[1:] + bin_edges[:-1])

        if hist.sum() == 0:
            return None, None

        if prob:
            hist /= hist.sum()

        return binp, hist


class PackingDefect2PMDA(ParallelAnalysisBase):
    def __init__(self, atomgroups, radii, nbins=600, bin_max=150, prefix='./', prob=True, leaflet='both'):
        u = atomgroups[0].universe
        self.N = 3000  
        self.dt = u.trajectory[0].dt 
        self.dx = 1 
        self.dy = 1  
        self.radii = radii
        self.nbins = nbins
        self.bin_max = bin_max
        self.prefix = prefix
        self.prob = prob
        self.leaflet = leaflet
        self.all_defect_data = {}
        self.protein_atoms = u.select_atoms("protein", updating=True)  
        self.bbox_data = {"x_min": 0, "x_max": 0, "y_min": 0, "y_max": 0}  

        super(PackingDefect2PMDA, self).__init__(u, atomgroups)

    def _prepare(self):
        pass
    
    def _single_frame(self, ts, atomgroups):

        ag = atomgroups[0] 
        protein_atoms = ag.universe.select_atoms("protein", updating=True)  # Select protein atoms
        dim = ts.dimensions.copy()  
        pbc = dim[0:3]  
        print('time: {:.3f}    pbc: {:.3f} {:.3f} {:.3f}'.format(ts.time/1000, pbc[0], pbc[1], pbc[2]))


        pbc_xy0 = np.array([pbc[0], pbc[1], 0])
        pbc_xyz = np.array([pbc[0], pbc[1], pbc[2]])
        aa = ag.universe.atoms
        aa.positions -= pbc_xy0 * np.floor(aa.positions / pbc_xyz)

        hz = np.average(ag.select_atoms('name P').positions[:, 2])

        # Calculate the bounding box for protein atoms
        protein_bbox = protein_atoms.bbox()
        x_min, y_min, _ = protein_bbox[0]
        x_max, y_max, _ = protein_bbox[1]

        x_min -= 10  
        x_max += 10
        y_min -= 10
        y_max += 10

        self.bbox_data["x_min"] = x_min
        self.bbox_data["x_max"] = x_max
        self.bbox_data["y_min"] = y_min
        self.bbox_data["y_max"] = y_max

        xarray = np.arange(0, pbc[0], self.dx)
        yarray = np.arange(0, pbc[1], self.dy)
        xx, yy = np.meshgrid(xarray, yarray)

        M = {'up': np.zeros_like(xx), 'dw': np.zeros_like(xx)}
        Z = {'up': np.zeros_like(xx), 'dw': np.zeros_like(xx)}
        Z['up'] += hz
        Z['dw'] += hz

        zlim = {'up': np.max(ag.positions[:, 2]), 'dw': np.min(ag.positions[:, 2])}

        # Calculate center of mass for phospholipids in each leaflet
        PL = {'up': ag.select_atoms('name P and prop z > %f' % hz).center_of_mass()[2],
              'dw': ag.select_atoms('name P and prop z < %f' % hz).center_of_mass()[2]}

        # Determine which atoms to analyze based on the specified leaflet
        atoms = {}
        if self.leaflet in ['both', 'up']:
            atoms['up'] = ag.select_atoms('prop z > %f' % (PL['up'] - 20))
        if self.leaflet in ['both', 'dw']:
            atoms['dw'] = ag.select_atoms('prop z < %f' % (PL['dw'] + 20))

        # Analyze each atom in the specified leaflets
        for l in atoms:
            for atom in atoms[l]:
                xatom, yatom, zatom = atom.position

                if l == 'up': 
                    assert zatom > PL[l] - 20, 'check Z pos'
                if l == 'dw': 
                    assert zatom < PL[l] + 20, 'check Z pos'

                radius, acyl = self.radii[atom.resname][atom.name]
                dxx = xx - xatom
                dxx -= pbc[0] * np.around(dxx/pbc[0])
                dyy = yy - yatom
                dyy -= pbc[1] * np.around(dyy/pbc[1])

                # Determine if the points in the meshgrid are within the interaction radius
                dist_meet = (np.sqrt(self.dx**2 + self.dy**2)/2 + radius)**2
                bAr = dxx ** 2 + dyy **2 < dist_meet

                # Exclude protein atoms if needed (commented out for now)
                # if atom.resname in self.protein_atoms.resnames:  
                #     continue  

                if acyl == -1:
                    M[l][bAr] = acyl
                    continue
                else:
                    bAnP = M[l] >= 0
                    if l == 'up':
                        baZ = zatom > Z[l]
                    else:
                        baZ = zatom < Z[l]
                    bA = bAr & bAnP & baZ
                    M[l][bA] = acyl
                    Z[l][bA] = zatom


        return M['up'], M['dw'], PL['up']+5, PL['dw']-5, dim
        # return M['up'].tolist(), M['dw'].tolist(), PL['up'] + 5, PL['dw'] - 5, dim



    def _conclude(self):
        print("Concluding...")
        # Initialize lists to store results
        Mup = []; Mdw = []; zlimup = []; zlimdw = []; dim = []

        # Aggregate results from each processed frame
        for r in self._results:
            for rr in r:
                if rr[0] is None:
                    continue
                Mup.append(rr[0])  # Append upper leaflet matrix
                Mdw.append(rr[1])  # Append lower leaflet matrix
                zlimup.append(rr[2])  # Append upper leaflet z-limit
                zlimdw.append(rr[3])  # Append lower leaflet z-limit
                dim.append(rr[4])  # Append frame dimensions

        # Set up a new Universe for storing defect information
        N = self.N  # The maximum number of defects
        df = Universe.empty(n_atoms=N,
                            n_residues=N,
                            atom_resindex=np.arange(N),
                            residue_segindex=[0] * N,
                            trajectory=True)

        df.add_TopologyAttr('resname', ['O'] * N)
        df.add_TopologyAttr('name', ['O'] * N)
        df.add_TopologyAttr('resid', np.arange(N) + 1)

        # Initialize the trajectory data for the new Universe
        nframes = len(dim)
        fac = np.zeros((nframes, N, 3))
        df.load_new(fac, order='fac')
        df.trajectory[0].dt = self.dt 

        for i, ts in enumerate(df.trajectory):
            df.trajectory[i].dimensions = dim[i]

        # Define defect types
        defects = ['PLacyl', 'TGglyc', 'TGacyl']

        # Initialize dictionaries to store defect universe and cluster data
        defect_uni = {}
        defect_clu = {}
        for d in defects:
            defect_uni[d] = df.copy()  # Create a copy of the universe for each defect type
            defect_clu[d] = []  # Initialize an empty list for each defect type

        # Define threshold values for each defect type
        defect_thr = {'PLacyl': 1, 'TGglyc': 2, 'TGacyl': 3}

        # Process each defect type
        for d in defects:
            for i, ts in enumerate(defect_uni[d].trajectory):
                num = 0

                # Identify defect locations
                bA = (Mup[i] == defect_thr[d]) 
                defect_clu[d].append(bA.astype(int))
                ind = np.where(bA)
                xs, ys = ind[1], ind[0]

                # Place defects in the universe
                for x1, y1 in zip(xs, ys):
                    if num >= self.N:
                        break
                    pos = np.array([x1, y1, zlimup[i]])
                    defect_uni[d].atoms[num].position = pos
                    num += 1

                # Repeat for the lower leaflet
                bA = (Mdw[i] == defect_thr[d])
                defect_clu[d].append(bA.astype(int))
                ind = np.where(bA)
                xs, ys = ind[1], ind[0]

                for x1, y1 in zip(xs, ys):
                    if num >= self.N:
                        break
                    pos = np.array([x1, y1, zlimdw[i]])
                    defect_uni[d].atoms[num].position = pos
                    num += 1


        output_base_dir = 'GRO_paper'  # Base directory for outputs
        for d in defects:
            # Use only the prefix to set the output directory
            output_dir = os.path.join(self.prefix, d)
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)

            u = defect_uni[d]

            for i, ts in enumerate(u.trajectory):
                self.protein_atoms.universe.trajectory[i]  

                combined_universe = mda.Merge(self.protein_atoms, u.atoms)

                # Update positions in the combined universe
                combined_universe.atoms.positions[len(self.protein_atoms):] = ts.positions
                combined_universe.atoms.positions[:len(self.protein_atoms)] = self.protein_atoms.positions
                combined_universe.trajectory.ts.dimensions = ts.dimensions

                # Construct file path for writing the combined universe
                output_filepath = os.path.join(output_dir, f"{d}_frame_{i}.gro")
                combined_universe.atoms.write(output_filepath)

