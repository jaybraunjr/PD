from pmda.parallel import ParallelAnalysisBase
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from defect import defects2
import numpy as np
from MDAnalysis import Universe
import MDAnalysis as mda
import warnings

warnings.filterwarnings("ignore")

if __name__ == '__main__':
    u = mda.Universe('../traj/6.6_2.gro')  # Initialize Universe with only the topology file
    radii_file = '../defect/types_radii.json'
    pd = defects2.PackingDefect2(radii_file)

    lipid = '../top/top_all36_lipid.rtf'
    TRIO  = '../top/TRIO.rtf'
    CHYO = '../top/CHYO.rtf'
    SAPI  = '../top/toppar_all36_lipid_inositol.str'

    radii = {'POPC': pd.read_top('POPC', lipid),
             'DOPE': pd.read_top('DOPE', lipid),
             'SAPI': pd.read_top('SAPI', SAPI),
             'TRIO': pd.read_top('TRIO', TRIO),
             'CHYO': pd.read_top('CHYO', CHYO)}

    MEMB = u.select_atoms('resname POPC DOPE SAPI TRIO CHYO')

    trajectory_files = ['../traj/rep3_skip100.xtc']  # Add more filenames as needed

    # Base output directory
    base_output_dir = '../GRO_paper'
    os.makedirs(os.path.join(base_output_dir), exist_ok=True)  

    for traj_file in trajectory_files:
        u.load_new(f'{traj_file}')

        traj_id = os.path.splitext(os.path.basename(traj_file))[0] 
        output_dir = os.path.join(base_output_dir,  traj_id)
        os.makedirs(output_dir, exist_ok=True)

        prefix = os.path.join(output_dir, '')  
        
        pdPMDA = defects2.PackingDefect2Sequential([MEMB], radii, prefix=prefix, leaflet='both')
        pdPMDA.process()
        # pdPMDA = defects2.PackingDefect2Sequential([MEMB], radii, prefix=prefix, leaflet='both')
        # pdPMDA.run(n_jobs=1)
