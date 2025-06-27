# scripts/run_defect.py
import os
import MDAnalysis as mda
from defect import defects2

def run_defect(top_file, traj_file, output_dir):
    radii_file = 'defect/types_radii.json'
    pd = defects2.PackingDefect2(radii_file)

    lipid = 'top/top_all36_lipid.rtf'
    TRIO  = 'top/TRIO.rtf'
    CHYO = 'top/CHYO.rtf'
    SAPI  = 'top/toppar_all36_lipid_inositol.str'

    radii = {
        'POPC': pd.read_top('POPC', lipid),
        'DOPE': pd.read_top('DOPE', lipid),
        'SAPI': pd.read_top('SAPI', SAPI),
        'TRIO': pd.read_top('TRIO', TRIO),
        'CHYO': pd.read_top('CHYO', CHYO)
    }

    u = mda.Universe(top_file)
    u.load_new(traj_file)
    MEMB = u.select_atoms('resname POPC DOPE SAPI TRIO CHYO')

    os.makedirs(output_dir, exist_ok=True)
    pdPMDA = defects2.PackingDefect2Sequential([MEMB], radii, prefix=output_dir, leaflet='both')
    pdPMDA.process()
