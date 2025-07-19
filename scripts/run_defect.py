# scripts/run_defect.py
import os
import json
import MDAnalysis as mda
from defect import defects2

def run_defect(top_file, traj_file, output_dir,
               radii_file='defect/types_radii.json',
               defect_config=None):
    """Run packing defect extraction.

    Parameters
    ----------
    top_file : str
        Path to the system topology (.gro).
    traj_file : str
        Path to the trajectory (.xtc/.dcd).
    output_dir : str
        Directory where output GRO files will be written.
    radii_file : str, optional
        JSON file mapping atom types to radii. Defaults to
        ``defect/types_radii.json``.
    defect_config : str, optional
        JSON file mapping defect names to integer thresholds. If
        provided, these values will be forwarded to
        :class:`PackingDefect2Sequential`.
    """

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
    defect_types = None
    defect_thresholds = None
    if defect_config:
        with open(defect_config) as f:
            defect_thresholds = json.load(f)
        defect_types = list(defect_thresholds.keys())

    pdPMDA = defects2.PackingDefect2Sequential(
        [MEMB], radii,
        prefix=output_dir,
        leaflet='both',
        defect_types=defect_types,
        defect_thresholds=defect_thresholds
    )
    pdPMDA.process()
