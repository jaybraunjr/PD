# Packing Defect Analysis
<img src="https://github.com/user-attachments/assets/be3a6ba1-96b4-40f5-9e44-fb903c30f052" width="900">


## Overview

This repository contains Python scripts and tools for identifying, analyzing, and visualizing packing defects in molecular dynamics (MD) simulations conducted using GROMACS. The primary goal is to detect, quantify, and plot defect patterns over simulation trajectories, providing insights into molecular packing irregularities.

## Features

- **Defect Detection**: Identify packing defects using trajectory analysis.
- **Visualization**: Plot defect frequencies and distributions for in-depth insights.
- **Radius-Based Analysis**: Includes tools for computing defect radii around proteins and comparing them over time.
<p>
  <img src="https://github.com/user-attachments/assets/32cbe599-252c-4f6d-b1e1-a49f109eb614" width="450" />
  <img src="https://github.com/user-attachments/assets/4e6aae24-3823-4c3d-8f9a-446fe4643c6b" width="210" />
</p>



This movie illustrates the dynamic interaction between a membrane protein and packing defects over the course of asimulation. 

https://github.com/user-attachments/assets/5f44d331-6c67-4bdd-aad1-71afb2ab3141

## Files

- **`run_defect.py`**: Executes the primary defect analysis pipeline.
- **`run_radius.py`**: Script to run radius-based defect analysis.
- **`defects2.py`**: Contains functions for calculating radii of packing defects.
- **`radius2.py`**: Contains functions for calculating radii of packing defects.


## Requirements

- Python 3.8+
- GROMACS
- MDAnalysis
- NumPy
- Matplotlib

  ## Installation

Clone the repository and install in editable mode so the `defect-tool` command
is available:

```bash
pip install -e .
```

This installs the required Python packages listed in `setup.py`. GROMACS must
be installed separately and available in your environment.

## CLI Usage

### Run packing defect extraction

```bash
defect-tool run-defect --top path/to/system.gro \
                      --traj path/to/trajectory.xtc \
                      --output-dir results/defects
```

### Run radius analysis

```bash
defect-tool run-radius --base-directory filtered_gro/ \
                      --output-base-dir processed/ \
                      --lipid-types PLacyl TGacyl TGglyc \
                      --frame-start 0 --frame-end 10 \
                      --protein-atom-count 626 \
                      --apply-protein-cutoff \
                      --cutoff-distance 1.5
```

These commands generate GRO files and plots inside the specified output
directories.

## Algorithm Overview

The defect detection workflow loads a GROMACS topology and trajectory with
MDAnalysis. For each frame:

1. A 2D grid representing the membrane leaflets is created.
2. Atom positions are mapped to the grid using residue-specific radii.
3. Connected empty regions on the grid correspond to packing defects and are
   saved as coordinates.
4. Resulting defect and protein atoms are written to GRO files for further
   analysis and visualization.

The radius workflow filters GRO files by distance from protein atoms and uses a
depth-first search to calculate contiguous defect cluster sizes across frames.

---
Happy analyzing!

