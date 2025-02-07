# Packing Defect Analysis in GROMACS Trajectories

## Overview

This repository contains Python scripts and tools for identifying, analyzing, and visualizing packing defects in molecular dynamics (MD) simulations conducted using GROMACS. The primary goal is to detect, quantify, and plot defect patterns over simulation trajectories, providing insights into molecular packing irregularities.

## Features

- **Defect Detection**: Identify packing defects using trajectory analysis.
- **Visualization**: Plot defect frequencies and distributions for in-depth insights.
- **Radius-Based Analysis**: Includes tools for computing defect radii around proteins and comparing them over time.

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

