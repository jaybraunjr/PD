import click
from scripts.run_defect import run_defect
from scripts.run_radius import run_radius

@click.group()
def cli():
    pass

@cli.command()
@click.option('--top', required=True, help="Topology file (.gro)")
@click.option('--traj', required=True, help="Trajectory file (.xtc/.dcd)")
@click.option('--output-dir', default='output', help="Output directory")
def run_defect_cmd(top, traj, output_dir):
    run_defect(top, traj, output_dir)

@cli.command()
@click.option('--base-directory', required=True, help="Input GRO directory")
@click.option('--output-base-dir', default='filtered_gro_paper', help="Filtered GRO output dir")
@click.option('--lipid-types', multiple=True, default=['PLacyl', 'TGacyl', 'TGglyc'])
@click.option('--frame-start', default=0)
@click.option('--frame-end', default=10)
@click.option('--protein-atom-count', default=626)
@click.option('--apply-protein-cutoff', is_flag=True)
@click.option('--cutoff-distance', default=1.5)
def run_radius_cmd(base_directory, output_base_dir, lipid_types, frame_start, frame_end,
                   protein_atom_count, apply_protein_cutoff, cutoff_distance):
    run_radius(base_directory, output_base_dir, list(lipid_types), frame_start, frame_end,
               protein_atom_count, apply_protein_cutoff, cutoff_distance)

if __name__ == '__main__':
    cli()
