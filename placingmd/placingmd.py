#!/usr/bin/env python
import sys
import os
import shutil
import argparse

import yaml
from meekovina import vina_dock

"""
# placingmd

Protein-Ligand Attached Complex Input Generator for Molecular Dynamics

## Requirements

* python: 3.7 or later
* numpy
* scipy
* pandas
* rdkit
* meeko
* meekovina
* vina 1.2.3 python API
* AutoDock-Vina 1.2.3 binary

## Install

- Install from github
pip install git+https://github.com/mkatouda/placingmd.git

- Local install
git clone https://github.com/mkatouda/placingmd.git
cd placingmd
pip install .

"""

def get_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="python script easy to use Autodock Vina basic docking simulation"
        #"- Ligand input from LIGAND file\n"
        #"meekovina -l LIGAND -r RECEPTOR -o OUTPUT -cx CENTER_X -cy CENTER_Y -cz CENTER_Z\n\n"
        #"- Ligand input from file and center of box is determined by the REFLIGAND file\n"
        #"meekovina -l LIGAND -r RECEPTOR -rl REFLIGAND -o OUTPUT\n\n"
        #"- Ligand input from SMILES string\n"
        #"meekovina --input_smiles INPUT_SMILES -r RECEPTOR -o OUTPUT -cx CENTER_X -cy CENTER_Y -cz CENTER_Z\n\n"
        #"- Ligand input from SMILES string and center of box is determined by the REFLIGAND file\n"
        #"meekovina --input_smiles INPUT_SMILES -r RECEPTOR -rl REFLIGAND -o OUTPUT\n"
    )
    parser.add_argument(
        "-i", "--input", type=str,
        help="yaml style input file, overwriting argument values",
    )
    parser.add_argument(
        "-l", "--ligand", type=str,
        help="ligand (PDBQT, MOL, SDF, MOL2, PDB)"
    )
    parser.add_argument(
        "-r", "--receptor", type=str,
        help="rigid part of the receptor (PDBQT)"
    )
    parser.add_argument(
        "-o", "--out", type=str,
        help="output models (PDBQT), the default is chosen based on the ligand file name"
    )
    parser.add_argument(
        "--input_smiles", type=str,
        help="SMILES string (Need to put the atom you want to extend at the end of the string)"
    )
    parser.add_argument(
        "-rl","--refligand", type=str,
        help="ligand (PDBQT, MOL, SDF, MOL2, PDB) to determine the box center"
    )
    parser.add_argument(
        "--vina_center_x", type=float,
        help="X coordinate of the center"
    )
    parser.add_argument(
        "--vina_center_y", type=float,
        help="Y coordinate of the center"
    )
    parser.add_argument(
        "--vina_center_z", type=float,
        help="Z coordinate of the center"
    )
    parser.add_argument(
        "--vina_size_x", type=float, default=22.5,
        help="size in the X dimension (Angstroms)"
    )
    parser.add_argument(
        "--vina_size_y", type=float, default=22.5,
        help="size in the Y dimension (Angstroms)"
    )
    parser.add_argument(
        "--vina_size_z", type=float, default=22.5,
        help="size in the Z dimension (Angstroms)"
    )
    parser.add_argument(
        "--vina_cpu", type=int, default=4,
        help="the number of CPUs to use (the default is to try to"
        "detect the number of CPUs or, failing that, use 1)"
    )
    parser.add_argument(
        "--vina_scoring", type=str, default='vina',
        help="force field name: vina(default), ad4, vinardo"
    )
    parser.add_argument(
        "--vina_seed", type=int, default=0,
        help="explicit random seed"
    )
    parser.add_argument(
        "--vina_exhaustiveness", type=int, default=8,
        help="exhaustiveness of the global search"
        "(roughly proportional to time): 1+"
    )
    parser.add_argument(
        "--vina_max_evals", type=int, default=0,
        help="number of evaluations in each MC run (if zero,"
        "which is the default, the number of MC steps is"
        "based on heuristics)"
    )
    parser.add_argument(
        "--vina_num_modes", type=int, default=9,
        help="maximum number of binding modes to generate"
    )
    parser.add_argument(
        "--vina_min_rmsd", type=int, default=1,
        help="minimum RMSD between output poses"
    )
    parser.add_argument(
        "--vina_energy_range", type=int, default=3,
        help="maximum energy difference between the best binding"
        "mode and the worst one displayed (kcal/mol)"
    )
    parser.add_argument(
        "--vina_spacing", type=float, default=0.375,
        help="grid spacing (Angstrom)"
    )
    parser.add_argument(
        "--vina_verbosity", type=int, default=1,
        help="verbosity (0=no output, 1=normal, 2=verbose)"
    )
    parser.add_argument(
        "--vina_score_only", action='store_true',
        help="evaluate the energy of the current pose or poses without strucutre optimization"
    )
    parser.add_argument(
        "--vina_local_only", action='store_true',
        help="evaluate the energy of the current pose or poses with local structure optimization"
    )
    parser.add_argument(
        "--vina_exec", type=str, default='lib',
        help="select AutoDock-Vina executer"
    )
    parser.add_argument(
        "--vina_bin_path", type=str, default='vina',
        help="AutoDock-Vina binary path"
    )
    parser.add_argument(
        "--vina_boxauto", action='store_true',
        help="boxauto"
    )
    parser.add_argument(
        "--vina_gybox_ratio", type=float, default=2.5,
        help="box ratio"
    )
    parser.add_argument(
        "-d", "--debug", action='store_true',
        help="debug mode"
    )
    args = parser.parse_args()

    print(args)

    if args.out is None:
        if args.ligand:
            args.out = os.path.splitext(os.path.basename(args.ligand))[0] + '_out.pdbqt'
        else:
            args.out = 'ligand_out.pdbqt'

    return args 

def set_config(args):
    # Read config yaml file
    if args.input is not None and os.path.isfile(args.input):
        with open(args.input, 'r') as f:
            conf = yaml.safe_load(f)
    else:
        conf = {}

    # Set up default config values from program arguments
    conf_def = vars(args).copy()
    del conf_def['input'], conf_def['debug']
    conf_def = {k: v for k, v in conf_def.items() if v is not None}
    [conf.setdefault(k, v) for k, v in conf_def.items()]

    return conf

def vina_dock_run(conf, debug=False):
    protein_pdbqt_path = conf['receptor']
    ligand_path = conf['ligand']
    ref_ligand_path = conf['refligand']
    vina_exec = conf['vina_exec']
    vina_bin_path = conf['vina_bin_path']
    vina_boxauto = conf['vina_boxauto']
    vina_gybox_ratio = conf['vina_gybox_ratio']
    vina_boxcenter = [conf['vina_center_x'], conf['vina_center_y'], conf['vina_center_z']]
    vina_boxsize = [conf['vina_size_x'], conf['vina_size_y'], conf['vina_size_z']]
    vina_scoring = conf['vina_scoring']
    vina_cpu = conf['vina_cpu']
    vina_seed = conf['vina_seed']
    vina_exhaustiveness = conf['vina_exhaustiveness']
    vina_max_evals = conf['vina_max_evals']
    vina_num_modes = conf['vina_num_modes']
    vina_min_rmsd = conf['vina_min_rmsd']
    vina_energy_range = conf['vina_energy_range']
    vina_spacing = conf['vina_spacing']
    vina_verbosity = conf['vina_verbosity']
    vina_score_only = conf['vina_score_only']
    vina_local_only = conf['vina_local_only']

    print('vina_boxcenter:', vina_boxcenter)

    scores = vina_dock(protein_pdbqt_path,
                       ligand_path,
                       ref_ligand_path,
                       vina_exec,
                       vina_bin_path,
                       vina_boxcenter,
                       boxauto=vina_boxauto,
                       boxsize=vina_boxsize,
                       gybox_ratio=vina_gybox_ratio,
                       scoring=vina_scoring,
                       cpu=vina_cpu,
                       seed=vina_seed,
                       exhaustiveness=vina_exhaustiveness,
                       max_evals=vina_max_evals,
                       num_modes=vina_num_modes,
                       min_rmsd=vina_min_rmsd,
                       energy_range=vina_energy_range,
                       spacing=vina_spacing,
                       verbosity=vina_verbosity,
                       score_only=vina_score_only,
                       local_only=vina_local_only,
                       debug=debug)

    return scores
       
def main():
    args = get_parser()
    debug = args.debug
    if debug: print(args)

    conf = set_config(args)

    print('======= Input configulations =======')
    for k, v in conf.items():
        print('{}: {}'.format(k, v))
    print('====================================')

    vina_dock_run(conf, debug=debug)

if __name__ == '__main__':
    main()
