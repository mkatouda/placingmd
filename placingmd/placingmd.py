#!/usr/bin/env python

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
* AutoDock Vina 1.2.3 binary

## Install

- Install from github
pip install git+https://github.com/mkatouda/placingmd.git

- Local install
git clone https://github.com/mkatouda/placingmd.git
cd placingmd
pip install .

"""

import sys
import os
import shutil
import argparse

import yaml
from meekovina import vina_dock
from stage3 import stage3_run
from stage3.config import ffligandsString, standardChargeMethods, chargeMethodsHelp, waterModels, boxTypes


def get_parser():
    class customHelpFormatter(argparse.ArgumentDefaultsHelpFormatter,
                              argparse.RawTextHelpFormatter):
        pass

    parser = argparse.ArgumentParser(
        #formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        #formatter_class=argparse.RawTextHelpFormatter,
        formatter_class=customHelpFormatter,
        description="python script easy to generate GROMACS input files of protein-ligand complex\n"
        "after performing Autodock Vina basic docking simulation"
    )
    parser.add_argument(
        '-i', '--inp', type=str,
        help = "yaml style input file, overwriting argument values"
    )
    parser.add_argument(
        '-l', '--ligand', type=str,
        help = "ligand (PDBQT, MOL, SDF, MOL2, PDB)"
    )
    parser.add_argument(
        '-s', '--input_smiles', type=str,
         help = 'use the specified smiles string as input '
        'instead of an input file (must be inside quotes).'
    )
    parser.add_argument(
        '-lr', '--refligand', type=str,
        help = "reference ligand (PDBQT, MOL, SDF, MOL2, PDB) to determine the box center"
    )
    parser.add_argument(
        '-rqt', '--receptor_pdbqt', type=str,
        help = "rigid part of the receptor (PDBQT)"
    )
    parser.add_argument(
        '-r', '--receptor', type=str,
        help = 'merge the created coordinates file (.gro) with an\n'
        'already existing coordinate file (.pdb or .gro), '
        'e.g. for combining \n'
        'ligand coordinates with protein coordinates. The generated topology\n'
        'will contain both the ligand and the protein. If a .gro file of the\n'
        'protein is provided and there exists a corresponding .top file that\n'
        'toplogy file will be used for the protein, otherwise a new topology\n'
        'file is generated.'
    )
    parser.add_argument(
        '-t', '--md_receptortopology',
        help = 'merge the created topology file (.top) with an\n'
        'already existing topology file.\n'
        'Must be used in combination with --receptor\n'
        'with a .gro file of the protein.'
    )
    parser.add_argument(
        '-o', '--out', type=str,
        help = 'basename of the output files (file extensions will be appended).'
    )
    parser.add_argument(
        '-v', '--verbose', action='store_true',
        help = 'verbose output.'
    )
    parser.add_argument(
        '--runvina', action='store_false',
        help = 'run AutoDock Vina docking.'
    )
    parser.add_argument(
        '--runmdinput', action='store_false',
        help = 'run MD input generator.'
    )

    parser.add_argument(
        '--vina_center_x', type=float,
        help = "X coordinate of the center"
    )
    parser.add_argument(
        '--vina_center_y', type=float,
        help = "Y coordinate of the center"
    )
    parser.add_argument(
        '--vina_center_z', type=float,
        help = "Z coordinate of the center"
    )
    parser.add_argument(
        '--vina_size_x', type=float, default=22.5,
        help = "size in the X dimension (Angstroms)"
    )
    parser.add_argument(
        '--vina_size_y', type=float, default=22.5,
        help = "size in the Y dimension (Angstroms)"
    )
    parser.add_argument(
        '--vina_size_z', type=float, default=22.5,
        help = "size in the Z dimension (Angstroms)"
    )
    parser.add_argument(
        '--vina_cpu', type=int, default=4,
        help = "the number of CPUs to use (the default is to try to"
        "detect the number of CPUs or, failing that, use 1)"
    )
    parser.add_argument(
        '--vina_scoring', type=str, default='vina',
        help = "force field name: vina(default), ad4, vinardo"
    )
    parser.add_argument(
        '--vina_seed', type=int, default=0,
        help = "explicit random seed"
    )
    parser.add_argument(
        '--vina_exhaustiveness', type=int, default=8,
        help = "exhaustiveness of the global search"
        "(roughly proportional to time): 1+"
    )
    parser.add_argument(
        '--vina_max_evals', type=int, default=0,
        help = "number of evaluations in each MC run (if zero,"
        "which is the default, the number of MC steps is"
        "based on heuristics)"
    )
    parser.add_argument(
        '--vina_num_modes', type=int, default=9,
        help = "maximum number of binding modes to generate"
    )
    parser.add_argument(
        '--vina_min_rmsd', type=int, default=1,
        help = "minimum RMSD between output poses"
    )
    parser.add_argument(
        '--vina_energy_range', type=int, default=3,
        help = "maximum energy difference between the best binding"
        "mode and the worst one displayed (kcal/mol)"
    )
    parser.add_argument(
        '--vina_spacing', type=float, default=0.375,
        help = "grid spacing (Angstrom)"
    )
    parser.add_argument(
        '--vina_verbosity', type=int, default=1,
        help = "verbosity in AutoDock Vina simulation (0=no output, 1=normal, 2=verbose)"
    )
    parser.add_argument(
        '--vina_score_only', action='store_true',
        help = "evaluate the energy of the current pose or poses without strucutre optimization"
    )
    parser.add_argument(
        '--vina_local_only', action='store_true',
        help = "evaluate the energy of the current pose or poses with local structure optimization"
    )
    parser.add_argument(
        '--vina_boxauto', action='store_true',
        help = 'enable automatic box determination algorithm'
    )
    parser.add_argument(
        '--vina_gybox_ratio', type=float, default=2.5,
        help = "scaling factor of radius of gyration to determine of docking box\n"
        "with automatic box determination algorithm"
    )
    parser.add_argument(
        '--vina_exec', type=str, default='lib',
        help = "select AutoDock Vina executer"
    )
    parser.add_argument(
        '--vina_bin_path', type=str, default='vina',
        help = "specify AutoDock Vina binary path"
    )
    parser.add_argument(
        '--md_ligand_id', type=str, default='0',
        help = 'target Ligand ID(s) to generate MD input, specified\n'
        'as a comma-separated string without spaces.'
    )
    parser.add_argument(
        '--md_ffligand', type=str, default='gaff',
        help = 'Force fields to generate parameters for, specified\n'
        'as a comma-separated string without spaces:\n'
        + ', '.join(ffligandsString)
    )
    parser.add_argument(
        '--md_ffprotein', type=str,
        help = 'force field of protein.'
    )
    parser.add_argument(
        '--md_calibration', type=str,
        help = 'modify van der Waals parameters according to specified\n'
        'calibration file.'
    )
    parser.add_argument(
        '--md_keep_ligand_name', action='store_true',
        help = 'Do not rename the ligand in the output files.\n'
        'When doing e.g. solvation or binding free energy\n'
        'it is convenient to always call the ligand the\n'
        'same thing - in this case "LIG". If this option\n'
        'is set the ligand name will not be changed to "LIG".\n'
        'If you need to assign parameters to e.g. co-factors\n'
        'it is good to keep their names to tell them apart\n'
        'from ligands.'
    )
    parser.add_argument(
        '--md_ph', type=float,
        help = 'Protonate the molecule according to this pH (float).\n'
        'This does not always give correct results. It is safer\n'
        'to provide correctly protonated input files.'
    )
    #parser.add_argument(
    #    '--md_virtualhydrogens', action='store_true',
    #    help = 'Turn hydrogens into virtual interaction sites to allow longer '
    #    'timesteps (experimental).'
    #)
    parser.add_argument(
        '--md_retain_charges', action='store_true',
        help = 'keep the mol2 charges.'
    )
    parser.add_argument(
        '--md_charge_method', type=str, default='am1bcc',
        help = 'Use the specified charge method for all force fields:\n'
        + '\n'.join(chargeMethodsHelp) + '\n'
    )
    parser.add_argument(
        '--md_charge_multiplier', type=float, default=1.0,
        help = 'multiply partial charges with this factor. Can only be used\n'
        'in combination with --charge_method.'
    )
    parser.add_argument(
        '--md_box_type', type=str, default='dodecahedron',
        help = 'Type of simulation box: ' 
        + ', '.join(boxTypes)
    )
    parser.add_argument(
        '--md_box_buffer', type=float, default=1.0,
        help = 'Buffer from the solute to the edge of the\n'
        'solvent box. Set to 0 to disable solvation (and ionisation).'
    )
    parser.add_argument(
        '--md_water_model', type=str,
        help = 'Solvent model to use in topology files. If not\n'
        'specified the solvent will not be specified in\n'
        'the topology. Suggested water models are:\n'
        + ', '.join(waterModels)
    )
    parser.add_argument(
        '--md_conc', type=str, default=0.0,
        help = 'Specify salt concentration (mol/liter).'
    )
    parser.add_argument(
        '--md_pname', type=str, default='NA',
        help = 'name of the positive counter ion in Solvent.'
    )
    parser.add_argument(
        '--md_nname', type=str, default='CL',
        help = 'name of the negative counter ion in Solvent.'
    )
    args = parser.parse_args()

    print(args)

    if args.out is None:
        if args.ligand:
            args.out = os.path.splitext(os.path.basename(args.ligand))[0]
        else:
            args.out = 'ligand_out'

    return args 

def set_config(args):
    # Read config yaml file
    if args.inp is not None and os.path.isfile(args.inp):
        with open(args.inp, 'r') as f:
            conf = yaml.safe_load(f)
    else:
        conf = {}

    # Set up default config values from program arguments
    conf_def = vars(args).copy()
    del conf_def['inp']
    [conf.setdefault(k, v) for k, v in conf_def.items()]

    if isinstance(conf['md_ligand_id'], str):
        conf['md_ligand_id'] = [int(s) for s in conf['md_ligand_id'].split(',')]
    elif isinstance(conf['md_ligand_id'], int):
        conf['md_ligand_id'] = [conf['md_ligand_id']]

    return conf

def vina_dock_main(conf):
    protein_pdbqt_path = conf['receptor_pdbqt']
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
    debug = conf['verbose']

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
       
def stage3_main(conf):

    if conf['runvina']:
        ligand_id = conf['md_ligand_id']
    else:
        ligand_id = [0]
        ligand = conf['ligand']
        output = conf['out']
    smiles = conf['input_smiles']
    ffligand = conf['md_ffligand']
    ffprotein = conf['md_ffprotein']
    calibration = conf['md_calibration']
    keep_ligand_name = conf['md_keep_ligand_name']
    ph = conf['md_ph']
    retain_charges = conf['md_retain_charges']
    charge_method = conf['md_charge_method']
    charge_method = conf['md_charge_method']
    charge_multiplier = conf['md_charge_multiplier']
    mergecoordinates = conf['receptor']
    mergetopology = conf['md_receptortopology']
    box_type = conf['md_box_type']
    box_buffer = conf['md_box_buffer']
    water = conf['md_water_model']
    conc = conf['md_conc']
    pname = conf['md_pname']
    nname = conf['md_nname']
    verbose = conf['verbose']

    for lid in ligand_id:
        if conf['runvina']:
            ligand = './{}_vinaout_{:02}.mol'.format(conf['out'], lid)
            output = '{}_vinaout_{:02}'.format(conf['out'], lid)
        print('ligand:', ligand)
        stage3_run(ligand, smiles, output, ffligand, ffprotein, calibration, 
                   keep_ligand_name, ph, retain_charges, charge_method, charge_multiplier,
                   mergecoordinates, mergetopology, box_type, box_buffer,
                   water, conc, pname, nname, verbose)

def main():
    args = get_parser()
    if args.verbose: print(args)

    conf = set_config(args)

    print('======= Input configulations =======')
    for k, v in conf.items():
        print('{}: {}'.format(k, v))
    print('====================================')

    if conf['runvina']:
        print('Run AutoDock Vina docking simulation')
        vina_dock_main(conf)

    if conf['runmdinput']:
        print('Run Gromacs MD input generator')
        stage3_main(conf)

if __name__ == '__main__':
    main()
