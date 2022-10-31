# placingmd

Protein-Ligand Attached Complex Input Generator for Molecular Dynamics 

## Licence

This package is distributed under the MIT License.

## Tutorial

You can try tutorial of placingmd to know how to install and run:  
https://colab.research.google.com/github/mkatouda/placingmd/blob/main/placingmd_tutorial_jp.ipynb

## Required software

placingmd also depends on a number of other softwares, some of which can should be installed seperately. They are:  

- Required  
1. AmberTools (https://ambermd.org/AmberTools.php )  
2. acpype (https://alanwilter.github.io/acpype/)
3. openbabel (http://openbabel.org/wiki/Main_Page )  
4. gromacs (https://www.gromacs.org/) 
5. AutoDock Vina (https://github.com/ccsb-scripps/AutoDock-Vina )
6. meeko (https://github.com/forlilab/Meeko )
7. meekovina (https://github.com/mkatouda/meekovina )
8. STaGE3 (https://github.com/mkatouda/stage3 )

- Optional  
9. MATCH (available after filling in form at http://brooks.chem.lsa.umich.edu/index.php?page=registerSoftware&subdir=articles/resources/software&link=%3Ca%20href=%22downloads/MATCH_RELEASE.tar.gz%22%3EVersion%201.000%3C/a%3E&name=MATCH ) 
10. Amsol (available after filling in form at http://t1.chem.umn.edu/license/form-user.html )  
11. Gaussian 16 or Gaussian 09 (https://gaussian.com/gaussian16/)  

## Installation (Required)

- Create conda virtual environment  
```
conda create -n py38-placingmd python=3.8  
conda activate py38-placingmd  
```

- Run the following command to install required conda packages  
```
conda install -c conda-forge pyyaml numpy ambertools acpype openbabel gromacs vina  
```

- Install meekovina and STaGE3 from github  
```
pip install git+https://github.com/mkatouda/meekovina.git
pip install git+https://github.com/mkatouda/stage3.git
```

- Install placingmd from github  
```
pip install git+https://github.com/mkatouda/placingmd.git
```

- Install placingmd from local repository  
```
git clone https://github.com/mkatouda/placingmd.git
cd placingmd
pip install .
```

## Installation (Optional)

- Resiter and download MATCH program from offical web cite and extract archive file
```
tar xzvf MATCH_RELEASE.tar.gz
```

- Set up enviroment variables for MATCH
```
export PerlChemistry=/path/to/MATCH_RELEASE/PerlChemistry
export MATCH=${HOME}/path/to/MATCH_RELEASE/MATCH
export PATH=${PATH}:${MATCH}/scripts
```

## Command usage

```
usage: placingmd [-h] [-i INP] [-l LIGAND] [-s INPUT_SMILES] [-lr REFLIGAND] 
                 [-rqt RECEPTOR_PDBQT] [-r RECEPTOR] [-t MD_RECEPTORTOPOLOGY]
                 [-o OUT] [-v] [--runvina] [--runmdinput]
                 [--vina_center_x VINA_CENTER_X]
                 [--vina_center_y VINA_CENTER_Y]
                 [--vina_center_z VINA_CENTER_Z] [--vina_size_x VINA_SIZE_X]
                 [--vina_size_y VINA_SIZE_Y] [--vina_size_z VINA_SIZE_Z]
                 [--vina_cpu VINA_CPU] [--vina_scoring VINA_SCORING]
                 [--vina_seed VINA_SEED]
                 [--vina_exhaustiveness VINA_EXHAUSTIVENESS]
                 [--vina_max_evals VINA_MAX_EVALS]
                 [--vina_num_modes VINA_NUM_MODES]
                 [--vina_min_rmsd VINA_MIN_RMSD]
                 [--vina_energy_range VINA_ENERGY_RANGE]
                 [--vina_spacing VINA_SPACING]
                 [--vina_verbosity VINA_VERBOSITY] [--vina_score_only]
                 [--vina_local_only] [--vina_boxauto]
                 [--vina_gybox_ratio VINA_GYBOX_RATIO] [--vina_exec VINA_EXEC]
                 [--vina_bin_path VINA_BIN_PATH] [--md_ligand_id MD_LIGAND_ID]
                 [--md_ffligand MD_FFLIGAND] [--md_ffprotein MD_FFPROTEIN]
                 [--md_calibration MD_CALIBRATION] [--md_keep_ligand_name]
                 [--md_ph MD_PH] [--md_retain_charges]
                 [--md_charge_method MD_CHARGE_METHOD]
                 [--md_charge_multiplier MD_CHARGE_MULTIPLIER]
                 [--md_box_type MD_BOX_TYPE] [--md_box_buffer MD_BOX_BUFFER]
                 [--md_water_model MD_WATER_MODEL] [--md_conc MD_CONC]
                 [--md_pname MD_PNAME] [--md_nname MD_NNAME]

optional arguments:
  -h, --help            show this help message and exit
  -i INP, --inp INP     yaml style input file, overwriting argument values (default: None)
  -l LIGAND, --ligand LIGAND
                        ligand (PDBQT, MOL, SDF, MOL2, PDB) (default: None)
  -s INPUT_SMILES, --input_smiles INPUT_SMILES
                        use the specified smiles string as input instead of an input file (must be inside quotes). (default: None)
  -lr REFLIGAND, --refligand REFLIGAND
                        reference ligand (PDBQT, MOL, SDF, MOL2, PDB) to determine the box center (default: None)
  -rqt RECEPTOR_PDBQT, --receptor_pdbqt RECEPTOR_PDBQT
                        rigid part of the receptor (PDBQT) (default: None)
  -r RECEPTOR, --receptor RECEPTOR
                        merge the created coordinates file (.gro) with an  
                        already existing coordinate file (.pdb or .gro), e.g. for combining
                        ligand coordinates with protein coordinates. The generated topology
                        will contain both the ligand and the protein. If a .gro file of the
                        protein is provided and there exists a corresponding .top file that
                        toplogy file will be used for the protein, otherwise a new topology
                        file is generated. (default: None)
  -t MD_RECEPTORTOPOLOGY, --md_receptortopology MD_RECEPTORTOPOLOGY
                        merge the created topology file (.top) with an
                        already existing topology file.
                        Must be used in combination with --receptor
                        with a .gro file of the protein. (default: None)
  -o OUT, --out OUT     basename of the output files (file extensions will be appended). (default: None)
  -v, --verbose         verbose output. (default: False)
  --runvina             run AutoDock Vina docking. (default: True)
  --runmdinput          run MD input generator. (default: True)
  --vina_center_x VINA_CENTER_X
                        X coordinate of the center (default: None)
  --vina_center_y VINA_CENTER_Y
                        Y coordinate of the center (default: None)
  --vina_center_z VINA_CENTER_Z
                        Z coordinate of the center (default: None)
  --vina_size_x VINA_SIZE_X
                        size in the X dimension (Angstroms) (default: 22.5)
  --vina_size_y VINA_SIZE_Y
                        size in the Y dimension (Angstroms) (default: 22.5)
  --vina_size_z VINA_SIZE_Z
                        size in the Z dimension (Angstroms) (default: 22.5)
  --vina_cpu VINA_CPU   the number of CPUs to use (the default is to try todetect the number of CPUs or,
                        failing that, use 1) (default: 4)
  --vina_scoring VINA_SCORING
                        force field name: vina(default), ad4, vinardo (default: vina)
  --vina_seed VINA_SEED
                        explicit random seed (default: 0)
  --vina_exhaustiveness VINA_EXHAUSTIVENESS
                        exhaustiveness of the global search(roughly proportional to time): 1+ (default: 8)
  --vina_max_evals VINA_MAX_EVALS
                        number of evaluations in each MC run (if zero,which is the default,
                        the number of MC steps isbased on heuristics) (default: 0)
  --vina_num_modes VINA_NUM_MODES
                        maximum number of binding modes to generate (default: 9)
  --vina_min_rmsd VINA_MIN_RMSD
                        minimum RMSD between output poses (default: 1)
  --vina_energy_range VINA_ENERGY_RANGE
                        maximum energy difference between the best bindingmode 
                        and the worst one displayed (kcal/mol) (default: 3)
  --vina_spacing VINA_SPACING
                        grid spacing (Angstrom) (default: 0.375)
  --vina_verbosity VINA_VERBOSITY
                        verbosity in AutoDock Vina simulation (0=no output, 1=normal, 2=verbose) (default: 1)
  --vina_score_only     evaluate the energy of the current pose or poses without strucutre optimization (default: False)
  --vina_local_only     evaluate the energy of the current pose or poses with local structure optimization (default: False)
  --vina_boxauto        enable automatic box determination algorithm (default: False)
 --vina_gybox_ratio VINA_GYBOX_RATIO
                        scaling factor of radius of gyration to determine of docking box
                        with automatic box determination algorithm (default: 2.5)
  --vina_exec VINA_EXEC
                        select AutoDock Vina executer (default: lib)
  --vina_bin_path VINA_BIN_PATH
                        specify AutoDock Vina binary path (default: vina)
  --md_ligand_id MD_LIGAND_ID
                        target Ligand ID(s) to generate MD input, specified
                        as a comma-separated string without spaces. (default: 0)
  --md_ffligand MD_FFLIGAND
                        Force fields to generate parameters for, specified
                        as a comma-separated string without spaces:
                        gaff, gaff2, cgenff (default: gaff)
  --md_ffprotein MD_FFPROTEIN
                        force field of protein. (default: None)
  --md_calibration MD_CALIBRATION
                        modify van der Waals parameters according to specified
                        calibration file. (default: None)
  --md_keep_ligand_name
                        Do not rename the ligand in the output files.
                        When doing e.g. solvation or binding free energy
                        it is convenient to always call the ligand the
                        same thing - in this case "LIG". If this option
                        is set the ligand name will not be changed to "LIG".
                        If you need to assign parameters to e.g. co-factors
                        it is good to keep their names to tell them apart
                        from ligands. (default: False)
  --md_ph MD_PH         Protonate the molecule according to this pH (float).
                        This does not always give correct results. It is safer
                        to provide correctly protonated input files. (default: None)
  --md_retain_charges   keep the mol2 charges. (default: False)
  --md_charge_method MD_CHARGE_METHOD
                        Use the specified charge method for all force fields:
                        am1bcc: AM1 with bond charge correction (antechamber)
                        am1bcc-pol: STaGE's own more polarized bond charge correction (antechamber)
                        mmff94: MMFF94 (Open Babel)
                        eem: electronegativity equalization method (Open Babel)
                        qeq: Assign QEq (charge equilibration) partial charges (Rappe and Goddard, 1991) (Open Babel)
                        qtpie: Assign QTPIE (charge transfer, polarization and equilibration) partial charges (Chen and Martinez, 2007) (Open Babel)
                        gaussian/hf: Hatree-Fock/6-31G(d) basis set followed by RESP (Gaussian)
                         (default: am1bcc)
  --md_charge_multiplier MD_CHARGE_MULTIPLIER
                        multiply partial charges with this factor. Can only be used
                        in combination with --charge_method. (default: 1.0)
  --md_box_type MD_BOX_TYPE
                        Type of simulation box: triclinic, cubic, dodecahedron, octahedron (default: dodecahedron)
  --md_box_buffer MD_BOX_BUFFER
                        Buffer from the solute to the edge of the
                        solvent box. Set to 0 to disable solvation (and ionisation). (default: 1.0)
  --md_water_model MD_WATER_MODEL
                        Solvent model to use in topology files. If not
                        specified the solvent will not be specified in
                        the topology. Suggested water models are:
                        opc, spce, tip4pew, spc, tip3p (default: None)
  --md_conc MD_CONC     Specify salt concentration (mol/liter). (default: 0.0)
  --md_pname MD_PNAME   name of the positive counter ion in Solvent. (default: NA)
  --md_nname MD_NNAME   name of the negative counter ion in Solvent. (default: CL)
```

## Exmaples of command line usage

### Protein-ligand binding system with AutoDock Vina docking simulation

Generates MD input files for water solvated protein-ligand complex using protein pdb file and ligand mol file as input.  
GAFF2 with AM1-BCC charge method and AMBER99SB-ILDN force fields are used for ligand and protein, respectively.  
Mininum distance of 1.0 nm from the solute to the edge of the cubic periodic box.  
TIP3P water model is used and the system is neutralized adding charged counter ions (K+) or (Cl-).  

```
placingmd -l ligand.mol -lr refligand.mol -rqt protein.pdbqt -r protein.pdb -o protein_ligand_solvated \
       --vina_size_x 20.0 --vina_size_y 20.0 --vina_size_z 20.0 --vina_num_modes 9 \
       --md_ffligand gaff2 --md_charge_method am1bcc --md_ffprotein amber99sb-ildn \
       --water_model tip3p -md_box_type cubic -md_box_buff 1.0 --md_conc 0.1 --md_pname K --md_nname CL
```

### Protein-ligand binding system without AutoDock Vina docking simulation

Generates MD input files for water solvated protein-ligand complex using protein pdb file and ligand mol file as input.  
GAFF2 with AM1-BCC charge method and AMBER99SB-ILDN force fields are used for ligand and protein, respectively.  
Mininum distance of 1.0 nm from the solute to the edge of the cubic periodic box.  
TIP3P water model is used and the system is neutralized adding charged counter ions (K+) or (Cl-).  

```
placingmd -l ligand.mol -r protein.pdb -o protein_ligand_solvated \  
       --md_ffligand gaff2 --md_charge_method am1bcc --md_ffprotein amber99sb-ildn \
       --water_model tip3p -md_box_type cubic -md_box_buff 1.0 --md_conc 0.1 --md_pname K --md_nname CL
```

## Exmaples of yaml input usage

### Protein-ligand binding system

Generates MD input files for water solvated protein-ligand complex using protein pdb file and ligand mol file as input.  
GAFF2 with AM1-BCC charge method and AMBER99SB-ILDN force fields are used for ligand and protein, respectively.  
Mininum distance of 1.0 nm from the solute to the edge of the cubic periodic box.  

Prepare input yaml file input.yml:

```
# System settings: required

ligand: './inputs/1iep_ligand.sdf'
refligand: './inputs/1iep_ligand.sdf'
receptor_pdbqt: './inputs/1iep_receptorH.pdbqt'
receptor: './inputs/1iep_receptorH.pdb'

# Job settings

out: '1iep_ligand'
verbose: True
runvina: True
runmdinput: True

# AutoDock Vina settings

vina_score_only: False
vina_center_x: 15.190
vina_center_y: 53.903
vina_center_z: 16.917
vina_size_x: 20.0
vina_size_y: 20.0
vina_size_z: 20.0
vina_cpu: 8
vina_seed: 1234
vina_exhaustiveness: 8
vina_num_modes: 9
vina_min_rmsd : 1
vina_energy_range : 3
vina_spacing : 0.375
vina_verbosity : 1
vina_exec : binary
vina_bin_path : 'vina'

# MD input settings

md_ligand_id : [0, 1, 2]
md_ffligand: 'gaff2'
md_ffpritein: 'amber99sb-ildn'
md_charge_method: 'am1bcc'
md_box_type: 'dodecahedron'
md_box_buff: 1.0
md_water_model: 'tip3p'
md_conc: 0.1
md_pname: 'NA'
md_nname: 'CL'
```

Then, run placingmd in command line:

```
placingmd -i input.yml
```

Keywards of yaml file are the same in the name of command line options.  
See above explanation of command line options.  

## Author

Michio Katouda (katouda@rist.or.jp)  
