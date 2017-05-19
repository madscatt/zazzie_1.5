# zazzie_1.5

## Branch for development of HPC deployment version of "zazzie" 

The following dependencies need to be installed prior to installation

### Python and required associated libraries:

python 2.7
numpy
scipy
sasmol

### Simulation and Modeling Engines:

NAMD
APBS
PDB2PQR
ROSETTA
PYROSETTA

### Compilation and run-time libraries:

GCC/G++
CUDA 

### EDIT THE FOLLOWING FILE TO PROVIDE INFORMATION REGARDING DEPENENCIES

trunk/sassie/util/sasconfig.py


### SASSIE INSTALLATION:

To install:

python installer.py

this will install sassie these additional dependencies that are provided in the distribution and compile CUDA and other required python extensions

