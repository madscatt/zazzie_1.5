# zazzie_1.5

## Branch for development of HPC deployment version of "zazzie" 

The following dependencies need to be installed prior to running the
installation script below.

### Python and required modules/libraries:

[python 2.7](https://www.python.org/downloads/release/python-2712/)

[numpy](https://www.scipy.org/scipylib/download.html)

[scipy](https://www.scipy.org/scipylib/download.html)

[sasmol](https://github.com/madscatt/sasmol)

[pandoc](http://pandoc.org)

### Simulation and modeling engines:

[CHARMM-Lite](http://charmm.chemistry.harvard.edu)

    * required for simulate/torsion angle molecular dynamics

[NAMD 2.12 multicore-cuda](http://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=NAMD)
   
    * required for simulate/energy minimization

### Compilation and run-time libraries:

[GCC/G++ with C++11 support](https://gcc.gnu.org/viewcvs/gcc/branches/)

[CUDA 8.0](https://developer.nvidia.com/cuda-downloads) 

### EDIT THE FOLLOWING FILE TO PROVIDE INFORMATION REGARDING DEPENENCIES

trunk/sassie/util/sasconfig.py

Missing dependencies and/or configuration errors will disable module(s) listed
for each depenency above.

### SASSIE INSTALLATION:

To install:

python installer.py

this will install sassie these additional dependencies that are provided in the distribution and compile CUDA and other required python extensions.


### Step by step installation instructions for various Linux OS versions can be found here: https://github.com/madscatt/zazzie_1.5/wiki

