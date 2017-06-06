import os

###     BEGIN SYSADMIN EDIT ###
###     BEGIN SYSADMIN EDIT ###
###     BEGIN SYSADMIN EDIT ###

## DEBUG LEVEL 

# __level__ = 'WARNING'
__level__ = 'DEBUG'

## INSTALLATION PATH

__installation_path__ = '/share/apps/local'
__installation_bin_path__ = '/share/apps/local/bin'

## CUDA DEFINITIONS ###

#__cuda__ = False
__cuda__ = True
__cuda_path__ = '/share/apps/local/cuda'

__cuda_lib_path__ = __cuda_path__ + '/lib64'
__cuda_include_path__ = __cuda_path__ + '/include'

__cuda_nvcc__ = __cuda_path__ + '/bin/nvcc'
__cuda_gpp__ = '/usr/bin/g++'

# total on-board memory for single GPU card (4799 MiB) Tesla K20m
__gpu_memory__ = 4799


## NAMD CONFIGURATION 

#(BASED ON NAMD 2.12 and default installation soft-link to __installation_bin_path__/namd)

# NAMD MULTICORE

__namd_run_command__ = 'namd/namd2'
__namd_run_additional_arguments__ = ''
__namd_completion_string__= 'WallClock:'


# THE FOLLOWING OPTIONS NEED A DIFFERENT NAMD INSTALLATION THAN THAT PROVIDED BY THE INSTALLATION INSTRUCTIONS

# NAMD LINUX OVER NETWORK (NOT USED IN DEFAULT INSTALLATION; NOT TESTED)

#__namd_run_command__ = 'namd/charmrun namd/namd2'
#__namd_run_additional_arguments__ = ''

# NAMD CUDA (NOT USED IN DEFAULT INSTALLATION; NOT TESTED)

#__namd_run_command__ = 'namd/namd2'
#__namd_run_additional_arguments__ = '+setcpuaffinity +devices 0,1'


###     END SYSADMIN EDIT ###
###     END SYSADMIN EDIT ###
###     END SYSADMIN EDIT ###

__bin_path__ = __installation_bin_path__


