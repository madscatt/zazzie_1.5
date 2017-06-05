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

###     END SYSADMIN EDIT ###
###     END SYSADMIN EDIT ###
###     END SYSADMIN EDIT ###

__bin_path__ = __installation_bin_path__


__arch__ = "cluster"
__arch__ = "linux"

