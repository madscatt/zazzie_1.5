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

__cuda_include_directory__ = '/share/apps/local/cuda/include'
__cuda_nvcc__ = '/share/apps/local/cuda/bin/nvcc'
__cuda_gpp__ = '/usr/bin/g++'

# total on-board memory for single GPU card (4799 MiB) Tesla K20m
__gpu_memory__ = 4799


###     END SYSADMIN EDIT ###
###     END SYSADMIN EDIT ###
###     END SYSADMIN EDIT ###


### THE FOLLOWING NEEDS TO BE CLEANED UP
### THE FOLLOWING NEEDS TO BE CLEANED UP
### THE FOLLOWING NEEDS TO BE CLEANED UP

#__bin_path__ = os.path.join(os.path.sep,*installation_bin_path)
__bin_path__ = __installation_bin_path__


### THE FOLLOWING NEEDS TO BE RE-WORKED WITH SETUP.PY IN MIND
### THE FOLLOWING NEEDS TO BE RE-WORKED WITH SETUP.PY IN MIND
### THE FOLLOWING NEEDS TO BE RE-WORKED WITH SETUP.PY IN MIND

__core_libraries_include__ = [os.path.join(os.path.sep, 'share', 'apps', 'local', 'core_libraries', 'include')]
__core_libraries_lib__ = [os.path.join(os.path.sep, 'share', 'apps', 'local', 'core_libraries', 'lib')]





