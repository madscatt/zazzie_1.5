# System imports
from distutils.core import *
import os, sys, platform

sys.path.append('../../../../util')
import sasconfig

cuda_dir = sasconfig.__cuda_path__
cuda_include_dir = sasconfig.__cuda_include_path__
cuda_lib_dir = sasconfig.__cuda_lib_path__

cpp_source_dir=os.path.join('..','cpp_and_cuda_source')
cuda_source_dir=os.path.join('..','cpp_and_cuda_source')

cpp_library_name = 'sascalc'
cuda_library_name = 'cudaSascalc'

# Third-party modules - we depend on numpy for everything
import numpy

from numpy.distutils.core import Extension, setup

# Obtain the numpy include directory.  This logic works across numpy versions.
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

if os.path.isdir(cuda_dir):
    cuda_driver = True
else:
    cuda_driver = False 

if os.path.isfile(os.path.join(cpp_source_dir,'lib','libsascalc.a')):
    cpp_lib = True
else:
    cpp_lib = False

if sasconfig.__cuda__:
    if os.path.isfile(os.path.join(cuda_source_dir,'lib','libcudaSascalc.a')):
        cuda_lib = True
    else:
        cuda_lib = False
else:
    cuda_lib = False
 
if not cpp_lib and not cuda_lib:
    print ("Either cpp or cuda lib needs to be pre-built")
    exit(0)
if cuda_lib and not cuda_driver:
    print ("Cuda lib found but no cuda driver detected")
    exit(0)

include_dir_names = [numpy_include]
library_dir_names = []
library_names = []
macros = []

cuda_driver = False
cuda_lib = False

if cpp_lib:
    include_dir_names.append(os.path.join(cpp_source_dir,'src')) 
    library_dir_names.append(os.path.join(cpp_source_dir,'lib'))
    library_names.append(cpp_library_name)
    macros.append(('CPP_LIB','1'))
if cuda_driver:
    include_dir_names.append(cuda_include_dir)
    include_dir_names.append(cuda_lib_dir)
    library_names.append('cuda')
    library_names.append('cudart')
    macros.append(('CUDA_DRIVER','1'))
if cuda_lib:
    include_dir_names.append(os.path.join(cuda_source_dir,'src')) 
    library_dir_names.append(os.path.join(cuda_source_dir,'lib'))
    library_names.append(cpp_library_name) #ZHL hack
    library_names.append(cuda_library_name)
    macros.append(('CUDA_LIB','1'))

if cuda_driver and cuda_lib:
    macros.append(('USE_CUDA','1'))
else:
    macros.append(('USE_CPU','1'))

# extension module
    sascalc_api = Extension(name="sassie.calculate.sascalc.sascalc_library.cpp_extension.sascalc_api",sources=['sascalc_api_extension.cpp'],
                    include_dirs = include_dir_names,
                    library_dirs = library_dir_names,
                    libraries = library_names,
                    define_macros = macros,
                   )

# setup
setup(  name        = "sascalc_api",
        description = "Module for sascalc_api",
        author      = "Joseph E. Curtis",
        ext_modules = [sascalc_api]
        )

