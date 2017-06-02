'''
    SASSIE  Copyright (C) 2011 Joseph E. Curtis
    This program comes with ABSOLUTELY NO WARRANTY
    This is free software, and you are welcome to redistribute it under certain
    conditions; see http://www.gnu.org/licenses/gpl-3.0.html for details.
'''
import os
import sys
from distutils.core import setup
from distutils      import sysconfig
#from distutils.extension import Extension
from numpy.distutils.core import Extension, setup

#       SETUP
#
#       12/01/2009      --      initial coding              :       jc
#       11/22/2014      --      adapted for 2.0             :       jc
#       05/25/2017      --      adapted for 1.5 HPC         :       jc
#
#LC      1         2         3         4         5         6         7
#LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''
      setup.py is the script to install and/or update sassie.

	installation: use the installer.py script (which calls this script)

	update: using the python version that was used to install python type:

	> sudo python setup.py build
	> sudo python setup.py install

'''

# Third-party modules - we depend on numpy for everything
import numpy

# Obtain the numpy include directory.  This logic works across numpy versions.

try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

import sys

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(name='sassie_2',
	version='1.99_rev_1777',
	author='Joseph E. Curtis',
	author_email='joseph.curtis@nist.gov',
	license='GPL 3',
	url='www.smallangles.net/sassie',
	platforms='Linux, Mac OS X',
	description=("A suite of programs to generate atomistic models of biological systems, calculate scattering observables, and compare results to experimental data"),
	long_description=read('README'),	
	classifiers=["Development Status :: 2.0 Release",
		"License :: OSI Approved :: GNU Public License 3",
		"Intended Audience :: Science/Research",
		"Natural Language :: English",
		"Operating System :: Linux :: MacOS :: MacOS X",
		"Programming Language :: Python :: C :: Fortran",
		"Topic :: Scientific/Engineering :: Chemistry :: Physics"],

	package_dir={'sassie':''},

	packages=['sassie','sassie/util','sassie/analyze','sassie/build','sassie/build/build_utilities','sassie/build/pdbscan','sassie/build/pdbrx','sassie/calculate','sassie/calculate/capriqorn','sassie/interface','sassie/simulate','sassie/simulate/monte_carlo','sassie/simulate/monte_carlo/extensions','sassie/simulate/monte_carlo/extensions/overlap','sassie/simulate/monte_carlo/monte_carlo_utilities','sassie/simulate/monte_carlo/monte_carlo_utilities/isopeptide_bond_torsion','sassie/simulate/monte_carlo/monte_carlo_utilities/single_stranded_nucleic_backbone_torsion','sassie/simulate/monte_carlo/monte_carlo_utilities/protein_backbone_torsion','sassie/simulate/monte_carlo/monte_carlo_utilities/double_stranded_nucleic','sassie/simulate/monte_carlo/monte_carlo_utilities/tamc_utilities','sassie/simulate/energy','sassie/simulate/energy/extensions/','sassie/simulate/energy/extensions/non_bonding','sassie/simulate/constraints','sassie/simulate/prody','sassie/tools','sassie/tools/align', 'sassie/tools/contrast_calculator','sassie/tools/data_interpolation','sassie/tools/extract_utilities','sassie/tools/merge_utilities','sassie/calculate/sascalc/sascalc_library', 'sassie/calculate/sascalc/sascalc_library/cpp_extension'],

	ext_modules=[
	Extension('sassie.simulate.energy.non_bonding_intermolecular',['sassie/simulate/energy/extensions/non_bonding/non_bonding_intermolecular.f']),
	Extension('sassie.simulate.energy.vdw',['sassie/simulate/energy/extensions/non_bonding/vdw.f']),
    Extension('sassie.simulate.monte_carlo.pairs',['sassie/simulate/monte_carlo/extensions/pairs/pairs.f'],include_dirs=[numpy_include]),
    Extension('sassie.simulate.monte_carlo.ooverlap',['sassie/simulate/monte_carlo/extensions/ooverlap/ooverlap.c'],include_dirs=[numpy_include]),
    Extension('sassie.simulate.monte_carlo.vdw_overlap',['sassie/simulate/monte_carlo/extensions/vdw_overlap/vdw_overlap.f'],include_dirs=[numpy_include]),
	Extension('sassie.simulate.monte_carlo.dna_overlap',['sassie/simulate/monte_carlo/extensions/dna_overlap/dna_overlap.f']),
	Extension('sassie.simulate.energy.electrostatics',['sassie/simulate/energy/extensions/non_bonding/electrostatics.f'])]
	)


