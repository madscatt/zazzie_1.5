'''
    SASSIE: Copyright (C) 2011 Joseph E. Curtis, Ph.D.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
import os
import sys
import locale
import time
import numpy
import subprocess
import sasmol.sasmol as sasmol
import sassie.util.sasconfig as sasconfig

from simtk.openmm.app import *
from simtk.openmm import *
import simtk.unit as unit
from sys import stdout, exit, stderr


#       OPENMM
#
#       11/08/2014       --      initial coding                	:       jc
#
# LC     1         2         3         4         5         6         7
# LC567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                     *      **
'''
        OPENMM is the module that contains the functions
        that are used to initialize and manage OpenMM simulations

        REFERENCE:

        P. Eastman, M. S. Friedrichs, J. D. Chodera, R. J. Radmer, C. M. Bruns,
        J. P. Ku, K. A. Beauchamp, T. J. Lane, L.-P. Wang, D. Shukla, T. Tye,
        M. Houston, T. Stich, C. Klein, M. R. Shirts, and V. S. Pande.
        "OpenMM 4: A Reusable, Extensible, Hardware Independent Library for
        High Performance Molecular Simulation." J. Chem. Theor. Comput.
        9(1):461-469. (2013)

'''


def initialize_system(variables):
    '''
        INITIALIZE_SYSTEM is the method that sets up a simulation
        whether it is a simple energy minimization or MD simulation.

        This method enables the usage model that the system is set up
        and defined using this method and returns the required objects
        that are used to perform the actual simulation using other methods
        in this file.  The idea is that one initializes a system and then
        later sends the initialization objects and new coordinates that
        are then processed using other methods.  This avoids redundant
        definitions and reduces overhead.

        INPUT:  variable descriptions:

        pdbfile:    input pdb file (reference)
        psffile:    name of psf file
        topfile:    path nand name of charmm topology file
        parmfile:   path nand name of charmm parameter file
        integrator: OpenMM integrator with parameters

        OUTPUT:

        simulation object (simTk)

    '''

    print '\n>>> initializing OpenMM system\n'
    pdb = PDBFile(variables.pdbfile)
    psf = CharmmPsfFile(variables.psffile)

    params = CharmmParameterSet(variables.topfile, variables.parmfile)
    system = psf.createSystem(params, nonbondedMethod=NoCutoff,
                              nonbondedCutoff=1 * unit.nanometer,
                              soluteDielectric=2.0, solventDielectric=80.0,
                              constraints=HBonds)
    integrator = variables.integrator
    simulation = Simulation(psf.topology, system, integrator)

    return simulation


def energy_minimization(simulation, mol, max_steps, energy_convergence=None):
    '''
    ENERGY_MINIMIZATION is the method that performs energy minimization
    on a pre-defined simulation openmm object using coordinates from
    the supplied SasMol object

    '''

    coor = mol.coor()[0]
    positions = unit.Quantity(coor, unit.angstroms)
    print '\n>>> starting OpenMM energy minimization\n'
    simulation.context.setPositions(positions)

    if not energy_convergence:
        simulation.minimizeEnergy(maxIterations=max_steps)
    else:
        simulation.minimizeEnergy(tolerance=energy_convergence,
                                  maxIterations=max_steps)

    positions = simulation.context.getState(getPositions=True).getPositions()

    coor = numpy.zeros((1, mol.natoms(), 3), numpy.float32)
    coor[0] = numpy.array(positions.value_in_unit(unit.angstroms))

    mol.setCoor(coor)
    print '\n>>> completed OpenMM energy minimization\n'

    return


def nvt_md(simulation, mol, number_of_steps, number_of_equilibration_steps):
    '''
    NVT_MD is the method that performs vacuum MD in canonical ensemble
    on a pre-defined simulation openmm object using coordinates from
    the supplied SasMol object

    '''

    coor = mol.coor()[0]
    positions = unit.Quantity(coor, unit.angstroms)
    print '\n>>> starting OpenMM NVT simulation\n'
    simulation.context.setPositions(positions)
    print '\n>>> equilibrating for ', number_of_equilibration_steps, ' steps'

    print_steps = int(number_of_equilibration_steps / 10.0)
    simulation.reporters.append(StateDataReporter(stdout, print_steps,
                                                  step=True,
                                                  potentialEnergy=True,
                                                  temperature=True))
    simulation.step(number_of_equilibration_steps)

    # simulation.reporters.append(DCDReporter('output.dcd', print_steps/10))

    print '\n>>> NVT simulation for ', number_of_steps, ' steps'
    simulation.step(number_of_steps)

    positions = simulation.context.getState(getPositions=True).getPositions()

    coor = numpy.zeros((1, mol.natoms(), 3), numpy.float32)
    coor[0] = numpy.array(positions.value_in_unit(unit.angstroms))

    mol.setCoor(coor)
    print '\n>>> completed OpenMM NVT simulation\n'

    return


class input_variables(object):

    def __init__(self, pdbfile=None, psffile=None, topfile=None,
                 parmfile=None):
        pass

if __name__ == '__main__':

    max_steps = 5000
    energy_convergence = None
    energy_convergence = 1.0 * unit.kilojoule / unit.mole  # the default value

    temperature = 300.0
    number_of_equilibration_steps = 1000
    number_of_nvt_steps = 1000
    # 1000000*0.002 = 2 ns total simulation time
    step_size = 0.002 * unit.picoseconds

    toppath = sasconfig.__bin_path__ + '/toppar/'
    testpath = ('/Users/curtisj/Desktop/august_sassie_development/svn_utk/'
                'sassie_1.0/trunk/testing/')
    testpath = '/home/curtisj/svn_utk/svn/sassie_1.0/trunk/testing/'
    testpath = ('/Users/curtisj/subversion_working_copies/svn_utk/sassie_2.0/'
                'trunk/testing/')
    variables = input_variables()
    variables.pdbfile = testpath + 'min3.pdb'
    variables.psffile = testpath + 'refgag.psf'
    variables.topfile = toppath + 'top_all27_prot_na.inp'
    variables.parmfile = toppath + 'par_all27_prot_na.inp'
    variables.integrator = LangevinIntegrator(temperature * unit.kelvin,
                                              1 / unit.picosecond, step_size)
    simulation = initialize_system(variables)

    mol = sasmol.SasMol(0)
    mol.read_pdb(variables.pdbfile)
    energy_minimization(simulation, mol, max_steps, energy_convergence)
    mol.write_pdb(testpath + 'new_min3.pdb', 0, 'w')
    nvt_md(simulation, mol, number_of_nvt_steps, number_of_equilibration_steps)
    mol.write_pdb(testpath + 'nvt_min3.pdb', 0, 'w')
