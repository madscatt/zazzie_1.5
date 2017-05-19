import sassie.sasmol.sasmol as sasmol
import sassie.simulate.openmm.openmm as omm
import os,sys,locale,time,numpy,subprocess
import sassie.sasconfig as sasconfig

from simtk.openmm.app import *
from simtk.openmm import *
import  simtk.unit as unit
from sys import stdout, exit, stderr

class my_variables(object):

    def __init__(self,pdbfile=None,psffile=None,topfile=None,parmfile=None):
           pass

max_steps = 5000
energy_convergence = None
energy_convergence = 1.0*unit.kilojoule/unit.mole  # the default value

temperature = 300.0
number_of_equilibration_steps = 1000
number_of_nvt_steps = 1000
step_size = 0.002*unit.picoseconds  # 1000000*0.002 = 2 ns total simulation time

toppath = sasconfig._bin_path+'/toppar/'
testpath = '/Users/curtisj/Desktop/august_sassie_development/svn_utk/sassie_1.0/trunk/testing/'
testpath = '/home/curtisj/svn_utk/svn/sassie_1.0/trunk/testing/'

variables = my_variables()
variables.pdbfile = testpath+'min3.pdb'
variables.psffile = testpath+'refgag.psf'
variables.topfile = toppath+'top_all27_prot_na.inp'
variables.parmfile = toppath+ 'par_all27_prot_na.inp'
variables.integrator = LangevinIntegrator(temperature*unit.kelvin, 1/unit.picosecond, step_size)
simulation = omm.initialize_system(variables)   

mol = sasmol.SasMol(0)
mol.read_pdb(variables.pdbfile)
omm.energy_minimization(simulation,mol,energy_convergence,max_steps)
mol.write_pdb(testpath+'new_min3.pdb',0,'w')
omm.nvt_md(simulation,mol,number_of_nvt_steps,number_of_equilibration_steps)
mol.write_pdb(testpath+'nvt_min3.pdb',0,'w')

