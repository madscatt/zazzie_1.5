import os, sys, shutil, glob
import multiprocessing
import sassie.interface.input_filter as input_filter
import sassie.interface.energy_minimization_filter as energy_minimization_filter
import sassie.simulate.energy_minimization.namd_minimize as namd_minimize
import sassie.util.sasconfig as sasconfig

class Drv():

   module = 'minimization'

   def run_me(self):

      #### BEGIN USER EDIT
      #### BEGIN USER EDIT
      #### BEGIN USER EDIT

      runname = 'run_0'
      infile = 'min3.pdb'
      infile = 'ten_mer.pdb'
      pdbfile = 'min3.pdb'
      pdbfile = 'ten_mer.pdb'
      outfile = 'min_run_0.dcd'
      nsteps = '1000'
      parmfile = sasconfig.__bin_path__+'/toppar/par_all27_prot_na.inp'
      psffile = 'refgag.psf'
      psffile = 'ten_mer.psf'
      ncpu = '2'
      keepout = '1'
      dcdfreq = '10'
      infiletype = 'pdb'
      md = '0'
      mdsteps = '1000'
      dielect = '80.0'
      temperature = '300.0'

      use_external_input_file = 'False'
      external_input_file = 'dum.txt'
      velocity_restart_file =  'False'
      extended_system_restart_file = 'False'
 
      #### END USER EDIT
      #### END USER EDIT
      #### END USER EDIT

      svariables={}

      svariables['runname']           = (runname,'string')
      svariables['infile']            = (infile,'string')
      svariables['pdbfile']           = (pdbfile,'string')
      svariables['outfile']           = (outfile,'string')
      svariables['nsteps']            = (nsteps,'int')
      svariables['parmfile']          = (parmfile,'string')
      svariables['psffile']           = (psffile,'string')
      svariables['ncpu']              = (ncpu,'int')
      #svariables['energyfile']       = (energyfile,'string')
      svariables['keepout']           = (keepout,'int')
      svariables['dcdfreq']           = (dcdfreq,'int')
      svariables['infiletype']        = ('dcd','string')

      svariables['md']                = (md,'int')
      svariables['mdsteps']           = (mdsteps,'int')
      svariables['dielect']           = (dielect,'float')
      svariables['temperature']       = (temperature,'float')

      svariables['use_external_input_file'] = (use_external_input_file, 'boolean')
      svariables['external_input_file'] = (external_input_file, 'string')
      svariables['velocity_restart_file'] = (velocity_restart_file, 'string')
      svariables['extended_system_restart_file'] = (extended_system_restart_file, 'string')
    
      #svariables['topfile']          = (str(topfile),'string')

      path = ''

      svariables['path']    = (str(path),'string')

      error,variables=input_filter.type_check_and_convert(svariables)

      if(len(error) != 0):
       	print 'error = ',error
        sys.exit()
      else:
      	error=energy_minimization_filter.check_minimize(variables)

       	if(len(error) != 0):
               	print 'error = ',error
         	sys.exit()

      runname=variables['runname'][0]

      import multiprocessing
      import shutil, os
      if os.path.exists(runname+'/'+self.module):
         shutil.rmtree(runname+'/'+self.module)

      txtQueue=multiprocessing.JoinableQueue()
      namd_minimize.minimize(variables,txtQueue)

      print 'back from namd_minimize'
    
      return
   
if __name__=='__main__':
   o=Drv()
   o.run_me()

