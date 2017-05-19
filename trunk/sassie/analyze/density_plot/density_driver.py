import sys
import sassie.interface.input_filter as input_filter
import sassie.analyze.density_plot.density_plot as density_plot

class Drv():

   module = 'density_plot'

   def run_me(self):

      #### BEGIN USER EDIT
      #### BEGIN USER EDIT
      #### BEGIN USER EDIT

      runname = 'run_0'
      dcdfile = './c7.dcd'
      pdbfile = './min3.pdb'
      nsegments = '1'
      ofile = 'test'
      xlength = '350.0'
      gridsp = '6.0'
      ylength = '350.0'
      basis = 'calpha'
      zlength = '350.0'
      nregions = '5'
      lowregions = '1, 124, 143, 355, 380'
      highregions = '123, 142, 353, 379, 420'
      equalweights = '1'
      weightfile = './data/chi_square_filter/x2lowweights.txt'
      save_occupancy = 'N'

      #### END USER EDIT
      #### END USER EDIT
      #### END USER EDIT

      svariables={}

      svariables['runname']      = (runname,'string')
      svariables['nsegments']      = (nsegments,'int')
      svariables['dcdfile']  	   = (dcdfile,'string')
      svariables['pdbfile']  	   = (pdbfile,'string')
      svariables['ofile']        = (ofile,'string')
      svariables['xlength']      = (xlength,'float')
      svariables['gridsp']       = (gridsp,'float')
      svariables['ylength']      = (ylength,'float')
      svariables['basis']    	   = (basis,'string')
      svariables['zlength']      = (zlength,'float')
      svariables['nregions']     = (nregions,'int')
      svariables['lowregions']   = (lowregions,'int_array')
      svariables['highregions']  = (highregions,'int_array')
      svariables['equalweights'] = (equalweights,'int')
      svariables['weightsfile']  = (weightfile,'string')
      svariables['save_occupancy']  = (save_occupancy,'string')

      error,self.variables=input_filter.type_check_and_convert(svariables)

      if(len(error)>0):
         print 'error = ',error
         sys.exit()

      runname=self.variables['runname'][0]

      import multiprocessing
      import shutil, os
      if os.path.exists(runname+'/'+self.module):
         shutil.rmtree(runname+'/'+self.module)

      txtQueue=multiprocessing.JoinableQueue()
      segvariables = [['5', '1,124,143,355,380', '123,142,353,379,420', 'CA', 'GAG']]
      density_plot.density(self.variables,segvariables,txtQueue)

   def verify_me(self):
      module = self.module

      import os, filecmp,glob
      result_path = './'+self.variables['runname'][0]+'/'+module+'/'
      expected_path = './expected_results/'+module+'/'
      result_files = glob.glob(result_path+'*')
      result_files = [os.path.split(result_files[i])[1] for i in range(len(result_files))]
      expected_files = glob.glob(expected_path+'*')
      expected_files = [os.path.split(expected_files[i])[1] for i in range(len(expected_files))]
      print '\nvariables: ',self.variables,'\n\nresult_files:   ',result_path, '\n', result_files, '\n\nexepected_files:',expected_path, '\n',expected_files

      from util import FileCmp
      flag = True
      for ifile in expected_files:
         if ifile in result_files:
            print '\ncomparing ',expected_path+ifile, result_path+ifile,
            flag = (flag and filecmp.cmp(expected_path+ifile, result_path+ifile))
            print '\n...to be...',flag
            if flag==False:
               return False
         else:
            return False
      return flag


if __name__=='__main__':
   o=Drv()
   o.run_me()
   #o.verify_me()
