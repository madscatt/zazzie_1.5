import sys
import sassie.interface.input_filter as input_filter
import sassie.interface.sld_mol_filter as sld_mol_filter
import sassie.calculate.sld_mol.sld_mol as sld_mol

class Drv():

   module = 'sld_mol'

   def run_me(self):

      # BEGIN USER EDIT
      # BEGIN USER EDIT
      # BEGIN USER EDIT

      runname = 'run_0'
      path = './'
      pdbfile = 'nef_nohis.pdb'
      dcdfile = 'a5.dcd'
      expdatafile = 'SNS_dDPPG_myrdNef_nSLD.txt'
      outputfile = 'tmp_smoothdist.txt'
      runtype = '0'
      bulksld = '-5.1e-07'
      xon = '0'
      numdregions = '1'
      lowres = '1'
      highres = '200'
      dbin = '0.2'
      width = '2.5'
      z0 = '25.0'
      zmin = '20.0'
      zmax = '30.0'
      zevalmin = '28.0'
      zevalmax = '145.0'
      A0 = '0.2'
      Amin = '0.3'
      Amax = '0.4'
      plotflag = '0'
      fit = '1'
      sldoffset = '0'

      # END USER EDIT
      # END USER EDIT
      # END USER EDIT

      svariables = {}

      svariables['runname'] = (runname, 'string')
      svariables['path'] = (path, 'string')
      svariables['pdbfile'] = (pdbfile, 'string')
      svariables['dcdfile'] = (dcdfile, 'string')
      svariables['expdatafile'] = (expdatafile, 'string')
      svariables['outputfile'] = (outputfile, 'string')
      svariables['runtype'] = (runtype, 'int')
      
      svariables['bulk_sld'] = (bulksld, 'float')
      svariables['xon'] = (xon, 'int')
      
      svariables['num_deut_regions'] = (numdregions, 'int')
      svariables['deut_low_res'] = (lowres, 'int_array')
      svariables['deut_high_res'] = (highres, 'int_array')

      svariables['sldfit'] = (fit, 'int')
      svariables['sldoffset'] = (sldoffset, 'float')
 
      svariables['dbin'] = (dbin, 'float')
      svariables['width'] = (width, 'float')
      svariables['zfit0'] = (z0, 'float')
      svariables['zfitmin'] = (zmin, 'float')
      svariables['zfitmax'] = (zmax, 'float')
      svariables['zevalmin'] = (zevalmin, 'float')
      svariables['zevalmax'] = (zevalmax, 'float')
      svariables['A0'] = (A0, 'float')
      svariables['Amin'] = (Amin, 'float')
      svariables['Amax'] = (Amax, 'float')
      
      # 0 == NO PLOT, 1 == matplotlib, 2 == Gnuplot
      svariables['plotflag'] = ('0', 'int')

      error, self.variables = input_filter.type_check_and_convert(svariables)

      if(len(error) > 0):
         print 'error = ', error
         sys.exit()

      error=sld_mol_filter.check_sld_mol(self.variables)

      if(len(error)>0):
         print 'error = ',error
         sys.exit()

      runname=self.variables['runname'][0]

      import multiprocessing
      import shutil, os
      if os.path.exists(runname+'/'+self.module):
         shutil.rmtree(runname+'/'+self.module)

      txtQueue=multiprocessing.JoinableQueue()
      sld_mol.sld_mol_main(self.variables,txtQueue)


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
            if ifile[-4:]=='.log':
               flag = (flag and FileCmp.cmp_skip(expected_path+ifile, result_path+ifile, [4, 55]))
            else:
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
