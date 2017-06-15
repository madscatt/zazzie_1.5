'''
Driver method to run the align module
'''

import sys
import os
import shutil
import time
import locale
#import sassie.tools.align.align as align
import sassie.interface.input_filter as input_filter
import sassie.interface.torsion_angle_md_filter as torsion_angle_md_filter
import multiprocessing

sys.path.append('../../util')
import sasconfig as sasconfig

def user_variables(self, **kwargs):

    ### BEGIN USER INPUT ###
    ### BEGIN USER INPUT ###
    ### BEGIN USER INPUT ###

    self.runname         = 'run_2'
    self.infile          = 'hiv1_gag_ma.pdb'
    self.pdbfile         = 'hiv1_gag_ma.pdb'
    self.outfile         = 'hiv1_gag_ma.dcd'
    self.nsteps          = '10'

    self.topfile         = sasconfig.__bin_path__ + '/toppar/top_all27_prot_na.inp'
    self.parmfile         = sasconfig.__bin_path__ + '/toppar/par_all27_prot_na.inp'
    self.keepout         = '1'
    self.dcdfreq         = '10'
    self.charmmexe       = sasconfig.__bin_path__ + '/charmm.exe'
    self.temperature     = '300.0'
    self.rgforce         = '0.0'
    self.rgvalue         = '0.0'

    self.dna_segnames    = ""

    self.number_flexible_segments = '1'

    self.pretamd_min_steps = '100'
    self.poll_frequency = '10'

    self.all_flexible_segnames=['MA']
    self.all_snumranges=['1']
    self.all_srlow=['114']
    self.all_srnum=['20']
    self.all_moltype=['protein']

    self.psegvariables = []
    for i in xrange(locale.atoi(self.number_flexible_segments)):
        self.psegvariables.append([self.all_flexible_segnames[i],self.all_snumranges[i],self.all_srlow[i],self.all_srnum[i],self.all_moltype[i]])

    self.path = ''

    ### END USER INPUT ###
    ### END USER INPUT ###
    ### END USER INPUT ###


def run_module(self, **kwargs):
    '''
    method to run the module and/or its input filter
    only the module input filter is run if kwargs is: test_filter=True
    method is defined outside the class so that it can be used 
    by other programs such as test_module and test_module_filter
    '''

    svariables = {}

    svariables['runname']           = (str(self.runname),'string')
    svariables['infile']            = (str(self.infile),'string')
    svariables['pdbfile']           = (str(self.pdbfile),'string')
    svariables['outfile']           = (str(self.outfile),'string')
    svariables['nsteps']            = (str(self.nsteps),'int')
    svariables['topfile']           = (str(self.topfile),'string')
    svariables['parmfile']          = (str(self.parmfile),'string')
    svariables['keepout']           = (str(self.keepout),'int')
    svariables['dcdfreq']           = (str(self.dcdfreq),'int')
    svariables['charmmexe']         = (str(self.charmmexe),'string')
    svariables['temperature']       = (str(self.temperature),'float')
    svariables['rgforce']           = (str(self.rgforce),'float')
    svariables['rgvalue']           = (str(self.rgvalue),'float')

    svariables['dna_segnames']      = (str(self.dna_segnames),'string')

    svariables['number_flexible_segments']  = (str(self.number_flexible_segments),'int')
    svariables['pretamd_min_steps'] = (str(self.pretamd_min_steps),'string')
    svariables['poll_frequency']    = (str(self.poll_frequency),'float')

    svariables['path']    = (str(self.path),'string')

    error, self.variables = input_filter.type_check_and_convert(svariables)

    if(len(error) > 0):
        print 'error = ', error
        return error

    error = torsion_angle_md_filter.check_torsion_angle_md(self.variables, self.psegvariables)

    if(len(error) > 0):
        print 'error = ', error
        return error

    runname = self.variables['runname'][0]

    if os.path.exists(os.path.join(runname, self.module)):
        shutil.rmtree(os.path.join(runname, self.module))

    txtQueue = multiprocessing.JoinableQueue()
    torsion_angle_md.tamd(self.variables,psegvariables,txtQueue)

class gui_mimic_torsion_angle_md():
    '''
    gui_mimic class contains the name of the module
    '''
    module = 'torsion_angle_md'

    def __init__(self):

        user_variables(self)
        run_module(self)


if __name__ == '__main__':

    test = False  # option to run with test variables not implemented in 1.0.
    paths = None

    run_gui = gui_mimic_torsion_angle_md()
    print "time used: ", time.time() - start
