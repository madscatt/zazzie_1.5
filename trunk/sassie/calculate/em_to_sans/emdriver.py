import sys

import sassie.calculate.em_to_sans.em_to_sans as em_to_sans
import sassie.interface.input_filter as input_filter
import sassie.interface.em_to_sans_filter as em_to_sans_filter
import multiprocessing

def runme():

    # BEGIN USER EDIT
    # BEGIN USER EDIT
    # BEGIN USER EDIT

    runname = 'run_0'
    emfiletype = '1' 
    inputpath = './'
    emdensityfile = 'p20a.mrc'
    angpix = '20.0, 20.0, 20.0'
    pdbfile = 'p20a.pdb'
    threshold = '0.1'
    sansfile = 'p20a.sub'
    npoints = '100'
    qmax = '0.3'

    # END USER EDIT
    # END USER EDIT
    # END USER EDIT

    svariables={}

    svariables['runname']       = (runname,'string')
    svariables['emfiletype']    = (emfiletype,'int')
    svariables['inputpath']     = (inputpath,'string')
    svariables['emdensityfile']     = (emdensityfile,'string')
    svariables['angpix']        = (angpix,'float_array')
    svariables['pdbfile']       = (pdbfile,'string')
    svariables['threshold']     = (threshold,'float')
    svariables['sansfile']      = (sansfile,'string')
    svariables['npoints']       = (npoints,'int')
    svariables['qmax']      = (qmax,'float')
    svariables['plotflag']      = ('0','int')

    error,variables=input_filter.type_check_and_convert(svariables)
	
    if(len(error)>0):
		print 'error = ',error
		sys.exit()

    error=em_to_sans_filter.check_em_to_sans(variables)

    if(len(error)>0):
        print 'error = ',error
        sys.exit()

    txtQueue=multiprocessing.JoinableQueue()
    em_to_sans.em(variables,txtQueue) 

if __name__=='__main__':
	runme()

