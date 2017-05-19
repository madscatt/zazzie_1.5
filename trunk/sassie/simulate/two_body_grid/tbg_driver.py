import sys

import sassie.simulate.two_body_grid.two_body_grid as two_body_grid
import sassie.interface.input_filter as input_filter
import sassie.interface.two_body_grid_filter as two_body_grid_filter
import multiprocessing

def runme():

    # BEGIN USER EDIT
    # BEGIN USER EDIT
    # BEGIN USER EDIT

    runname = 'run_0'
    path='./'
    pdbmol1='iase1.pdb'
    pdbmol2='iase2.pdb'
    ofile='output.dcd'
    accpos='0'
    pos='-20, -20, -20'
    trans='2, 2, 2'
    dtrans='20, 20, 20'
    theta='2, 2, 2'
    dtheta='45, 45, 45'
    basis='CA'
    cutoff='3.0'
    lowrg='0.0'
    highrg='400.0'
    zflag='0'
    zcutoff='0.0'
    cflag='0'
    confile='constraints.txt'
        
    nexsegments1 = '3'
    nsegments1 = 'INT1, INT2, LED1'
    reslow1 = '271, 271, 94'
    numcont1 = '18, 18, 31'
    nexsegments2 = '3'
    nsegments2 = 'INT3, INT4, LED2'
    reslow2 = '271, 271, 94'
    numcont2 = '18, 18, 31'
    
    # END USER EDIT
    # END USER EDIT
    # END USER EDIT

    svariables={}

    svariables['runname']       = (runname,'string')
    svariables['path']      = (path,'string')
    svariables['pdbmol1']       = (pdbmol1,'string')
    svariables['pdbmol2']       = (pdbmol2,'string')
    svariables['ofile']     = (ofile,'string')
    svariables['accpos']        = (accpos,'int')
    svariables['pos']           = (pos,'float_array')
    svariables['trans']     = (trans,'int_array')
    svariables['dtrans']        = (dtrans,'float_array')
    svariables['theta']     = (theta,'int_array')
    svariables['dtheta']        = (dtheta,'float_array')
    svariables['basis']     = (basis,'string')
    svariables['cutoff']        = (cutoff,'float')
    svariables['lowrg']     = (lowrg,'float')
    svariables['highrg']        = (highrg,'float')
    svariables['zflag']     = (zflag,'int')
    svariables['zcutoff']       = (zcutoff,'float')
    svariables['cflag']     = (cflag,'int')
    svariables['confile']       = (confile,'string')
    svariables['nexsegments1']  = (nexsegments1,'int')
    svariables['nsegments1']    = (nsegments1,'string')
    svariables['reslow1']       = (reslow1,'int_array')
    svariables['numcont1']      = (numcont1,'int_array')
    svariables['nexsegments2']  = (nexsegments2,'int')
    svariables['nsegments2']    = (nsegments2,'string')
    svariables['reslow2']       = (reslow2,'int_array')
    svariables['numcont2']      = (numcont2,'int_array')

    error,variables=input_filter.type_check_and_convert(svariables)
	
    if(len(error)>0):
		print 'error = ',error
		sys.exit()

    error = two_body_grid_filter.check_input_values(variables)

    if(len(error)>0):
        print 'error = ',error
        sys.exit()

    txtQueue=multiprocessing.JoinableQueue()
    two_body_grid.two_body_grid(variables,txtQueue)

if __name__=='__main__':
	runme()

