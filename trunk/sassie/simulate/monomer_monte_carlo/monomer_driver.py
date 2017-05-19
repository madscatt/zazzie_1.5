import sys

def runme():

#### BEGIN USER EDIT
#### BEGIN USER EDIT
#### BEGIN USER EDIT

	runname='run_0'
        dcdfile='run_0.dcd'
        path=''
        pdbfile='min3.pdb'
        trials='50'
        goback='10'
        temp='300.0'
        moltype='protein'
        numranges='5'
        dtheta='30.0,30.0,30.0,30.0,30.0'
        reslow='123,278,354,378,408'
        numcont='21,5,24,11,4'
        lowres1='284'
        highres1='350'
        basis='calpha'
        cutoff='3.0'
        lowrg='0.0'
        highrg='400.0'
        zflag='0'
        zcutoff='0.0'
        nbflag='0'
        nbscale='1.0'
        psffilepath='./'
        psffilename='refgag.psf'
        parmfilepath='/usr/local/bin/toppar/'
        parmfilename='par_all27_prot_na.inp'
        
        directedmc = '0.0'
        cflag = '0'
        confile = 'dum.txt'
        plotflag = '0'
        seed = '0, 123'

#### END USER EDIT
#### END USER EDIT
#### END USER EDIT

        svariables={}

        svariables['runname']           = (runname,'string')
        svariables['dcdfile']           = (dcdfile,'string')
        svariables['path']              = (path,'string')
        svariables['pdbfile']           = (pdbfile,'string')
        svariables['trials']            = (trials,'int')
        svariables['goback']            = (goback,'int')
        svariables['temp']              = (temp,'float')
        svariables['moltype']           = (moltype,'string')
        svariables['numranges']         = (numranges,'int')
        svariables['dtheta']            = (dtheta,'float_array')
        svariables['reslow']            = (reslow,'int_array')
        svariables['numcont']           = (numcont,'int_array')
        svariables['lowres1']           = (lowres1,'int')
        svariables['highres1']          = (highres1,'int')
        svariables['basis']             = (basis,'string')
        svariables['cutoff']            = (cutoff,'float')
        svariables['lowrg']             = (lowrg,'float')
        svariables['highrg']            = (highrg,'float')
        svariables['zflag']             = (zflag,'int')
        svariables['zcutoff']           = (zcutoff,'float')
        svariables['nonbondflag']       = (nbflag,'int')
        svariables['nonbondscale']      = (nbscale,'float')
        svariables['psffilepath']       = (psffilepath,'string')
        svariables['psffilename']       = (psffilename,'string')
        svariables['parmfilepath']      = (parmfilepath,'string')
        svariables['parmfilename']      = (parmfilename,'string')

        svariables['cflag']             = (cflag,'int')
        svariables['confile']           = (confile,'string')
        svariables['plotflag']          = (plotflag,'int')
        svariables['directedmc']        = (directedmc,'float')
        svariables['seed']              = (seed,'int_array')

#	sys.path.append('/home/curtisj/lib/python/')

	import sassie.simulate.monomer_monte_carlo.dihedral_monte_carlo as dihedral
        import sassie.interface.input_filter as input_filter
        #import sassie.interface.generate_filter as generate_filter

        error,variables=input_filter.type_check_and_convert(svariables)

        #eflag=0 ; monflag=1
        #error=generate_filter.check_protein(variables,eflag,monflag)

	if(len(error)>0):
		print 'error = ',error
		sys.exit()

       	runname=variables['runname'][0]

        import multiprocessing

        txtQueue=multiprocessing.JoinableQueue()
        dihedral.dihedralgenerate(variables,txtQueue)

if __name__=='__main__':
	runme()

