import sys
import multiprocessing
import sassie.simulate.complex_monte_carlo.nmer_dihedral as nmer_dihedral
import sassie.interface.input_filter as input_filter
import sassie.interface.complex_monte_carlo_filter as complex_monte_carlo_filter

def runme():

    #### BEGIN USER EDIT
    #### BEGIN USER EDIT
    #### BEGIN USER EDIT

    runname='run_0'
    dcdfile='run_0.dcd'
    path=''
    pdbfile='fram601.pdb'
    trials='10'
    goback='10'
    temp='300.0'
    nsegments='2'
    segbasis='CA, CA'
    npsegments='2'
    flpsegname='ENDA, ENDB'
    seglow='95, 95'
    seghigh='110, 110'
    lowrg='20.0'
    highrg='185.0'
    zflag='0'
    zcutoff='0.0'
    cflag='0'
    confile = 'constraints.txt'
    plotflag = '0'
    directedmc = '0'
    seed = '0,123'

    allsnumranges = ["2","2"]
    allsith = ["30.0, 30.0", "30.0, 30.0"]
    allsrlow = ["123,277","123,277"]
    allsrnum = ["21,5","21,5"]
    allmoltype = ["protein","protein"]

    psegvariables = []
    print 'asnr = ',allsnumranges
    print 'len asnr = ',len(allsnumranges)

    for i in xrange(len(allsnumranges)):
        print 'i = ', i
        psegvariables.append([allsnumranges[i],allsith[i],allsrlow[i],allsrnum[i],allmoltype[i]])

    path='./'

    #### END USER EDIT
    #### END USER EDIT
    #### END USER EDIT

    svariables={}

    plotflag = '0'

    svariables={}

    svariables['runname'] = (str(runname),'string')

    svariables['dcdfile']           = (str(dcdfile),'string')
    svariables['path']              = (str(path),'string')
    svariables['pdbfile']           = (str(pdbfile),'string')
    svariables['trials']            = (str(trials),'int')
    svariables['goback']            = (str(goback),'int')
    svariables['temp']              = (str(temp),'float')

    svariables['nsegments']         = (str(nsegments),'int')
    svariables['segbasis']          = (str(segbasis),'string')
    svariables['npsegments']        = (str(npsegments),'int')
    svariables['flpsegname']        = (str(flpsegname),'string')
    svariables['seglow']            = (str(seglow),'int_array')
    svariables['seghigh']           = (str(seghigh),'int_array')

    svariables['lowrg']             = (str(lowrg),'float')
    svariables['highrg']            = (str(highrg),'float')

    svariables['zflag']             = (str(zflag),'int')
    svariables['zcutoff']           = (str(zcutoff),'float')
    svariables['cflag']             = (str(cflag),'int')
    svariables['confile']           = (str(confile),'string')
    svariables['plotflag']          = (str(plotflag),'int')
    svariables['directedmc']        = (str(directedmc),'float')
    svariables['seed']              = ('0,123', 'int_array') # set this to '1,123' if you want to set the seed or '0,123' if not

    error = []

    error,variables=input_filter.type_check_and_convert(svariables)

    if(len(error)>0):
        print 'error = ',error
        sys.exit()

    monflag = '1' ; eflag = 0

    error=complex_monte_carlo_filter.check_complex(variables,psegvariables,eflag,monflag,no_file_check="true")
	
    if(len(error)>0):
		print 'error = ',error
		sys.exit()

    txtQueue=multiprocessing.JoinableQueue()
    nmer_dihedral.dihedralgenerate(variables,psegvariables,txtQueue)

if __name__=='__main__':
	runme()

