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
import string
import time
import subprocess
from write_apbs_input import *
import sassie.sasmol.sasmol as sasmol
import sassie.sasconfig as sasconfig

#       APBS
#
#       12/05/2004       --      initial coding                	:	jc
#       01/02/2011       --      added sasmol support 		:	jc
#       08/26/2011       --      adapted for mdx 			:	jc
#       06/16/2012       --      adapted for namd v. 2.9 		:	jc
#       09/10/2012       --      adapted for apbs			:	jc
#
# LC      1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''
      APBS is the module that contains the functions
      that are used to run a series of electrostatic calculations
	on a set of structures in a supplied pdb/dcd file.

      This module is called from APBS in the main 
      GUI through the graphical_apbs.py script.

        REFERENCE:

      Baker NA, Sept D, Joseph S, Holst MJ, McCammon JA. Electrostatics of 
	nanosystems: application to microtubules and the ribosome. 
	Proc. Natl. Acad. Sci. USA 98, 10037-10041 2001. 

	M. Holst and F. Saied, Multigrid solution of the Poisson-Boltzmann equation. 
	J. Comput. Chem. 14, 105-113, 1993.

	M. Holst and F. Saied, Numerical solution of the nonlinear Poisson-Boltzmann 
	equation: Developing more robust and efficient methods. 
	J. Comput. Chem. 16, 337-364, 1995.

	M. Holst, Adaptive numerical treatment of elliptic systems on manifolds. 
	Advances in Computational Mathematics 15, 139-191, 2001. 

	R. Bank and M. Holst, A New Paradigm for Parallel Adaptive Meshing Algorithms. 
	SIAM Review 45, 291-323, 2003.

'''


def unpack_variables(variables):

    runname = variables['runname'][0]
    infile = variables['infile'][0]
    pdbfile = variables['pdbfile'][0]
    outfile = variables['outfile'][0]
    temperature = variables['temperature'][0]
    ph = variables['ph'][0]
    ion_charge = variables['ion_charge'][0]
    ion_conc = variables['ion_conc'][0]
    ion_radius = variables['ion_radius'][0]
    manual_flag = variables['manual_flag'][0]
    manual_file = variables['manual_file'][0]

    #energyfile      = variables['energyfile'][0]
    #keepout         = variables['keepout'][0]

    return runname, infile, pdbfile, outfile, temperature, ph, ion_charge, ion_conc, ion_radius, manual_flag, manual_file


def print_failure(message, txtOutput):

    txtOutput.put("\n\n>>>> RUN FAILURE <<<<\n")
    txtOutput.put(">>>> RUN FAILURE <<<<\n")
    txtOutput.put(">>>> RUN FAILURE <<<<\n\n")
    txtOutput.put(message)

    return


def rename_his(m1):

    natoms = m1.natoms()
    resname = m1.resname()
    new_resname = []
    for i in xrange(natoms):
        this_resname = resname[i]
        if(this_resname == 'HSE' or this_resname == 'HSD' or this_resname == 'HSP'):
            new_resname.append('HIS')
        else:
            new_resname.append(this_resname)

    m1.setResname(new_resname)

    return


def apbs_driver(variables, txtOutput):
    '''
    APBS_DRIVER is the function to read in variables from GUI input and 
    used to run a series of apbs calculations
    on a set of structures in a supplied pdb/dcd file.

    INPUT:  variable descriptions:

    runname:          run_name 
    infile:           input pdb or dcd filename
    pdbfile:          input pdb file (reference)
    temperature:	  temperature of simulation

    OUTPUT:

    txtOutput:        TK handler for output to GUI textbox

    files stored in ~/run_name/apbs directory:
    outfile:          output filename 

    '''

    runname, infile, pdbfile, outfile, temperature, ph, ion_charge, ion_conc, ion_radius, manual_flag, manual_file = unpack_variables(
        variables)

    keepout = 1
    dcdfreq = 1

    path = runname + '/apbs/'
    print 'path = ', path
    print 'infile = ', infile

    vers = 'version 0.1 : 09/10/12 : jc'
    direxist = os.path.exists(path)
    if(direxist == 0):
        try:
            result = os.system('mkdir -p ' + path)
        except:
            message = 'can not create project directory: ' + path
            message += '\nstopping here\n'
            print_failure(message, txtOutput)
        if(result != 0):
            message = 'can not create project directory: ' + path
            message += '\nstopping here\n'
            print_failure(message, txtOutput)

    m1 = sasmol.SasMol(0)
    m2 = sasmol.SasMol(0)
    m1.read_pdb(pdbfile)
    m2.read_pdb(pdbfile, fastread=True)

    rename_his(m1)
    rename_his(m2)

    try:
        if(infile[-3:] == 'dcd'):
            infiletype = 'dcd'
        elif(infile[-3:] == 'pdb'):
            infiletype = 'pdb'

    except:
        message = 'input filename is a PDB or DCD file but it must end with ".pdb" or ".dcd" '
        message += ' :  stopping here'
        print_failure(message, txtOutput)

    print 'infiletype = ', infiletype

    if(infiletype == 'dcd'):
        min_max = m2.calc_minmax_all_steps(infile)
        dcdfile = m1.open_dcd_read(infile)
        nf = dcdfile[2]
    else:
        m1.read_pdb(infile)
        nf = m1.coor()[:, 0, 0].shape[0]
        min_max = m2.calc_minmax_all_steps(infile, pdb='pdb')

    print 'number of frames = ', nf

    print 'min_max = ', min_max

    maximum_dimensions = [min_max[1][0] - min_max[0][0],
                          min_max[1][1] - min_max[0][1], min_max[1][2] - min_max[0][2]]
    print 'min_max = ', min_max
    print 'maximum_dimensions = ', maximum_dimensions

    # ttxt=time.ctime()
    ttxt = time.asctime(time.gmtime(time.time()))
    st = ''.join(['=' for x in xrange(60)])

    txtOutput.put("\n%s \n" % (st))
    txtOutput.put("DATA FROM RUN: %s \n\n" % (ttxt))

    final_energy = []
    coorlist = []
    for i in range(nf):
        print 'apbs calculation for frame ', i + 1, ' of ', nf
        print 'apbs calculation for frame ', i + 1, ' of ', nf
        print 'apbs calculation for frame ', i + 1, ' of ', nf
        print 'writing temporary PDB file'

        if(infiletype == 'dcd'):
            m1.read_dcd_step(dcdfile, i)
            m1.write_pdb(path + 'junk.pdb', 0, 'w')
        else:
            m1.write_pdb(path + 'junk.pdb', i, 'w')

        print 'writing temporary APBS input file'
        if(i < 9):
            istr = '0000' + str(i + 1)
        elif(i < 99):
            istr = '000' + str(i + 1)
        elif(i < 999):
            istr = '00' + str(i + 1)
        elif(i < 9999):
            istr = '0' + str(i + 1)
        elif(i < 99999):
            istr = str(i + 1)
        else:
            print 'wow, man!'
            istr = str(i + 1)

        thisdcd = path + 'min_' + istr + '.dcd'

        if(manual_flag == 0):
            inputfilename = 'junk.in'
            write_apbs_input(maximum_dimensions, temperature,
                             inputfilename, ion_charge, ion_conc, ion_radius)
        else:
            inputfilename = manual_file

        print 'starting apbs calculation ( nfiles = ', nf, ')'
        ttime = time.ctime()
        runstring = vers + ' : ' + outfile + ' run stated at : ' + ttime
        print runstring
        ncpu = 1
        bin_path = sasconfig._bin_path
        if(ncpu == 1):
            print 'starting pdb2pqr calculation number: ', istr
            #run_pdb2pqr = 'python /usr/local/bin/pdb2pqr/pdb2pqr.py --ff=charmm --with-ph='+str(ph)+' -v '+path+'junk.pdb junk.pqr >& pdb2pqr.out'
            run_pdb2pqr = 'python ' + bin_path + 'pdb2pqr.py --ff=charmm --with-ph=' + \
                str(ph) + ' -v ' + path + 'junk.pdb junk.pqr >& pdb2pqr.out'

            os.system(run_pdb2pqr)

            print 'starting apbs calculation number: ', istr
            #nst='/usr/local/bin/apbs junk.in >& junk.out &'
            nst = bin_path + '/apbs junk.in >& junk.out &'

            p = subprocess.Popen(nst, shell=True, executable='/bin/bash')
            sts = os.waitpid(p.pid, 0)[1]
            print 'p.pid = ', p.pid
            thisjob = str(int(p.pid) + 1)

        run = 1
        esteps = 0
        while(run == 1):
            # time.sleep(5)
            lsst = 'ls junk.out | grep -c "junk.out" '
            lsfile = os.popen(lsst, 'r').readlines()
            stls = string.split(lsfile[0])
            nstls = locale.atoi(stls[0])
            if(nstls > 0):
                tout2 = os.popen(
                    'tail -15 junk.out | grep "Thanks for using"', 'r').readlines()

            if(len(tout2) > 0):
                print 'finished apbs calculation'
                run = 0

        fraction_done = (float(i + 1) / float(nf))
        progress_string = 'COMPLETED ' + \
            str(i + 1) + ' of ' + str(nf) + ' : ' + \
            str(fraction_done * 100.0) + ' % done'
        print('%s\n' % progress_string)
        print('%s\n' % progress_string)
        report_string = 'STATUS\t' + str(fraction_done)
        txtOutput.put(report_string)

        print 'finished run'

        mvst = 'mv io.mc ' + path + 'apbs_' + istr + '_io.mc'
        os.system(mvst)
        mvst = 'mv pot.dx ' + path + 'apbs_' + istr + '_pot.dx.mc'
        os.system(mvst)
        mvst = 'mv pdb2pqr.out ' + path + 'apbs_' + istr + '_pdb2pqr.dat'
        os.system(mvst)
        mvst = 'mv ' + path + 'junk.pdb ' + path + 'apbs_' + istr + '.pdb'
        os.system(mvst)
        mvst = 'mv junk.out ' + path + 'apbs_' + istr + '.out'
        os.system(mvst)
        mvst = 'mv junk.pqr ' + path + 'apbs_' + istr + '.pqr'
        os.system(mvst)
        mvst = 'mv junk.propka ' + path + 'apbs_' + istr + '.propka'
        os.system(mvst)
#		mvst = 'mv junk-input.p '+path+'apbs_input.p.'+istr+'.pqr'
#		os.system(mvst)
        mvst = 'mv junk.in ' + path + 'apbs_' + istr + '.in'
        os.system(mvst)

        #os.system('mv energy_results.out '+path+'energy_results_'+istr+'.out')

    if(infiletype == 'dcd'):
        m1.close_dcd_read(dcdfile[0])

    txtOutput.put("Total number of frames = %d\n\n" % (nf))
    txtOutput.put("output energies saved to : %s\n" % ('./' + path))
    txtOutput.put("\n%s \n" % (st))
    time.sleep(0.5)
    print 'APBS IS DONE'

    return()

if __name__ == '__main__':

    # BEGIN USER EDIT
    # BEGIN USER EDIT
    # BEGIN USER EDIT

    runname = 'run_0'
    pdbfile = 'ten_mer.pdb'
    infile = 'ten_mer_two_frames.dcd'
    outfile = 'apbs.dat'
    ph = '5.5'
    temperature = '300.0'
    ion_conc = '0.15'
    ion_charge = '1.0'
    ion_radius = '1.62'
    manual_flag = '0'
    manual_file = 'test_input_file.txt'

# END USER EDIT
# END USER EDIT
# END USER EDIT

    svariables = {}

    svariables['runname'] = (runname, 'string')
    svariables['infile'] = (infile, 'string')
    svariables['pdbfile'] = (pdbfile, 'string')
    svariables['outfile'] = (outfile, 'string')
    svariables['ph'] = (ph, 'float')
    svariables['temperature'] = (temperature, 'float')
    svariables['ion_charge'] = (ion_charge, 'float')
    svariables['ion_conc'] = (ion_conc, 'float')
    svariables['ion_radius'] = (ion_radius, 'float')
    svariables['manual_flag'] = (manual_flag,'int' )
    svariables['manual_file'] = (manual_file, 'string')

    import sassie.interface.input_filter as input_filter
    error, variables = input_filter.type_check_and_convert(svariables)

    if(len(error) > 0):
        print 'error = ', error
        sys.exit()

    runname = variables['runname'][0]

    import multiprocessing
    import shutil
    import os
    if os.path.exists(runname + '/apbs'):
        shutil.rmtree(runname + '/apbs')

    txtQueue = multiprocessing.JoinableQueue()
    apbs_driver(variables, txtQueue)
