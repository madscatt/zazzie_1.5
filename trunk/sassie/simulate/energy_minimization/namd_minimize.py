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
import time
import numpy
import subprocess
from write_namd_input import *
from sassie.simulate.energy_minimization.prepend_namd_input import *
import sasmol.sasmol as sasmol
import sassie.util.sasconfig as sasconfig

#       NAMD_MINIMIZE
#
#       12/05/2004       --      initial coding                	:       jc
#       01/02/2011       --      added sasmol support 		:       jc
#       08/26/2011       --      adapted for mdx 		:       jc
#       06/16/2012       --      adapted for namd v. 2.9 	:       jc
#
# LC      1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''
        MINIMIZE is the module that contains the functions
        that are used to run a series of energy minimization calculations
	on a set of structures in a supplied pdb/dcd file.

        This module is called from Structure Minization in the main 
        GUI through the graphical_minimize.py script.

        REFERENCE:

        J. C. Phillips et al.
        Journal of Computational Chemistry  26  1781-1802  (2005)

'''


def unpack_variables(variables):

    runname = variables['runname'][0]
    infile = variables['infile'][0]
    pdbfile = variables['pdbfile'][0]
    outfile = variables['outfile'][0]
    nsteps = variables['nsteps'][0]
    parmfile = variables['parmfile'][0]
    psffile = variables['psffile'][0]
    ncpu = variables['ncpu'][0]
    keepout = variables['keepout'][0]
    dcdfreq = variables['dcdfreq'][0]
    infiletype = variables['infiletype'][0]

    md = variables['md'][0]
    mdsteps = variables['mdsteps'][0]
    dielect = variables['dielect'][0]
    temperature = variables['temperature'][0]
    use_external_input_file = variables['use_external_input_file'][0]
    external_input_file = variables['external_input_file'][0]

    velocity_restart_file = variables['velocity_restart_file'][0]
    extended_system_restart_file = variables['extended_system_restart_file'][0]

    return runname, infile, pdbfile, outfile, nsteps, parmfile, psffile, ncpu, infiletype, keepout, dcdfreq, md, mdsteps, dielect, temperature, use_external_input_file, external_input_file, velocity_restart_file, extended_system_restart_file


def print_failure(message, txtOutput):

    txtOutput.put("\n\n>>>> RUN FAILURE <<<<\n")
    txtOutput.put(message)
    st = ''.join(['=' for x in xrange(60)])
    txtOutput.put("\n%s \n" % (st))
    time.sleep(1.5)

    return


def minimize(variables, txtOutput):
    '''
    MINIMIZE is the function to read in variables from GUI input and 
    used to run a series of energy minimization calculations
    on a set of structures in a supplied pdb/dcd file.

    INPUT:  variable descriptions:

            runname:              	run_name 
            infile:               	input pdb or dcd filename
            pdbfile:               	input pdb file (reference)
            nsteps:               	number of steps
            parmfile:              	path nand name of topology file
                    psffile: 		name of psf file
                    ncpu			number of cpu to use
                    keepout: 		keep output files (0==no, 1==yes)
                    dcdfreq: 		save individual dcd frequency

    OUTPUT:

            txtOutput:              TK handler for output to GUI textbox

            files stored in ~/run_name/energy_minimization directory:
            outfile:                   output filename (dcd usually)

    '''

    runname, infile, pdbfile, outfile, nsteps, parmfile, psffile, ncpu, infiletype, keepout, dcdfreq, md, mdsteps, dielect, temperature, use_external_input_file, external_input_file, velocity_restart_file, extended_system_restart_file = unpack_variables(
        variables)

    path = runname + '/energy_minimization/'
    print 'path = ', path
    print 'infile = ', infile

    vers = 'version 0.7 : 06/16/12 : jc'
    direxist = os.path.exists(path)
    if(direxist == 0):
        try:
            result = os.system('mkdir -p ' + path)
            os.system('cp ' + psffile + ' ' + path)
            os.system('cp ' + pdbfile + ' ' + path)
        except:
            message = 'can not create project directory: ' + path
            message += '\nstopping here\n'
            print_failure(message, txtOutput)
        if(result != 0):
            message = 'can not create project directory: ' + path
            message += '\nstopping here\n'
            print_failure(message, txtOutput)

    m1 = sasmol.SasMol(0)
    m1.read_pdb(pdbfile)

    print 'infiletype = ', infiletype

    try:
        if(infile[-3:] == 'dcd'):
            infiletype = 'dcd'
        elif(infile[-3:] == 'pdb'):
            infiletype = 'pdb'

    except:
        message = 'input filename is a PDB or DCD file but it must end with ".pdb" or ".dcd" '
        message += ' :  stopping here'
        print_failure(message, txtOutput)

    if(infiletype == 'dcd'):
        dcdfile = m1.open_dcd_read(infile)
        nf = dcdfile[2]
    else:
        m1.read_pdb(infile)
        nf = m1.coor()[:, 0, 0].shape[0]

    print 'number of frames = ', nf

    # ttxt=time.ctime()
    ttxt = time.asctime(time.gmtime(time.time()))
    st = ''.join(['=' for x in xrange(60)])

    txtOutput.put("\n%s \n" % (st))
    txtOutput.put("DATA FROM RUN: %s \n\n" % (ttxt))

    dcdlist = []
    coorlist = []
    for i in range(nf):
        print '\nminimizing frame ', i + 1, ' of ', nf
        print 'minimizing frame ', i + 1, ' of ', nf
        print 'minimizing frame ', i + 1, ' of ', nf
        print 'writing temporary PDB file'

        if(infiletype == 'dcd'):
            m1.read_dcd_step(dcdfile, i)
            m1.write_pdb(path + 'junk.pdb', 0, 'w')
        else:
            m1.write_pdb(path + 'junk.pdb', i, 'w')

        print 'writing temporary NAMD input file'
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
        dcdlist.append(thisdcd)

        if use_external_input_file:
            write_namd_input('temp.inp', str(nsteps), str(
                dcdfreq), path + 'junk.pdb', psffile, thisdcd, parmfile, md, mdsteps, dielect, temperature)
        else:
            prepend_namd_input('temp.inp', path + 'junk.pdb', psffile, thisdcd, parmfile,
                               external_input_file, velocity_restart_file, extended_system_restart_file)

        print 'starting minimization ( nfiles = ', nf, ')'
        ttime = time.ctime()
        runstring = vers + ' : ' + outfile + ' run stated at : ' + ttime
        print runstring
        print 'starting namd minimization'
        #bin_path = sasconfig._bin_path
        bin_path = sasconfig.__bin_path__ + os.sep
        if(ncpu == 1):
            #nst='/usr/local/bin/namd/namd2 temp.inp >& junk.out &'
            nst = bin_path + 'namd/namd2 temp.inp >& junk.out &'
            print '> nst = ', nst
        else:
            print '\n> using ' + str(ncpu) + ' cpus\n'
            #nst='/usr/local/bin/namd/charmrun +p'+str(ncpu)+' /usr/local/bin/namd/namd2 temp.inp >& junk.out &'
            #if sasconfig.__cluster__:
            if sasconfig.__arch__ == "cluster" :
                nst = bin_path + 'namd/charmrun +p' + \
                    str(ncpu) + ' ++remote-shell ssh ' + bin_path + \
                    'namd/namd2 temp.inp >& junk.out &'
            else:
                nst = bin_path + 'namd/charmrun +p' + \
                    str(ncpu) + ' ' + bin_path + \
                    'namd/namd2 temp.inp >& junk.out &'

        p = subprocess.Popen(nst, shell=True, executable='/bin/bash')
        sts = os.waitpid(p.pid, 0)[1]
        print 'p.pid = ', p.pid
        thisjob = str(int(p.pid) + 1)

        run = 1
        esteps = 0
        while(run == 1):
            time.sleep(10)
            lsst = 'ls junk.out | grep -c "junk.out" '
            lsfile = os.popen(lsst, 'r').readlines()
            stls = string.split(lsfile[0])
            nstls = locale.atoi(stls[0])
            if(nstls > 0):
                if(ncpu == 1 or not sasconfig.__arch__ == "cluster"):
                    tout2 = os.popen(
                        'tail -15 junk.out | grep "Program finished"', 'r').readlines()
                else:
                    tout2 = os.popen(
                        'tail -15 junk.out | grep "WallClock:"', 'r').readlines()

                fail_check = os.popen(
                    'tail -200 junk.out | grep "Abort"', 'r').readlines()
                if len(fail_check) == 0:
                    fail_check = os.popen(
                        'tail -200 junk.out | grep "FATAL"', 'r').readlines()
                if len(fail_check) == 0:
                    fail_check = os.popen(
                        'tail -200 junk.out | grep "ERROR: Exiting prematurely"', 'r').readlines()
                if(ncpu > 1 and sasconfig.__arch__ == "cluster"):
                    tout2 = os.popen(
                        'tail -15 junk.out | grep "WallClock:"', 'r').readlines()

            if(len(tout2) > 0):
                print 'finished minimization'
                run = 0
            if(len(fail_check) > 0):
                message = '\n>>>> investigate error in junk.out file <<<<\n\n'
                message += "".join(fail_check) + '\n'
                print_failure(message, txtOutput)
                return()

        time.sleep(10)
        fraction_done = (float(i + 1) / float(nf))
        progress_string = 'COMPLETED ' + \
            str(i + 1) + ' of ' + str(nf) + ' : ' + \
            str(fraction_done * 100.0) + ' % done'
        print('%s\n' % progress_string)
        print('%s\n' % progress_string)
        report_string = 'STATUS\t' + str(fraction_done)
        txtOutput.put(report_string)

        if(keepout == 1):
            os.system('mv junk.out ' + path + 'min_' + istr + '.out')
        else:
            os.system('rm -f junk.out')
        os.system('rm -f junk.coor junk.xs* junk.vel ')
        try:
            os.system('rm -f ' + path + 'junk.pdb ')
        except:
            print '\n> could not move junk.pdb'

        os.system('rm -f ' + path + '*.BAK')

        print 'thisdcd = ', thisdcd
        temp_mol = sasmol.SasMol(0)
        temp_mol.read_pdb(pdbfile, fastread=True)

        header = temp_mol.open_dcd_read(thisdcd)
        temp_mol.close_dcd_read(header[0])

        ndcdfiles = header[2]

        # get the last configuration from the dcd file
        if ndcdfiles > 1:
            temp_mol.read_dcd(thisdcd)
            nframes = temp_mol.number_of_frames()
            natoms = temp_mol.natoms()
            coor = numpy.zeros((1, natoms, 3), numpy.float32)
            coor[0, :, :] = temp_mol.coor()[nframes - 1]
            os.system('rm -f ' + thisdcd)
            temp_mol.setCoor(coor)
            temp_mol.write_dcd(thisdcd)
            # temp_mol.write_pdb(thisdcd+'.pdb',0,'w')
        if ndcdfiles < 1:
            print 'ndcdfiles = ', ndcdfiles
            message = 'Did not save any dcd files.  Decrease dcd write frequency?'
            message + ' :  stopping here'
            print_failure(message, txtOutput)
            sys.exit()

    print '\n> finished minimizing all frames\n'

    if(infiletype == 'dcd'):
        m1.close_dcd_read(dcdfile[0])

    final_mol = sasmol.SasMol(0)
    final_molw = sasmol.SasMol(1)
    final_mol.read_pdb(pdbfile, fastread=True)
    final_molw.read_pdb(pdbfile, fastread=True)

    finaldcdfile = final_molw.open_dcd_write(path + outfile)

    for i in range(len(dcdlist)):
        print 'i = ', i
        final_mol.read_dcd(dcdlist[i])
        final_molw.setCoor(final_mol.coor())
        final_molw.write_dcd_step(finaldcdfile, 0, i + 1)

    final_molw.write_pdb(path + outfile + '.pdb', 0, 'w')

    final_molw.close_dcd_write(finaldcdfile)

    rmcmd = 'rm -f '
    for i in range(len(dcdlist)):
        rmcmd = rmcmd + dcdlist[i] + ' '
    print 'rmcmd = ', rmcmd
    os.system(rmcmd)
    os.system('mv temp.inp ' + path)

    txtOutput.put("Total number of frames = %d\n\n" % (nf))
    txtOutput.put("Minimized structures saved to : %s\n" % ('./' + path))
    txtOutput.put("\n%s \n" % (st))
    time.sleep(0.5)
    print 'NAMD MINIMIZATION IS DONE'
 
    return()

if __name__ == '__main__':

    ### BEGIN USER EDIT ###
    ### BEGIN USER EDIT ###
    ### BEGIN USER EDIT ###

    runname = 'run_0'
    infile = 'ten_mer.pdb'
    pdbfile = 'ten_mer.pdb'
    outfile = 'min_ten_mer.dcd'
    nsteps = '100'
    parmfile = '/share/apps/local/bin/toppar/par_all27_prot_na.inp'
    psffile = 'ten_mer.psf'
    ncpu = '2'
    keepout = '1'
    dcdfreq = '20'
    infiletype = 'pdb'
    md = '0'
    mdsteps = '10'
    dielect = '80'
    temperature = '300.0'

    use_external_input_file = 'True'
    external_input_file = 'external_input_2.inp'

    velocity_restart_file = 'False'
    extended_system_restart_file = 'False'

    ### END USER EDIT ###
    ### END USER EDIT ###
    ### END USER EDIT ###

    svariables = {}

    svariables['runname'] = (runname, 'string')
    svariables['infile'] = (infile, 'string')
    svariables['pdbfile'] = (pdbfile, 'string')
    svariables['outfile'] = (outfile, 'string')
    svariables['nsteps'] = (nsteps, 'int')
    svariables['parmfile'] = (parmfile, 'string')
    svariables['psffile'] = (psffile, 'string')
    svariables['ncpu'] = (ncpu, 'int')
    svariables['keepout'] = (keepout, 'int')
    svariables['dcdfreq'] = (dcdfreq, 'int')
    svariables['infiletype'] = (infiletype, 'string')
    svariables['md'] = (md, 'int')
    svariables['mdsteps'] = (mdsteps, 'int')
    svariables['dielect'] = (dielect, 'float')
    svariables['temperature'] = (temperature, 'float')

    svariables['use_external_input_file'] = (
        use_external_input_file, 'boolean')
    svariables['external_input_file'] = (external_input_file, 'string')

    svariables['velocity_restart_file'] = (velocity_restart_file, 'string')
    svariables['extended_system_restart_file'] = (
        extended_system_restart_file, 'string')

    import sassie.interface.input_filter as input_filter
    import sassie.interface.minimize_filter as minimize_filter

    error, variables = input_filter.type_check_and_convert(svariables)

    if(len(error) != 0):
        print 'error = ', error
        sys.exit()
    else:
        error = minimize_filter.check_minimize(variables)
        if(len(error) != 0):
            print 'error = ', error
            sys.exit()

    import multiprocessing

    txtQueue = multiprocessing.JoinableQueue()
    process = multiprocessing.Process(
        target=minimize, args=(variables, txtQueue))
    process.start()
