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
from write_tamd_input import *
import sassie.sasmol.sasmol as sasmol

#       TAMD
#
#       12/05/2004       --      initial coding            	:       jc
#       01/02/2011       --      added sasmol support 	:       jc
#       08/26/2011       --      adapted for mdx 		:       jc
#       06/16/2012       --      adapted for namd v. 2.9 	:       jc
#       12/08/2013       --      adapted for tamd 		:       jc
#
# LC      1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''
	TAMD is the module that contains the functions
    	that are used to run a series of tamd dynamics calculations
	on a set of structures in a supplied pdb/dcd file.

	This module is called from TAMD in the main 
	GUI through the graphical_tamd.py script.

	REFERENCE:


'''


def unpack_variables(variables):

    runname = variables['runname'][0]
    infile = variables['infile'][0]
    pdbfile = variables['pdbfile'][0]
    outfile = variables['outfile'][0]
    nsteps = variables['nsteps'][0]

    topfile = variables['topfile'][0]
    parmfile = variables['parmfile'][0]
    keepout = variables['keepout'][0]
    dcdfreq = variables['dcdfreq'][0]
    charmmexe = variables['charmmexe'][0]
    temperature = variables['temperature'][0]
    rgforce = variables['rgforce'][0]
    rgvalue = variables['rgvalue'][0]

    dna_segnames = variables['dna_segnames'][0]

    number_flexible_segments = variables['number_flexible_segments'][0]

    pretamd_min_steps = variables['pretamd_min_steps'][0]
    poll_frequency = variables['poll_frequency'][0]

    return runname, infile, pdbfile, outfile, nsteps, topfile, parmfile, keepout, dcdfreq, charmmexe, temperature, rgforce, rgvalue, dna_segnames, number_flexible_segments, pretamd_min_steps, poll_frequency


def print_failure(message, txtOutput):

    txtOutput.put("\n\n>>>> RUN FAILURE <<<<\n")
    txtOutput.put(">>>> RUN FAILURE <<<<\n")
    txtOutput.put(">>>> RUN FAILURE <<<<\n\n")
    txtOutput.put(message)

    return


def process_input_variables(psegvariables, segname_molecules, filter_flag):

    all_flexible_segnames = []
    all_snumranges = []
    all_srlow = []
    all_srnum = []
    all_moltype = []
    all_segnames = []
    all_resid = []

    if not filter_flag:

        for mol in segname_molecules:
            all_segnames.append(mol.segname()[0])
            all_resid.append(mol.resid()[0])

    offset = 0

    for i in range(len(psegvariables)):
        if not filter_flag:

            this_flexible_segname = psegvariables[i][0]
            for j in xrange(len(all_segnames)):
                if this_flexible_segname == all_segnames[j]:
                    offset = - all_resid[j] + 1

        all_flexible_segnames.append(psegvariables[i][0])
        all_snumranges.append(psegvariables[i][1])
        if filter_flag:
            all_srlow.append(psegvariables[i][2])
        else:
            this_srlow = psegvariables[i][2]
            tmp_list = [int(val.strip()) + offset for val in this_srlow.split(',')]
            new_list = ', '.join(str(x) for x in tmp_list) 
            all_srlow.append(new_list)
        all_srnum.append(psegvariables[i][3])
        all_moltype.append(psegvariables[i][4])

    all_numranges = []

    for i in range(len(all_snumranges)):
        nr = locale.atoi(all_snumranges[i])
        all_numranges.append(nr)

    all_rlow = []
    all_rnum = []

    for i in range(len(all_srlow)):
        linrlow = string.split(all_srlow[i], ',')
        linrnum = string.split(all_srnum[i], ',')
        rlow = []
        rnum = []
        for k in range(len(linrlow)):
            trlow = locale.atoi(linrlow[k])
            trnum = locale.atoi(linrnum[k])
            rlow.append(trlow)
            rnum.append(trnum)
        all_rlow.append(rlow)
        all_rnum.append(rnum)

    print 'all_flexible_segnames = ', all_flexible_segnames
    print 'all_numranges = ', all_numranges
    print 'all_rlow = ', all_rlow
    print 'all_rnum = ', all_rnum
    print 'all_moltype = ', all_moltype

    flexible_segment_variables = {}

    for i in xrange(len(all_flexible_segnames)):

        flexible_segment_variables[all_flexible_segnames[i]] = [
            all_numranges[i], all_rlow[i], all_rnum[i], all_moltype[i]]

#	flexible_segment_variables['number_of_ranges'] = all_numranges
#	flexible_segment_variables['residue_lows'] = all_rlow
#	flexible_segment_variables['residue_numcont'] = all_rnum
#	flexible_segment_variables['moltype'] = all_moltype

#>>> fsn = {}
#>>> s = 'seg1'
#>>> nr = 2
#>>> lr = [6,99]
#>>> fsn[s] = [nr,lr]

    return flexible_segment_variables


def get_segname_molecules(m1):

    frame = 0

    m1.initialize_children()

    segnames = m1.segnames()
    segname_masks = m1.segnames_mask()

    number_of_segnames = len(segnames)

    segname_molecules = []

    for i in xrange(number_of_segnames):
        this_mol = sasmol.SasMol(0)
        error = m1.copy_molecule_using_mask(this_mol, segname_masks[i], frame)
        segname_molecules.append(this_mol)

    print 'segnames = ', segnames

    return segname_molecules, segname_masks


def renumber_residues(mol):

      # renumber resid to start at 1

    resid = mol.resid()

    number = []
    resid_array = []
    count = 1
    for i in xrange(len(resid)):
        this_resid = resid[i]
        if(i == 0):
            last_resid = this_resid
        else:
            if(this_resid != last_resid):
                count += 1
        resid_array.append(count)
        last_resid = this_resid
        number.append(i + 1)

    mol.setResid(resid_array)
    mol.setIndex(number)

    return


def update_segname_molecules(m1, segname_molecules, segname_masks, frame, path):

    #	error,coor=self.get_coor_using_mask(frame,mask)

    temp_pdb_files = []

    for i in xrange(len(segname_molecules)):
        error, coor = m1.get_coor_using_mask(frame, segname_masks[i])
        segname_molecules[i].setCoor(coor)

        renumber_residues(segname_molecules[i])

        temp_pdb_name = path + 'temp_' + str(i) + '.pdb'
        segname_molecules[i].write_pdb(temp_pdb_name, 0, 'w')
        temp_pdb_files.append(temp_pdb_name)

    return temp_pdb_files


def fix_moltypes(mol, dna_segnames):

    dna_segnames = dna_segnames.split(',')

    # if "None" not in dna_segnames:
    if len(dna_segnames) > 0:
        print 'len(dna_segnames) = ', len(dna_segnames)
        print 'dna_segnames = ', dna_segnames
        for i in xrange(len(mol)):
            if(mol[i].segname()[0] in dna_segnames):
                moltype = ['dna'] * mol[i].natoms()
                mol[i].setMoltype(moltype)
    else:
        print '>> no DNA containing segments in system'

    return


def tamd(variables, psegvariables, txtOutput):
    '''
    TAMD is the method to read in variables from GUI input and 
    used to run a series of tamd calculations
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

            files stored in ~/run_name/torsion_angle_md directory:
            outfile:                   output filename (dcd usually)

    '''

    #flexible_segment_variables = process_input_variables(psegvariables)

    runname, infile, pdbfile, outfile, nsteps, topfile, parmfile, keepout, dcdfreq, charmmexe, temperature, rgforce, rgvalue, dna_segnames, number_flexible_segments, pretamd_min_steps, poll_frequency = unpack_variables(
        variables)

    path = runname + '/torsion_angle_md/'
    print 'path = ', path
    print 'infile = ', infile
    outfile = path + outfile

    vers = 'version 0.1 : 12/08/13 : jc'
    direxist = os.path.exists(path)
    if(direxist == 0):
        try:
            result = os.system('mkdir -p ' + path)
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

    residue = m1.resid()
    number_of_residues = residue[-1] - residue[0] + 1

    print 'number_of_residues = ', number_of_residues

    print 'infile = ', infile

    infiletype = os.path.splitext(infile)[1][1:]

    print 'infiletype = ', infiletype

    if(infiletype == 'dcd'):
        dcdfile = m1.open_dcd_read(infile)
        nf = dcdfile[2]
    elif(infiletype == 'pdb'):
        m1.read_pdb(infile)
        nf = m1.coor()[:, 0, 0].shape[0]
    else:
        st = 'only .dcd or .pdb suffix allowed for infile'
        print_failure(st, txtOutput)
        sys.exit()

    print 'number of frames = ', nf

    segname_molecules, segname_masks = get_segname_molecules(m1)

    flexible_segment_variables = process_input_variables(
        psegvariables, segname_molecules, False)

    fix_moltypes(segname_molecules, dna_segnames)

    #ttxt = time.ctime()
    ttxt = time.asctime( time.gmtime( time.time() ) ) 
    st = ''.join(['=' for x in xrange(60)])

    txtOutput.put("\n%s \n" % (st))
    txtOutput.put("DATA FROM RUN: %s \n\n" % (ttxt))

    dcdlist = []
    coorlist = []
    for i in range(nf):
        print '\ntamd dynamics on frame ', i + 1, ' of ', nf
        print 'tamd dynamics on frame ', i + 1, ' of ', nf
        print 'tamd dynamics on frame ', i + 1, ' of ', nf
        print 'writing temporary PDB file'

        if(infiletype == 'dcd'):
            m1.read_dcd_step(dcdfile, i)
            temp_pdb_files = update_segname_molecules(
                m1, segname_molecules, segname_masks, 0, path)
        else:
            temp_pdb_files = update_segname_molecules(
                m1, segname_molecules, segname_masks, i, path)

        print 'writing temporary TAMD input file'
        istr = str(i + 1).zfill(5)  # 99999 maximum number of frames

        thisdcd = path + 'tamd_dyn_' + istr + '.dcd'
        dcdlist.append(thisdcd)
        if (rgvalue == "-1"):
            my_rg = str(m1.calcrg(0))
        else:
            my_rg = rgvalue

        write_tamd_input(m1, segname_molecules, 'temp.inp', nsteps, dcdfreq, temp_pdb_files, topfile, parmfile,
                         thisdcd, temperature, rgforce, my_rg, residue, flexible_segment_variables, pretamd_min_steps)

        print 'starting tamd dynamics ( nfiles = ', nf, ')'
        ttime = time.ctime()
        runstring = vers + ' : ' + outfile + ' run stated at : ' + ttime
        print runstring
        print 'starting tamd dynamics'

        nst = charmmexe + ' < temp.inp >& junk.out &'
        print '> nst = ', nst

        p = subprocess.Popen(nst, shell=True, executable='/bin/bash')
        sts = os.waitpid(p.pid, 0)[1]
        print 'p.pid = ', p.pid
        thisjob = str(int(p.pid) + 1)

        run = 1
        esteps = 0
        run_fail = False
        fail_check = None

        while(run == 1):
            time.sleep(poll_frequency)
            lsst = 'ls junk.out | grep -c "junk.out" '
            lsfile = os.popen(lsst, 'r').readlines()
            stls = string.split(lsfile[0])
            nstls = locale.atoi(stls[0])
            if(nstls > 0):
                tout2 = os.popen(
                    'tail -15 junk.out | grep "NORMAL TERMINATION"', 'r').readlines()

                fail_check = os.popen(
                    'tail -15 junk.out | grep "ABNORMAL TERMINATION"', 'r').readlines()

            if(len(tout2) > 0):
                print 'finished tamd run'
                print 'tout2 = ', tout2
                run = 0

            if(len(fail_check) > 0):
                print 'CHARMM execution has failed: check output file junk.out for details'
                print 'CHARMM execution has failed: check output file junk.out for details'
                print 'CHARMM execution has failed: check output file junk.out for details'
                run_fail = True
                run = 0

        if not run_fail:
            fraction_done = (float(i + 1) / float(nf))
            progress_string = 'COMPLETED ' + \
                str(i + 1) + ' of ' + str(nf) + ' : ' + \
                str(fraction_done * 100.0) + ' % done'
            print('%s\n' % progress_string)
            print('%s\n' % progress_string)
            report_string = 'STATUS\t' + str(fraction_done)
            print 'report_string = ', report_string
            sys.stdout.flush()

            txtOutput.put(report_string)

        else:
            os.system('mv junk.out ' + path + 'min_' + istr + '.out')
            report_string = 'STATUS\t' + str(100.0)
            txtOutput.put(report_string)
            txtOutput.put(
                "RUN FAILED : read contents of min_%s.out file in : %s\n" % (istr, path))
            txtOutput.put(
                "RUN FAILED : read contents of min_%s.out file in : %s\n" % (istr, path))
            txtOutput.put(
                "RUN FAILED : read contents of min_%s.out file in : %s\n" % (istr, path))
            txtOutput.put(
                "Files from the partial run saved to : %s\n" % ('./' + path))
            txtOutput.put("\n%s \n" % (st))
            time.sleep(1)
            break

        if(keepout == 1):
            os.system('mv junk.out ' + path + 'min_' + istr + '.out')
        elif(i == 0):
            os.system('mv junk.out ' + path + 'min_' + istr + '.out')
        else:
            os.system('rm -f junk.out')

        print 'thisdcd = ', thisdcd
        temp_mol = sasmol.SasMol(0)
        temp_mol.read_pdb('tamd_output.pdb', fastread=True)

        header = temp_mol.open_dcd_read(thisdcd)
        temp_mol.close_dcd_read(header[0])

        ndcdfiles = header[2]

        # get the last configuration from the dcd file
        if ndcdfiles > 1:
            if(keepout != 1):
                print 'keeping only the final frame from this trajectory'
                temp_mol.read_dcd(thisdcd)
                nframes = temp_mol.number_of_frames()
                natoms = temp_mol.natoms()
                coor = numpy.zeros((1, natoms, 3), numpy.float32)
                coor[0, :, :] = temp_mol.coor()[nframes - 1]
                temp_mol.setCoor(coor)
                os.system('rm -f ' + thisdcd)
                temp_mol.write_dcd(thisdcd)
        if ndcdfiles < 1:
            print 'ndcdfiles = ', ndcdfiles
            message = 'Did not save any dcd files.  Decrease dcd write frequency?'
            message + ' :  stopping here'
            print_failure(message, txtOutput)
            sys.exit()

    if not run_fail:
        print '\n> finished dynamics all frames\n'
        print 'writing final data to tamd directory'

        if(infiletype == 'dcd'):
            m1.close_dcd_read(dcdfile[0])

        final_mol = sasmol.SasMol(0)
        final_molw = sasmol.SasMol(1)
        final_mol.read_pdb("tamd_output.pdb", fastread=True)
        final_molw.read_pdb("tamd_output.pdb", fastread=True)

        finaldcdfile = final_molw.open_dcd_write(outfile)

        for i in range(len(dcdlist)):
            print 'i = ', i
            sys.stdout.flush()
            final_mol.read_dcd(dcdlist[i])

            nframes = final_mol.number_of_frames()
            natoms = final_mol.natoms()
            coor = numpy.zeros((1, natoms, 3), numpy.float32)
            coor[0, :, :] = final_mol.coor()[nframes - 1]
            final_molw.setCoor(coor)
            final_molw.write_dcd_step(finaldcdfile, 0, i + 1)

        final_molw.close_dcd_write(finaldcdfile)

        if(keepout != 1):
            rmcmd = 'rm -f '
            for i in range(len(dcdlist)):
                rmcmd = rmcmd + dcdlist[i] + ' '
            print 'rmcmd = ', rmcmd
            os.system(rmcmd)

        txtOutput.put("Total number of frames = %d\n\n" % (nf))
        txtOutput.put("Structures saved to : %s\n" % ('./' + path))
        txtOutput.put("\n%s \n" % (st))

    try:
        if os.path.isfile('cluster_*.str'):
            os.system('mv cluster_*.str ' + path)
    except:
        pass

    os.system('mv temp.inp ' + path)
    os.system('mv tamd.tree ' + path)
    os.system('mv tamd_loops.rst ' + path)
    os.system('mv tamd_output.psf ' + path)
    os.system('mv tamd_output.pdb ' + path)

    time.sleep(2.5)
    print 'TAMD IS DONE'

    return()
