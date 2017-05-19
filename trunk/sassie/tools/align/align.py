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
import string
import locale
import bisect
import time
import numpy
import sasmol.sasmol as sasmol

#       ALIGN
#
#	12/17/2004	--	adapted from trehalose project for gag modeling :	jc
#	12/20/2004	--	align N-terminal of CA in 1L6N and 1E6J		:	jc
#	10/16/2005	--	generic align structures 			:	jc
#	03/11/2009	--	modified to allow for different sized molecules	:	jc
#	01/01/2010	--	added sasmol support				:	jc
#	07/29/2011	--	adapted to use sasop & sasmath classes		:	jc
# LC      1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''
        ALIGN is the module that overlaps molecules from a dcd/pdb file
	onto another molecule over a given basis.  The two molecule types
	do not need to be the same but the number of basis atoms used for
	the overlap do need to be identical.

        This module is called from Align Frames from the main
        GUI through the graphical_align.py script.

	REFERENCE:

	W. Kabsch
    	Acta Crystallog. sect. A  32  922-923  (1976)

    	W. Kabsch
    	Acta Crystallog. sect. A  34  827-828  (1978)

'''


def print_failure(message, txtOutput):

	txtOutput.put("\n\n>>>> RUN FAILURE <<<<\n")
	txtOutput.put(">>>> RUN FAILURE <<<<\n")
	txtOutput.put(">>>> RUN FAILURE <<<<\n\n")
	txtOutput.put(message)

	return


def unpack_variables(variables):

      	runname = variables['runname'][0]
      	path = variables['path'][0]
      	infile = variables['infile'][0]
      	pdbmol1 = variables['pdbmol1'][0]
      	pdbmol2 = variables['pdbmol2'][0]
      	ofile = variables['ofile'][0]
      	basis1 = variables['basis1'][0]
      	basis2 = variables['basis2'][0]
      	lowres1 = variables['lowres1'][0]
      	lowres2 = variables['lowres2'][0]
      	highres1 = variables['highres1'][0]
      	highres2 = variables['highres2'][0]
      	ebasis1 = variables['ebasis1'][0]
      	ebasis2 = variables['ebasis2'][0]
        zflag = variables['zflag'][0]
        zcutoff = variables['zcutoff'][0]

	return runname, path, infile, pdbmol1, pdbmol2, basis1, lowres1, highres1, basis2, lowres2, highres2, ofile, ebasis1, ebasis2, zflag, zcutoff


def write_frame_to_file(m2, nf, frame, alignpath, ofile, outtype, dcdoutfile, intype):

	if(outtype == 'dcd'):
		if(intype == 'dcd'):
			m2.write_dcd_step(dcdoutfile, 0, frame)
		elif(intype == 'pdb'):
			# m2.write_dcd_step(dcdoutfile,frame,frame)
                        if(nf == 1):
			    m2.write_dcd_step(dcdoutfile, 0, frame)
                        else:
			    m2.write_dcd_step(dcdoutfile, frame, frame)

	elif(outtype == 'pdb'):
		if(nf == 1):
			m2.write_pdb(alignpath + ofile, 0, 'w')
		elif(nf > 1 and intype == 'dcd'):
			m2.write_pdb(alignpath + ofile, 0, 'a')
		elif(nf > 1 and intype == 'pdb'):
			m2.write_pdb(alignpath + ofile, frame, 'a')


def check_output_type(ofile):

	if(ofile[-3:] == 'dcd'):
		print '\noutput file will be a DCD file'
		outtype = 'dcd'
	elif(ofile[-3:] == 'pdb'):
		print '\noutput file will be a PDB file'
		outtype = 'pdb'
	else:
		outtype = 'dcd'
		message = 'output filename ' + ofile + \
		    ' needs to end in either ".pdb" (1 frame) or ".dcd" (1 or more frames)\n'
		message += ' :  writing output file as a ' + ofile + '.dcd\n'
		print '\n\n', message, '\n\n'
		ofile = ofile + '.dcd'

	return outtype, ofile


def align(variables, txtOutput):
        '''
        ALIGN is the function to read in variables from GUI input and
       	overlap the molecules in a dcd/pdb file onto the coordinates of
	    a reference pdb structure over a given basis.

                runname: 	            project name
                path:                   input/output filepath
                pdbmol1:                reference pdb (mol 1)
                pdbmol2:                input pdb file (mol 2)
                infile:                 input (pdb or dcd) filename (mol 2)
                basis1:                 basis for molecule 1
                basis2:                 basis for molecule 2
                lowres1:                low residue for overlap molecule 1
                highres1:               high residue for overlap molecule 1
                lowres2:                low residue for overlap molecule 2
                highres2:               high residue for overlap molecule 2
                ebasis1:		            extra basis statement molecule 1
                ebasis2:		            extra basis statement molecule 2
                zflag:                  flag for zcutoff 
                zcutoff:                zcutoff value for infile (frame is excluded if it has ANY atoms with a z-value less than the cutoff)

        OUTPUT:

                files stored in "runname"/align directory

		        ofile:			        output filename
		        ofile*.minmax:		    text file with min & max dimensions

        '''

	runname, path, infile, pdbmol1, pdbmol2, basis1, lowres1, highres1, basis2, lowres2, highres2, ofile, ebasis1, ebasis2, zflag, zcutoff = unpack_variables(
	    variables)

	alignpath = runname + '/align/'
	direxist = os.path.exists(alignpath)
	if(direxist == 0):
		os.system('mkdir -p ' + alignpath)

	print 'runname = ', runname

	minmaxfile = ofile + '.minmax'
	mmfile = open(alignpath + minmaxfile, 'w')

	# ttxt=time.ctime()
	ttxt = time.asctime(time.gmtime(time.time()))

	st = ''.join(['=' for x in xrange(60)])

	txtOutput.put("\n%s \n" % (st))
	txtOutput.put("DATA FROM RUN: %s \n\n" % (ttxt))

	m1 = sasmol.SasMol(0)
	m2 = sasmol.SasMol(1)

	m1.read_pdb(path + pdbmol1, check_zero_coor=True)
	m2.read_pdb(path + pdbmol2, check_zero_coor=True)

	try:
		if(infile[-3:] == 'dcd'):
			dcdfile = m2.open_dcd_read(path + infile)
			nf2 = dcdfile[2]
			intype = 'dcd'
		elif(infile[-3:] == 'pdb'):
			m2.read_pdb(path + infile)
			nf2 = m2.number_of_frames()
			intype = 'pdb'
	except:
		message = 'input filename is a PDB or DCD file but it must end with ".pdb" or ".dcd" '
		message += ' :  stopping here'
		print_failure(message, txtOutput)

	outtype, ofile = check_output_type(ofile)

	dcdoutfile = 'None'
	if(outtype == 'dcd'):
		dcdoutfile = m2.open_dcd_write(alignpath + ofile)

 	txtOutput.put("Total number of frames = %d\n\n" % (nf2))

	mass1 = m1.mass()
	mass2 = m2.mass()

	name1 = m1.name()
	name2 = m2.name()

	basis_filter_1_a = '(name[i] == "' + basis1 + '" and (resid[i] >= ' + \
	                     str(lowres1) + ' and resid[i] <= ' + str(highres1) + '))'
	basis_filter_2_a = '(name[i] == "' + basis2 + '" and (resid[i] >= ' + \
	                     str(lowres2) + ' and resid[i] <= ' + str(highres2) + '))'

	if(ebasis1 == "None"):
		basis_filter_1 = basis_filter_1_a
	else:
		basis_filter_1 = ebasis1 + ' and ' + basis_filter_1_a

	if(ebasis2 == "None"):
		basis_filter_2 = basis_filter_2_a
	else:
		basis_filter_2 = ebasis2 + ' and ' + basis_filter_2_a

	error, mask1 = m1.get_subset_mask(basis_filter_1)
	error, mask2 = m2.get_subset_mask(basis_filter_2)

	sub_m1 = sasmol.SasMol(2)
	error = m1.copy_molecule_using_mask(sub_m1, mask1, 0)

	sub_m2 = sasmol.SasMol(3)
	error = m2.copy_molecule_using_mask(sub_m2, mask2, 0)

	com_sub_m1 = sub_m1.calccom(0)
	sub_m1.center(0)
	coor_sub_m1 = sub_m1.coor()[0]

	minx = []; miny = []; minz = []
	maxx = []; maxy = []; maxz = []

	saved = 0
	for i in xrange(nf2):
		if(intype == 'dcd'):
			m2.read_dcd_step(dcdfile, i)
			m2.center(0)
			minmax = m2.calcminmax_frame(0)
			error, sub_m2.coor = m2.get_coor_using_mask(0, mask2)
			sub_m2.setCoor(sub_m2.coor)
			com_sub_m2 = sub_m2.calccom(0)
			sub_m2.center(0)
			coor_sub_m2 = sub_m2.coor[0]
			m2.align(0, coor_sub_m2, com_sub_m2, coor_sub_m1, com_sub_m1)
		elif(intype == 'pdb'):
			m2.center(i)
			minmax = m2.calcminmax_frame(i)
			error, sub_m2.coor = m2.get_coor_using_mask(i, mask2)
			sub_m2.setCoor(sub_m2.coor)
			com_sub_m2 = sub_m2.calccom(0)
			sub_m2.center(0)
			coor_sub_m2 = sub_m2.coor[0]
			m2.align(i, coor_sub_m2, com_sub_m2, coor_sub_m1, com_sub_m1)

		if zflag:
			coor = m2.coor()[0, :, 2]
			if numpy.alltrue(numpy.greater_equal(coor, zcutoff)):
			    write_frame_to_file(m2, nf2, saved + 1, alignpath,ofile, outtype, dcdoutfile, intype)
			    saved += 1
		else:
			write_frame_to_file(m2, nf2, i + 1, alignpath, ofile,outtype, dcdoutfile, intype)

		minx.append(minmax[0][0]) ; miny.append(minmax[0][1]) ; minz.append(minmax[0][2])
		maxx.append(minmax[1][0]) ; maxy.append(minmax[1][1]) ; maxz.append(minmax[1][2])

		if(((i+1)%(float(nf2)/10.0)==0 or (nf2<10))):
			fraction_done = (float(i+1)/float(nf2))
			progress_string='\nCOMPLETED '+str(i+1)+' of '+str(nf2)+' : '+str(fraction_done*100.0)+' % done'
			print('%s\n' % progress_string)
			report_string='STATUS\t'+str(fraction_done)
			txtOutput.put(report_string)

	if(intype == 'dcd'):
		m2.close_dcd_read(dcdfile[0])
	if(outtype == 'dcd'):
		m2.close_dcd_write(dcdoutfile)

	min_x = numpy.min(minx) ; min_y = numpy.min(miny) ; min_z = numpy.min(minz)
	max_x = numpy.max(maxx) ; max_y = numpy.max(maxy) ; max_z = numpy.max(maxz)

	mmfile.write("%s\n" % ("#min_x,min_y,min_z,max_x,max_y,max_z"))
	mmfile.write("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n" % (min_x,min_y,min_z,max_x,max_y,max_z))
	mmfile.close()
	
	txtOutput.put("minimum x = %lf\t maximum x = %lf -> range: %lf Angstroms\n" % (min_x,max_x,(max_x-min_x)))
	txtOutput.put("minimum y = %lf\t maximum y = %lf -> range: %lf Angstroms\n" % (min_y,max_y,(max_y-min_y)))
	txtOutput.put("minimum z = %lf\t maximum z = %lf -> range: %lf Angstroms\n\n" % (min_z,max_z,(max_z-min_z)))

	if zflag:
		print '\nAligned data (nf=%i) were written to %s\n' % (saved,'./'+alignpath+ofile)
		txtOutput.put('\nAligned data (nf=%i) were written to %s\n' % (saved,'./'+alignpath+ofile))
		if saved == 0:
			print '\nNO FRAMES WERE SAVED TO FILE DUE TO SUPPLIED Z-CUTOFF\n'
			txtOutput.put('\nNO FRAMES WERE SAVED TO FILE DUE TO SUPPLIED Z-CUTOFF\n')

	else:
		print '\nAligned data (nf=%i) were written to %s\n' % (nf2,'./'+alignpath+ofile)
		txtOutput.put('\nAligned data (nf=%i) were written to %s\n' % (nf2,'./'+alignpath+ofile))
	
	print '\nDimension data were written to %s\n' % ('./'+alignpath+minmaxfile)
	txtOutput.put("\nDimension data were written to %s\n" % ('./'+alignpath+minmaxfile))
	txtOutput.put("\n%s \n" %(st))
	# print 'ALIGN2 IS DONE'
	time.sleep(1.5)

	return() 

