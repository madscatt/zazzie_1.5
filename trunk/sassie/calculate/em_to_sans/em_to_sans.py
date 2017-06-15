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
import os,string,os,locale,sys,struct,array,time,platform
try:
    import Gnuplot,Gnuplot.PlotItems, Gnuplot.funcutils
except:
    pass

import numpy
import sassie.util.sasconfig as sasconfig

#       EM_TO_SANS
#
#	12/28/06	--	initial coding			  :	jc
#	03/13/09	--	adapted subroutine to SASSIE GUI  :	jc
#	03/31/09	--	added MRC file type input	  :	jc
#
#
#LC      1         2         3         4         5         6         7
#LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''
        EM_TO_SANS is the module that reads a three-dimensional file
	(either Gaussian Cube or a binary MRC electron density file) and 
	then represent each voxel with an occupation > threshold
	as a scattering center represented as beads (C-atoms), 
	that is then written as PDB/xyz files
	and the SANS profile is calculated using the "scat" program

        This module is called from EM to SANS from the main
        GUI through the graphical_em_to_sans.py script.

'''


def print_failure(message,txtOutput):

        txtOutput.put("\n\n>>>> RUN FAILURE <<<<\n")
        txtOutput.put(">>>> RUN FAILURE <<<<\n")
        txtOutput.put(">>>> RUN FAILURE <<<<\n\n")
        txtOutput.put(message)

        return

def read_mrcfile(filename,ofile,path,threshold,angpix,txtOutput):
	'''
        READ_MRCFILE is the function to read a 3D binary MRC file
	'''

	if(len(angpix)!=3):
		print 'angpix = ',angpix
		print '\nerror in angstrom to pixel resolution input\nstopping\n'
		message='error in angstrom to pixel resolution input: '+str(angpix)
		message+='stopping here'
		print_failure(message,txtOutput)
	
	xangpix=angpix[0] ; yangpix=angpix[1] ; zangpix=angpix[2]


	bit = platform.architecture()[0]

	
	ang2au=1.0/0.5291772108
	au2ang=1.0/ang2au
	outfile=open(path+'/'+ofile,'w')
	outfile2=open(path+'/'+ofile+'.xyz','w')

	input=open(filename,'rb').read()

        if(bit != '64bit'):## @NOTE to ZHL: I do not understand
	    start,stop=0,struct.calcsize('llllllllllllllll')
	    nx,ny,nz,mode,nxstart,nystart,nzstart,mx,my,mz,cella1,cella2,cella3,cellb1,cellb2,cellb3=struct.unpack('llllllllllllllll',input[start:stop])	
        else:## @NOTE to ZHL: I do not understand
            #### FIX HERE ####
	    start,stop=0,struct.calcsize('iiiiiiiiiiffffff')
	    nx,ny,nz,mode,nxstart,nystart,nzstart,mx,my,mz,cella1,cella2,cella3,cellb1,cellb2,cellb3=struct.unpack('iiiiiiiiiiffffff',input[start:stop])	

	print 'start = ',start
	print 'stop = ',stop
	print 'nx = ',nx
	print 'ny = ',ny
	print 'nz = ',nz
	print 'mode/type = ',mode
	print 'nxstart = ',nxstart
	print 'nystart = ',nystart
	print 'nzstart = ',nzstart
	print 'mx = ',mx
	print 'my = ',my
	print 'mz = ',mz
	print 'cella1 = ',cella1
	print 'cella2 = ',cella2
	print 'cella3 = ',cella3
	print 'cellb1 = ',cellb1
	print 'cellb2 = ',cellb2
	print 'cellb3 = ',cellb3
	if(bit != '64bit'): ## @NOTE to ZHL: I do not understand
	    start=1024 ; fsize=struct.calcsize('l')
	else: ## @NOTE to ZHL: I do not understand
	    start=1024 ; fsize=struct.calcsize('f')
	stop=start+fsize
	print 'start = ',start
	print 'stop = ',stop
	nxgp=mx ; nygp=my ; nzgp=mz
	zmin=-mz/2.0 ; ymin=-my/2.0 ; xmin=-mx/2.0
	xgs=(2.0*abs(xmin)/nxgp)*xangpix ; xgs2=xgs/2.0
	ygs=(2.0*abs(ymin)/nygp)*yangpix; ygs2=ygs/2.0
	zgs=(2.0*abs(zmin)/nzgp)*zangpix ; zgs2=zgs/2.0
	zmina=zmin*zangpix ; ymina=ymin*yangpix ; xmina=xmin*xangpix	
	ii=0
	rbin=numpy.zeros((nxgp*nygp*nzgp,3),numpy.float32)
	for i in range(mz):
		thisz=zmina+zgs2+i*zgs
		for j in range(my):
			thisy=ymina+ygs2+j*ygs
			for k in range(mx):
				thisx=xmina+xgs2+k*xgs
				val=struct.unpack('f',input[start:stop])
				start=stop
				stop=stop+fsize
	 			if(val[0]>threshold):
					rbin[ii][0]=thisx	
					rbin[ii][1]=thisy	
					rbin[ii][2]=thisz	
					ii=ii+1
	
	numi=ii

	#st1="ATOM      1  C   CAR     1     "
	st1="ATOM  "
	st2="  C   DUM"
	st3="    "
	st4="  0.00  0.00      DUM"
#ATOM      2  HT1 LYS     1     -13.587 -15.185  12.469  0.32  0.60      7RSA
	outfile.write("REMARK\n")
	for i in range(numi):

		if(i>9999):
			resid = 9999
		else:
			resid = i+1

		if(i>99999):
			index = 99999
		else:
			index = i+1

		outfile.write('%6s%5i%s %5i%s%8.3f%8.3f%8.3f%s\n' % (st1,index,st2,resid,st3,rbin[i][0],rbin[i][1],rbin[i][2],st4))
		outfile2.write('%f\t%f\t%f\t%f\n' % (rbin[i][0],rbin[i][1],rbin[i][2],1.0))
	outfile.close()	
	outfile2.close()	
	
	return

def wait(str=None, prompt='Plot will clear in 10 seconds ...\n'):
	'''
        WAIT is the function to wait for user input to proceed.
	'''
	
	if str is not None:
		print str

	try:
		if(platform.system() == "Linux"):
			import curses
			stdscr = curses.initscr()
			stdscr.addstr('press a key to continue')
			c = stdscr.getch()
			curses.endwin()
	except:
		time.sleep(1)

def read_cubefile(emdensityfile,ofile,path,threshold):
	'''
        READ_CUBEFILE is the function to read a 3D binary Gaussian cube file
	'''
	
	ang2au=1.0/0.5291772108
	au2ang=1.0/ang2au
	
	outfile=open(path+'/'+ofile,'w')
	outfile2=open(path+'/'+ofile+'.xyz','w')
	infile=open(path+'/'+emdensityfile,'r').readlines()

	line0=infile[0]		#calpha_JC CUBE FILE
	line1=infile[1]		#OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z
	line2=infile[2]		#431     -236.215766     -236.215766     -236.215766
	line3=infile[3]		#41      11.338357       0.000000        0.000000
	line4=infile[4]		#41      0.000000        11.338357       0.000000
	line5=infile[5]		#41      0.000000        0.000000        11.338357	
	lin2=string.split(line2)
	lin3=string.split(line3)
	lin4=string.split(line4)
	lin5=string.split(line5)

	natom=locale.atoi(lin2[0])
	xmin=au2ang*locale.atof(lin2[1])	
	ymin=au2ang*locale.atof(lin2[2])	
	zmin=au2ang*locale.atof(lin2[3])	

	nxgp=locale.atoi(lin3[0])
	nygp=locale.atoi(lin4[0])
	nzgp=locale.atoi(lin5[0])

	print 'natom = ',natom
	print 'xmin = ',xmin,' ymin = ',ymin,' zmin = ',zmin
	print 'nxgp = ',nxgp,' nygp = ',nygp,' nzgp = ',nzgp 

	cbin=numpy.zeros((nxgp,nygp,nzgp),numpy.float32)
	rbin=numpy.zeros((nxgp*nygp*nzgp,3),numpy.float32)
	
	xgs=(2.0*abs(xmin)/nxgp) ; xgs2=xgs/2.0
	ygs=(2.0*abs(ymin)/nygp) ; ygs2=ygs/2.0
	zgs=(2.0*abs(zmin)/nzgp) ; zgs2=zgs/2.0

	print 'xgs = ',xgs,' ygs = ',ygs, ' zgs = ',zgs

	skip=natom+6 
	line=infile[skip] ; thisline=skip
	lin=string.split(line) 
	ii=0	
        for i in range(nxgp):
		thisx=xmin+xgs2+i*xgs
                for j in range(nygp):
			thisy=ymin+ygs2+j*ygs
			loc=0
                        for k in range(nzgp):
				thisz=zmin+ygs2+k*zgs
				l=len(lin)
				cbin=locale.atof(lin[loc])
				stop=0	
				if(i==nxgp-1 and j==nygp-1 and k==nzgp-1):
					stop=1
				if(((k%6)==5 or k==nzgp-1) and stop==0) :
                                	thisline=thisline+1
					line=infile[thisline]	
					lin=string.split(line)
					loc=0
				else:
					loc=loc+1

	 			if(cbin>threshold):
					rbin[ii][0]=thisx	
					rbin[ii][1]=thisy	
					rbin[ii][2]=thisz	
					ii=ii+1

	numi=ii

	#st1="ATOM      1  C   CAR     1     "
	st1="ATOM  "
	st2="  C   DUM"
	st3="    "
	st4="  0.00  0.00      DUM"
#ATOM      2  HT1 LYS     1     -13.587 -15.185  12.469  0.32  0.60      7RSA
	outfile.write("REMARK\n")
	for i in range(numi):
		outfile.write('%6s%5i%s %5i%s%8.3f%8.3f%8.3f%s\n' % (st1,i+1,st2,i+1,st3,rbin[i][0],rbin[i][1],rbin[i][2],st4))
		outfile2.write('%f\t%f\t%f\t%f\n' % (rbin[i][0],rbin[i][1],rbin[i][2],1.0))
	outfile.close()	
	outfile2.close()	
	
	return

def write_inputfile(inp,ofile,sansfile,npoints,qmax):
	'''
        WRITE_INPUTFILE is the function to write the input file for scat
	'''
	
	outfile=open(inp,'w')
	outfile.write('%i\n%s\n%s\n%i\n%f\n' % (1,ofile,sansfile,npoints,qmax))
	outfile.close()	

	return

def run_scat(inp):
	'''
        RUN_SCAT is the function to run the scat program
	'''
        bin_path = sasconfig.__bin_path__+os.sep
	
	currcommand=bin_path+'scat.exe < '+inp
	#currcommand='/usr/local/bin/scat.exe < '+inp
	os.system(currcommand)

def read_dat(ifile):
	'''
        READ_DAT is the function to read data from the input file
	'''

	infile=open(ifile,'r').readlines()
	nl=len(infile)
	data=[]
	for i in range(nl):
		lin=string.split(infile[i])
		x=locale.atof(lin[0]) ; y=locale.atof(lin[1]) ; z=locale.atof(lin[2])
		data.append([x,y,z])

	return data

def unpack_variables(variables):

        runname         = variables['runname'][0]
        emfiletype      = variables['emfiletype'][0]
        inputpath       = variables['inputpath'][0]
        emdensityfile   = variables['emdensityfile'][0]
        angpix          = variables['angpix'][0]
        pdbfile         = variables['pdbfile'][0]
        threshold       = variables['threshold'][0]
        sansfile        = variables['sansfile'][0]
        npoints         = variables['npoints'][0]
        qmax            = variables['qmax'][0]
        plotflag 		= variables['plotflag'][0]

	return runname,emfiletype,inputpath,emdensityfile,angpix,pdbfile,threshold,sansfile,npoints,qmax,plotflag

def em(variables,txtOutput):

        '''
        EM is the function to read in three-dimensional voxel data and
	calculate the SANS profile through the binary program scat.

        INPUT:  variable descriptions:

                runname:        project name
                emfiletype:     em file type (0=cube,1=mrc)
                inputpath:      input file path
                emdensityfile:  em filename
                angpix:         mrc angstroms/pixel resolution (x,y,z)
                threshold:      treshold cutoff
                npoints:        number of points in sans calculation
                qmax:           q-max in sans calculation

        OUTPUT:

                pdbfile:        output filename (pdb)
                
                files stored in ./"runname"/em_to_sans/ directory:

                sansfile*.sub:  output sans profile
                sansfile*.pr:   output p(r)
                dum.inp:        input file written for scat
                sansfile*.pdb:  pdb file of coordinates used for scat
                sansfile*.xyz:  xyz file of coordinates used for scat
        '''
	runname,emfiletype,inputpath,emdensityfile,angpix,pdbfile,threshold,sansfile,npoints,qmax,plotflag=unpack_variables(variables)

        empath=runname+'/em_to_sans'
        direxist=os.path.exists(empath)
        if(direxist==0):
                try:
                        result=os.system('mkdir -p '+empath)
                except:
                        message='can not create project directory: '+empath
                        message+='\nstopping here\n'
                        print_failure(message,txtOutput)
                if(result!=0):
                        message='can not create project directory: '+empath
                        message+='\nstopping here\n'
                        print_failure(message,txtOutput)

        print 'runname = ',runname
        
	#ttxt=time.ctime()
	ttxt=time.asctime( time.gmtime( time.time() ) ) 
        st=''.join(['=' for x in xrange(60)])

        txtOutput.put("\n%s \n" %(st))
        txtOutput.put("DATA FROM RUN: %s \n\n" %(ttxt))

	inp='dum.inp'
	if(emfiletype==0):
		read_cubefile(emdensityfile,pdbfile,inputpath,threshold)
	elif(emfiletype==1):
		read_mrcfile(emdensityfile,pdbfile,inputpath,threshold,angpix,txtOutput)
	else:
		print 'wrong file type entered'
		print '0==Gaussian Cube file'
		print '1==MRC file'
		print 'stopping now'
                message='wrong file type entered: '+emfiletype
                message+='\nstopping here\n'
                print_failure(message,txtOutput)


	fraction_done = 0.5 
	report_string='STATUS\t'+str(fraction_done)
	txtOutput.put(report_string)

	write_inputfile(inp,pdbfile+'.xyz',sansfile,npoints,qmax)
	run_scat(inp)
	
	fraction_done = 1.0 
	report_string='STATUS\t'+str(fraction_done)
	txtOutput.put(report_string)


	if(plotflag == 1):
		graph = Gnuplot.Gnuplot(debug=1)
        	graph.clear()
        	graph('set title "SANS Profile"')
        	graph.xlabel('Q (1/A)')
        	graph.ylabel('I(Q)')
        	graph('set logscale x')
        	graph('set logscale y')

		iqdat=read_dat(sansfile)

		graph.plot(Gnuplot.Data(iqdat,using='1:2 w lp ps 4',title='Calculated SANS Profile'))

	for i in range(len(sansfile)):
		char=sansfile[i]
		if(char=='.'):
			pos=i
	prfile=sansfile[0:pos]+'.pr'
	print 'prfile = ',prfile

	prdat=read_dat(prfile)

	print 'pdbfile = ',pdbfile
	mvst='mv -f '+prfile+' '+pdbfile+' dum.inp '+sansfile+' '+pdbfile+'.xyz '+empath
	print 'mvst = ',mvst
	os.system(mvst)

        txtOutput.put("Data stored in directory: %s\n\n" % ('./'+empath))
        txtOutput.put("\n%s \n" %(st))
	time.sleep(1)

	if(plotflag == 1):
        	graph2 = Gnuplot.Gnuplot(debug=1)
        	graph2.clear()
        	graph2('set title "P(r) Profile"')
        	graph2.xlabel('r (Angstrom)')
        	graph2.ylabel('P(r)')

        	graph2.plot(Gnuplot.Data(prdat,using='1:2 w lp ps 4',title='Calculated P(r) Profile'))
	        
	if(plotflag == 1):
		wait('\n')

	print 'EM TO SANS IS DONE'
	
	return

