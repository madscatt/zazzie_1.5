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
import sys
import os
import string
import numpy
import time 	
import sasmol.sasmol as sasmol
import sassie.simulate.two_body_grid.poverlap as poverlap
import sassie.simulate.constraints.constraints as constraints

import pprint,copy

#	TWO-BODY GRID
#
#	ADAPTED FROM GRID(IASE)
#
#	12/15/09	--	initial coding 				:	jc	
#	12/16/09	--	hard coded for (integrase)		:	jc
#	08/19/11	--	connected to gui				:	jc
#
#LC	 1         2         3         4         5         6         7
#LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''

		TWO-BODY GRID is the program to move on molecule on a grid
	
		Molecule 1 is the reference molecule.

		Molecule 2 is the molecule to be moved on the grid 
			
		The molecules must have "residue" fields
	
		INPUT:
	
	 	filename1:	Is the PDB file for molecule 1 
	
	 	filename2:	Is the PDB file for molecule 2 (the molecule to be moved)

	 	outfile:	Name of file with the new coordinates for molecule 2
	
	 	txtOutput:	Object to send textual information back to the GUI


		OUTPUT:
	
	 	outfile:	A file with the aligned coordinates
		cartoutfile:	A file with x,y,z,thetax,thetay,thetaz coords of molecule 2

	 	txtOutput: Object text sent to GUI
	
'''

def print_failure(message,txtOutput):

	txtOutput.put("\n\n>>>> RUN FAILURE <<<<\n")
	txtOutput.put(">>>> RUN FAILURE <<<<\n")
	txtOutput.put(">>>> RUN FAILURE <<<<\n\n")
	txtOutput.put(message)
	
	return

def report_status(this_trial,total_number_of_trials,txtOutput):

                #if(((i+1)%(float(trials)/100.0)==0 or (trials<10))):
	fraction_done = (float(this_trial)/float(total_number_of_trials))
#	progress_string='\nCOMPLETED '+str(i+1)+' of '+str(trials)+' : '+str(fraction_done*100.0)+' % done'
#	print('%s\n' % progress_string)
#	print accepted,' configurations accepted out of ',nsteps,(float(accepted)/nsteps)*100.0,' %\n\n'
	report_string='STATUS\t'+str(fraction_done)
	txtOutput.put(report_string)

	return

def unpack_variables(variables):

	runname		=	variables['runname'][0]
	path 			=	variables['path'][0]
	pdbmol1		=	variables['pdbmol1'][0]	
	pdbmol2		=	variables['pdbmol2'][0]
	ofile			=	variables['ofile'][0]
	accpos		=	variables['accpos'][0]	
	pos			=	variables['pos'][0]		
	trans			=	variables['trans'][0]	
	dtrans		=	variables['dtrans'][0]	
	theta			=	variables['theta'][0]	
	dtheta		=	variables['dtheta'][0]	
	basis			=	variables['basis'][0]	
	cutoff		=	variables['cutoff'][0]	
	lowrg			=	variables['lowrg'][0]	
	highrg		=	variables['highrg'][0]	
	zflag			=	variables['zflag'][0]	
	zcutoff		=	variables['zcutoff'][0]	
	cflag			=	variables['cflag'][0]	
	confile		=	variables['confile'][0]	
	nexsegments1	=	variables['nexsegments1'][0]
	nsegments1		=	variables['nsegments1'][0]	
	reslow1		=	variables['reslow1'][0]	
	numcont1		=	variables['numcont1'][0]	
	nexsegments2	=	variables['nexsegments2'][0]
	nsegments2		=	variables['nsegments2'][0]	
	reslow2		=	variables['reslow2'][0]	
	numcont2		=	variables['numcont2'][0]	

	return runname,path,pdbmol1,pdbmol2,ofile,accpos,pos,trans,dtrans,theta,dtheta,basis,cutoff,lowrg,highrg,zflag,zcutoff,cflag,confile,nexsegments1,nsegments1,reslow1,numcont1,nexsegments2,nsegments2,reslow2,numcont2

def euler_rotation(m,phi,theta,psi):
   c11 = numpy.cos(theta)*numpy.cos(psi)
   c12 = numpy.cos(phi)*numpy.sin(psi) + numpy.sin(phi)*numpy.sin(theta)*numpy.cos(psi)
   c13 = numpy.sin(phi)*numpy.sin(psi) - numpy.cos(phi)*numpy.sin(theta)*numpy.cos(psi)
   c21 = -numpy.cos(theta)*numpy.sin(psi)
   c22 = numpy.cos(phi)*numpy.cos(psi) - numpy.sin(phi)*numpy.sin(theta)*numpy.sin(psi)
   c23 = numpy.sin(phi)*numpy.cos(psi) + numpy.cos(phi)*numpy.sin(theta)*numpy.sin(psi)
   c31 = numpy.sin(theta)
   c32 = -numpy.sin(phi)*numpy.cos(theta)
   c33 = numpy.cos(phi)*numpy.cos(theta)
   C = numpy.matrix([[c11,c12,c13],[c21,c22,c23],[c31,c32,c33]])
   newcoor = numpy.array(m.coor()[0]*C)
   m.setCoor(numpy.array([newcoor]))


def molgrid(m1,m2,m3,ofile,genpaths,accpos,pos,trans,dtrans,theta,dtheta,cutoff,basis,mask_array1,mask_array2,zflag,zcutoff,cflag,mask_a_array,mask_b_array,distance_array,type_array,txtOutput):

	frame=0

	m1.calccom(frame) ; com1=m1.com()
	m2.calccom(frame) ; com2=m2.com()
	m3.calccom(frame) ; com3=m3.com()
	
	mass1=m1.mass() ; totalmass1=m1.totalmass()
	mass2=m1.mass() ; totalmass2=m2.totalmass()
	mass3=m3.mass() ; totalmass3=m3.totalmass()

	temp_m3 = sasmol.SasMol(4)
	error = temp_m3.merge_two_molecules(m1,m2)

	print 'total mass1 = ',totalmass1
	print 'total mass2 = ',totalmass2
	print 'total mass3 = ',totalmass3

	print 'mass 1 + mass 2 = ',totalmass1+totalmass2

	numx=trans[0] ; numy=trans[1] ; numz=trans[2]
	delx=dtrans[0] ; dely=dtrans[1] ; delz=dtrans[2]
	
	if(accpos == 1):
		lowx=pos[0] ; lowy=pos[1] ; lowz=pos[2]
	else:
		lowx=com2[0] ; lowy=com2[1] ; lowz=com2[2]	

	numthx=theta[0] ; numthy=theta[1] ; numthz=theta[2]
	delthx=dtheta[0] ; delthy=dtheta[1] ; delthz=dtheta[2]
	
	lowthx = 0.0 ; lowthy = 0.0 ; lowthz = 0.0

	error,coor1 = m1.get_coor_using_mask(frame,mask_array1)
	rm1=sasmol.SasMol(3)
	rm1.setCoor(coor1)

	error,coor2 = m2.get_coor_using_mask(frame,mask_array2)
	rm2=sasmol.SasMol(4)
	rm2.setCoor(coor2)

	error = m1.copy_molecule_using_mask(rm1,mask_array1,frame)
	error = m2.copy_molecule_using_mask(rm2,mask_array2,frame)

	natoms1=m1.natoms() ; natoms2=m2.natoms()
	totalatoms=natoms1+natoms2
	coor=numpy.zeros((1,totalatoms,3),numpy.float)

	for i in range(natoms1):
		coor[frame,i:]=m1.coor()[frame,i,:]	

	cartoutfile = genpaths+'gridrun_output.txt'
	cartout=open(cartoutfile,'w')
	cartout.write("# COM positions and angular values FOR ACCEPTED STRUCTURES\n")
	cartout.write("# x	y	z	thetax	thetay	thetaz\n")
	rad2deg=180.0/numpy.pi
	failtally=0 ; ntrials=0; cyclesum=0.0 ; cycles=0
	start=time.time()    ###	TIME	Timing

	taccepted = 0

	dcdoutfile = m3.open_dcd_write(genpaths+ofile)

	total_number_of_trials = numx*numy*numz*numthx*numthy*numthz

#	check dcd file size

	temporary_dcdoutfile = m3.open_dcd_write('/tmp/temp.dcd')

	m3.write_dcd_step(temporary_dcdoutfile,0,1)
	
	file_size_mb = (os.path.getsize("/tmp/temp.dcd")/1048576.)*total_number_of_trials
	file_size_gb = (os.path.getsize("/tmp/temp.dcd")/1073741824.)*total_number_of_trials

	m3.close_dcd_write(temporary_dcdoutfile)
	rmst = 'rm -f /tmp/temp.dcd'
	os.system(rmst)

	print '\nTOTAL NUMBER OF TRIALS = ',total_number_of_trials

	print '\nMAXIUMUM FINAL DCD FILE SIZE (MB) = ',file_size_mb	
	print 'MAXIUMUM FINAL DCD FILE SIZE (GB) = ',file_size_gb

	print '\n\nSTARTING GRID MOVES\n\n'
	for i in xrange(numx):
		print '\nstarting i = ', i+1,' of ',numx
		print 'ntrials = ',ntrials
		print 'naccepted = ',taccepted
		print 'nfail = ',failtally,'\n'
		if(ntrials>0):
			print '\npercent acceptance = ',((ntrials-failtally)/float(ntrials))*100.0,'\n'					
		tx=lowx+delx*i	
		for j in xrange(numy):
			ty=lowy+dely*j
			for k in xrange(numz):
				if(ntrials==0):
					lasttime=start
#
####	OPEN	Timing of overlap functions
#
				now=time.time()    ###	TIME	Timing
				cycles=cycles+1
				cyclesum=cyclesum+(now-lasttime)
				#print "time = ",now-lasttime, ": avg = ",cyclesum/cycles
				lasttime=now
				tz=lowz+delz*k
				thxcount=0 ; thisthx=0 ; thisthy=0 ; thisthz=0
				for angx in xrange(numthx):	
					axis='x'
					if(thxcount==0):
						theta=lowthx	
						thxcount=1
					else:
						theta=delthx
					thisthx=thisthx+theta
					m2.center(frame)	
					m2.rotate(frame,axis,theta)
					m2.moveto(frame,[tx,ty,tz])
					thycount=0 
					for angy in xrange(numthy):	
						axis='y'
						if(thycount==0):
							theta=lowthy	
							thycount=1
						else:
							theta=delthy
						thisthy=thisthy+theta
						m2.center(frame)	
						m2.rotate(frame,axis,theta)
						m2.moveto(frame,[tx,ty,tz])
						thzcount=0 
						for angz in xrange(numthz):	
							ntrials=ntrials+1
							check=0
							axis='z'
							if(thzcount==0):
								theta=lowthz	
								thzcount=1
							else:
								theta=delthz
							thisthz=thisthz+theta
							m2.center(frame)	
							m2.rotate(frame,axis,theta)
							m2.moveto(frame,[tx,ty,tz])
							report_status(ntrials,total_number_of_trials,txtOutput)
							error,coor2 = m2.get_coor_using_mask(frame,mask_array2)
							rm2.setCoor(coor2)

	#						#check=poverlap.aboverlap(rm1.coor(),rm2.coor(),rm1.natoms(),rm2.natoms(),self._cut)
	#						#check=poverlap.saboverlap(rm1.coor(),rm2.coor(),self._cut)
							check=poverlap.faboverlap(rm1.coor(),rm2.coor(),cutoff)
	
							if(check==0 and zflag==1):
								zee=rm2.coor[0,:,2]
								zcheck=numpy.alltrue(numpy.greater_equal(zee,zval))
								if(zcheck==0):
									check=1
							if(check==0):
								for jj in range(natoms2):
									coor[frame,jj+natoms1,:]=m2.coor()[frame,jj,:]
								temp_m3.setCoor(coor)
								con_check=0		
								if(cflag==1):
									con_check = constraints.check_constraints(temp_m3,mask_a_array,mask_b_array,distance_array,type_array)
								if(con_check==0):	
									m3.setCoor(coor)		
									m3.write_dcd_step(dcdoutfile,0,taccepted+1)
									taccepted += 1
									cartout.write("%f\t%f\t%f\t%f\t%f\t%f\n" % (tx,ty,tz,thisthx*rad2deg,thisthy*rad2deg,thisthz*rad2deg))
									cartout.flush()
								else:
									failtally=failtally+1
							else:
								failtally=failtally+1

	report_status(ntrials,total_number_of_trials,txtOutput)
	time.sleep(1)
	
	#print '\n\nFINAL STATS\n\n'
	#print 'ntrials = ',ntrials
	#print 'naccepted = ',taccepted
	#print 'nfail = ',failtally

	if(taccepted>0):
	    txtOutput.put('percent acceptance = '+str(((ntrials-failtally)/float(ntrials))*100.0))
	    txtOutput.put("\nConfigurations and statistics saved in %s directory\n\n" % ('./'+genpaths))

	else:
	    txtOutput.put('percent acceptance = '+str(((ntrials-failtally)/float(ntrials))*100.0))
	    txtOutput.put("\n NO ACCEPTED MOVES\n\n Statistics saved in %s directory\n\n" % (genpaths))
                                            
	lineintxtOutput=''.join(['=' for x in xrange(60)])
	txtOutput.put("\n%s \n" %(lineintxtOutput))
	time.sleep(1.0)
                                            	
	cartout.close()
	m3.close_dcd_write(dcdoutfile)

	time.sleep(2)
	return

def two_body_grid(variables,txtOutput):

	runname,path,pdbmol1,pdbmol2,ofile,accpos,pos,trans,dtrans,theta,dtheta,basis,cutoff,lowrg,highrg,zflag,zcutoff,cflag,confile,nexsegments1,nsegments1,reslow1,numcont1,nexsegments2,nsegments2,reslow2,numcont2 = unpack_variables(variables)

	if(runname[-1]=='/'):
		lin=len(runname)
		runname=runname[:lin-1]

	direxist=os.path.exists(runname)
	if(direxist==0):
		os.system('mkdir -p '+runname)

	genpath=runname+'/two_body_grid'
	genpaths=genpath+'/'
	direxist=os.path.exists(genpath)
	if(direxist==0):
		os.system('mkdir -p '+genpath)

	m1=sasmol.SasMol(0)
	m2=sasmol.SasMol(1)
	m3=sasmol.SasMol(2)

	m1.read_pdb(pdbmol1)
	m2.read_pdb(pdbmol2)

	error = m3.merge_two_molecules(m1,m2)

	if(error!=[]):
		print 'ERROR:'+error[0]
		print 'ERROR:'+error[0]
		print 'ERROR:'+error[0]

	m3.write_pdb(genpaths+ofile+'.pdb',0,'w')

	cpst = 'cp '+pdbmol1+' '+genpaths
	os.system(cpst)
	cpst = 'cp '+pdbmol2+' '+genpaths
	os.system(cpst)

	frame=0

	mm1=m1.calcminmax()
	mm2=m2.calcminmax()

	lineintxtOutput=''.join(['=' for x in xrange(60)])
	#ttxt=time.ctime()
	ttxt=time.asctime( time.gmtime( time.time() ) ) 
	txtOutput.put("\n%s \n" %(lineintxtOutput))
	txtOutput.put("DATA FROM RUN: %s \n\n" % (ttxt))
    
	
	'''	
	print 'mm1 = ',mm1
	print 'mm2 = ',mm2
	'''

	# set overlap basis for each molecule

	segment_names_1 = string.split(nsegments1,',')
	segment_names_2 = string.split(nsegments2,',')

	'''
	print 'segment_names_1 = ',segment_names_1
	print 'segment_names_2 = ',segment_names_2
	'''

	if(nexsegments1 > 0):

		for i in xrange(nexsegments1):
			if(i==0):	
				basis1st = "(name[i] == 'CA' and not (segname[i] == '"+segment_names_1[i]+"' and ( resid[i] >= "+str(reslow1[i])+" and resid[i] <= "+str(reslow1[i]+numcont1[i])+")))"
			else:
				basis1st = basis1st + " or (name[i] == 'CA' and not (segname[i] == '"+segment_names_1[i]+"' and ( resid[i] >= "+str(reslow1[i])+" and resid[i] <= "+str(reslow1[i]+numcont1[i])+" )))"

	else:
		basis1st = "name[i] == 'CA'"
	
	'''
	print 'basis1st = ',basis1st
	'''

	if(nexsegments2 > 0):

		for i in xrange(nexsegments2):
			if(i==0):	
				basis2st = "(name[i] == 'CA' and not (segname[i] == '"+segment_names_2[i]+"' and ( resid[i] >= "+str(reslow2[i])+" and resid[i] <= "+str(reslow2[i]+numcont2[i])+")))"
			else:
				basis2st = basis2st + " or (name[i] == 'CA' and not (segname[i] == '"+segment_names_2[i]+"' and ( resid[i] >= "+str(reslow2[i])+" and resid[i] <= "+str(reslow2[i]+numcont2[i])+" )))"

	else:
		basis2st = "name[i] == 'CA'"
	
	'''
	print 'basis2st = ',basis2st
	'''

	error,mask_array1 = m1.get_subset_mask(basis1st)
	error,mask_array2 = m2.get_subset_mask(basis2st)

	if(cflag == 1):
		filter_flag = 0
		error,constraint_basis1_array, constraint_basis2_array,distance_array,type_array = constraints.read_constraints(m3,confile,filter_flag)

		mask_a_array = [] ; mask_b_array = []

		for i in xrange(len(distance_array)):
			print constraint_basis1_array[i]
			print constraint_basis2_array[i]
			print distance_array[i]
			print type_array[i]

			error,local_mask_a_array = m3.get_subset_mask(constraint_basis1_array[i])
			error,local_mask_b_array = m3.get_subset_mask(constraint_basis2_array[i])
	
			mask_a_array.append(local_mask_a_array)	
			mask_b_array.append(local_mask_b_array)	

	else:
		mask_a_array = [] ; mask_b_array = []
		distance_array = [] ; type_array = []

	molgrid(m1,m2,m3,ofile,genpaths,accpos,pos,trans,dtrans,theta,dtheta,cutoff,basis,mask_array1,mask_array2,zflag,zcutoff,cflag,mask_a_array,mask_b_array,distance_array,type_array,txtOutput)
			
	return






