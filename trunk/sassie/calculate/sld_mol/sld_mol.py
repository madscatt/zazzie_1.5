'''
    SASSIE: Copyright (C) 2011 Joseph E. Curtis, Ph.D., Hirsh Nanda, Ph.D.

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

import sys,string,math,os,time,multiprocessing,time
import numpy
import scipy.optimize
import sassie.calculate.sld_mol.modified_scipy_interpolate as interpolate
import sasmol.sasmol as sasmol
import Gnuplot,Gnuplot.PlotItems, Gnuplot.funcutils
#import matplotlib.pyplot as plt

#       SLD_MOL
#
#	11/10/2004	--	initial coding				:	jc
#	12/10/2004	--	prototype for alpha sassie		:	jc
#	12/06/2011	--	re-written with new volumetric code	:	hn
#	12/07/2011	--	added SASMOL support			:	jc
#
#LC      1         2         3         4         5         6         7
#LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''
	SLD_MOL is the module that calculates the scattering length density profile 
	from a dcd/pdb file

	This module is called from Scattering Length Density Profile from the main
	GUI through the graphical_sld_mol.py script.

	This method compares an experimentally derived SLD profile with heavy atom 
	distribution from a pdb or dcd file containing protein structure(s).
	It performs a fit allowing a normalization factor, z-pos & constant shift.
	The SLD profile is convolved with a gauss of user defined width to mimic instrument
	resolution and roughness and outputs a text file with frame number and fit_error.
'''

def print_failure(message,txtOutput):

	txtOutput.put("\n\n>>>> RUN FAILURE <<<<\n")
	txtOutput.put(">>>> RUN FAILURE <<<<\n")
	txtOutput.put(">>>> RUN FAILURE <<<<\n\n")
	txtOutput.put(message)

	return

def unpack_variables(variables):

	runname 	= variables['runname'][0]
	path		= variables['path'][0]             
	pdbfile	= variables['pdbfile'][0]         
	dcdfile	= variables['dcdfile'][0]        
	expdatafile	= variables['expdatafile'][0]
	outputfile	= variables['outputfile'][0]

	runtype	= variables['runtype'][0]
	bulk_sld	= variables['bulk_sld'][0]
	xon		= variables['xon'][0]
	
	num_deut_regions= variables['num_deut_regions'][0]
	deut_low_res	= variables['deut_low_res'][0]
	deut_high_res	= variables['deut_high_res'][0]
	
	dbin		= variables['dbin'][0]
	width		= variables['width'][0]

	sldfit	= variables['sldfit'][0]		
	sldoffset	= variables['sldoffset'][0]		

	zfit0		= variables['zfit0'][0]
	zfitmin	= variables['zfitmin'][0]
	zfitmax	= variables['zfitmax'][0]
	zevalmin	= variables['zevalmin'][0]
	zevalmax	= variables['zevalmax'][0]
	A0		= variables['A0'][0]	
	Amin		= variables['Amin'][0]
	Amax		= variables['Amax'][0]
	plotflag	= variables['plotflag'][0]
	
	return runname,path,pdbfile,dcdfile,expdatafile,outputfile,plotflag,runtype,bulk_sld,xon,num_deut_regions,deut_low_res,deut_high_res,dbin,width,sldfit,sldoffset,zfit0,zfitmin,zfitmax,zevalmin,zevalmax,A0,Amin,Amax 

def plot_me(plot_type,plot_variables):
	
	if(plot_type == 'gnuplot'):
		graph = plot_variables[0]
		sld_data = plot_variables[1]
		graph.plot(Gnuplot.Data(sld_data,using='1:2 w p ps 1',title='SLD'))

	elif(plot_type == 'matplotlib'):
		pass

	return

def read_experimental_sld_file(path,expdatafile,plotflag,graph,sldoffset):

	z_value=[]
	sld_value=[]
	variance_value = []
	inputfile = open(expdatafile,'r').readlines()

	test_dat = []
	
	for i in range(len(inputfile)):
		lin = string.split(inputfile[i])
		if(lin[0][0] != "#" and len(lin) >= 2):
			z_value.append(float(lin[0])+sldoffset); sld_value.append(float(lin[1]))
			test_dat.append([float(lin[0])+sldoffset,float(lin[1])])
			if(len(lin) == 3):
				variance_value.append(float(lin[2]))

	print 'interpolating data'
	
	interpolated_data = interpolate.interp1d(z_value,sld_value,'cubic')

	if(len(variance_value) == len(sld_value)):	
		interpolated_variance_data = interpolate.interp1d(z_value,variance_value,'cubic')
	else:
		interpolated_variance_data = None

	fnew = interpolated_data(z_value)

	if(plotflag == 1):
		plt.plot(z_value,sld_value,'-', z_value,fnew,'o')
		raw_input('Please press return to continue...\n')
	elif(plotflag == 2):
		plot_type = 'gnuplot'
		plot_variables = [graph,test_dat]	
		plot_me(plot_type,plot_variables)

	return z_value,interpolated_data,interpolated_variance_data

def update_coor(m1,residues,resids_mask,frame):

	error,new_coor=m1.get_coor_using_mask(frame,resids_mask)		

	residues.setCoor(new_coor)

	return

def print_status(step,total_steps,txtOutput):

	fraction_done = (float(step+1)/float(total_steps))
	progress_string='\nCOMPLETED '+str(step+1)+' of '+str(total_steps)+' : '+str(fraction_done*100.0)+' % done'
	print('%s\n' % progress_string)
	report_string='STATUS\t'+str(fraction_done)
	txtOutput.put(report_string)
	
	return

def residue_com(m1,residues,frame):

	all_residue_com = []

	for i in xrange(m1.number_of_resids()):
		resids_mask = m1.resids_mask()[i]
		update_coor(m1,residues[i],resids_mask,frame)
		com = residues[i].calccom(0)
		all_residue_com.append(com)

	return numpy.array(all_residue_com)

def setup_sld_calculation(m1,xon,deuterated_residues):
	'''
		xon: xray (xon=1) or neutron (xon=0)
		deuterated_residues: list of resids for deuterated residues
	'''

	amino_acid_properties = m1.amino_acid_sld()

	resid = m1.resid()

	m1.initialize_children()

	residues = m1.init_child("resids")

	number_of_resids = m1.number_of_resids()

	isotopic_residues = numpy.zeros(number_of_resids,numpy.short)
  
	residue_sld_hydrogen = numpy.zeros(number_of_resids)
	residue_sld_deuterium = numpy.zeros(number_of_resids)
	residue_sld_electron = numpy.zeros(number_of_resids)
	volumes = numpy.zeros(number_of_resids)
	radius = numpy.zeros(number_of_resids)
 
	hydrogen_scattering_length = -3.7409
	deuterium_scattering_length = 6.671

	first_residue = resid[0]

	for i in deuterated_residues: 
		isotopic_residues[i-first_residue]=1

	for i in range(number_of_resids):
		residue_name=residues[i].resname()[0]
		residue_value=amino_acid_properties[residue_name]
		volumes[i]= residue_value[0]
		if(xon==1):   
			residue_sld_electron[i] = residue_value[1]
		elif(xon==0):
			residue_sld_hydrogen[i] = (residue_value[2+isotopic_residues[i]]+residue_value[4]*hydrogen_scattering_length)
			residue_sld_deuterium[i] = (residue_value[2+isotopic_residues[i]]+residue_value[4]*deuterium_scattering_length)

	radius = numpy.power(volumes*3/(4*math.pi),1.0/3.0)

	return radius,volumes,residue_sld_hydrogen,residue_sld_deuterium,residue_sld_electron,residues
 
def calc_sld(m1,dcdinputfile,frame,dbin,xon,runtype,minmaxz,radius,volumes,residue_sld_hydrogen,residue_sld_deuterium,residue_sld_electron,residues,txtOutput):
	'''
		frames: list of frame numbers to calculate SLD
		xon: xray (xon=1) or neutron (xon=0)
		runtype: 0 == average sld over entire ensemble, 1 == individual sld
		radius, volumes, residue_sld_hydrogen, residue_sld_deuterium, 
		residue_sld_electron, residues all generated in setup_sld_calculation
		txtOutput object is required to communicate with GUI
	'''

	number_of_resids = m1.number_of_resids()
	com_zpos = numpy.array(()) #Becomes a 2D array [frame,resid]
	
	minz = minmaxz[0]
	maxz = minmaxz[1]
           
	binnum = int((maxz-minz)/dbin)+1

	slprofileE = numpy.zeros((binnum))
	slprofileH = numpy.zeros((binnum))
	slprofileD = numpy.zeros((binnum))
	volprofile = numpy.zeros((binnum))

	if(xon==0):   
		fact=1.0E-5
	elif(xon==1):
		fact=2.8E-5

	if(runtype == 1):
		number_of_frames = 1
	else:
		number_of_frames = dcdinputfile[2]

	if(runtype == 0):
		message = 'CALCULATING SLD FOR STRUCTURES' 
		print "\n",message	
		txtOutput.put("\t"+message)

	for f in xrange(number_of_frames):

		if(runtype == 1):
			m1.read_dcd_step(dcdinputfile,frame)
		else:
			m1.read_dcd_step(dcdinputfile,f)

		all_residue_com = residue_com(m1,residues,0)
		
		com_zpos_tmp = numpy.zeros((number_of_resids))
	
		for i in range(number_of_resids):
			com_zpos_tmp[i] = all_residue_com[i][2]

		com_zpos=numpy.append(com_zpos,com_zpos_tmp)
		
		if(runtype == 0):
			if((((f+1.0)%(number_of_frames/100.0))==0) or (f==number_of_frames-1) or (f==0)):
				print_status(f,number_of_frames,txtOutput)
	
	com_zpos=numpy.reshape(com_zpos,(-1,number_of_resids))

	if(runtype == 0):
		message = 'FILLING PROFILE ARRAYS' 
		print "\n",message	
		txtOutput.put("\t"+message)

	for f in xrange(number_of_frames):
		if(runtype == 0):
			if((((f+1.0)%(number_of_frames/100.0))==0) or (f==number_of_frames-1) or (f==0)):
				print_status(f,number_of_frames,txtOutput)

		for i in numpy.arange(binnum):
			local_a  = dbin*math.pi*(radius**2-((minz+i*dbin)-com_zpos[f])**2)
			local_a = local_a.clip(min=0.0)
			volprofile[i] = numpy.sum(local_a)
			if(xon==1):
				slprofileE[i] = numpy.sum(residue_sld_electron*(local_a/volumes))
			elif(xon==0):
				slprofileH[i] = numpy.sum(residue_sld_hydrogen*(local_a/volumes))
				slprofileD[i] = numpy.sum(residue_sld_deuterium*(local_a/volumes))
		
	if(xon==0): 
		return volprofile/number_of_frames, slprofileH*fact/number_of_frames, slprofileD*fact/number_of_frames
	elif(xon==1): 
		return volprofile/number_of_frames, slprofileE*fact/number_of_frames, slprofileD*fact/number_of_frames

def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string   
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"

    s=numpy.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]

    if window == 'flat': #moving average
        w=numpy.ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')

    y=numpy.convolve(w/w.sum(),s,mode='valid')
    return y

def smooth_spectra(width,hvydist,dbin):

	newdist = smooth(hvydist,10,'blackman')

	return numpy.r_[newdist]

def calc_error(lFitParam, dbin, zmaxprot, zeval, fnvfprot, fnsldprot, bulk_sld, interpolated_experimental_sld, interpolated_variance_data, txtOutput):
	zfit, Afit= lFitParam
	err=[]
	zevalmin = zeval[0]; zevalmax = zeval[1]
	zmaxpr = zfit+zmaxprot-dbin
	npoints = 0
	for z in numpy.arange(zevalmin,zevalmax,dbin):
		if zfit>z:
			if(interpolated_variance_data == None):
				err.append(abs(interpolated_experimental_sld(z)-bulk_sld)*1e6)
			else:
				#err.append(((interpolated_experimental_sld(z)-SLDz)**2)/(interpolated_variance_data(z)**1))
				err.append(((interpolated_experimental_sld(z)-bulk_sld)**2)/(interpolated_variance_data(z)**1))
				npoints += 1
		elif zfit<=z and zmaxpr>z:
			vfprot = fnvfprot(z-zfit); sldprot = fnsldprot(z-zfit)
			SLDz = Afit*(vfprot*sldprot) + Afit*(1-vfprot)*bulk_sld+(1-Afit)*bulk_sld
			if(interpolated_variance_data == None):
				err.append(abs(interpolated_experimental_sld(z)-SLDz)*1e6)
			else:
				#err.append(((interpolated_experimental_sld(z)-SLDz)**2)/(interpolated_variance_data(z)**1))
				err.append(((interpolated_experimental_sld(z)-bulk_sld)**2)/(interpolated_variance_data(z)**1))
				npoints += 1
		elif zmaxpr<=z:
			if(interpolated_variance_data == None):
				err.append(abs(interpolated_experimental_sld(z)-bulk_sld)*1e6)
			else:
				#err.append(((interpolated_experimental_sld(z)-SLDz)**2)/(interpolated_variance_data(z)**1))
				err.append(((interpolated_experimental_sld(z)-bulk_sld)**2)/(interpolated_variance_data(z)**1))
				npoints += 1
		else:
			message = 'failed in error calculation\n'+str(zfit)
			print_failure(message,txtOutput)


	if(npoints>1):
		sum_error = sum(err)/(npoints-1)
	else:
		sum_error = sum(err)
			
	print zfit, Afit, sum_error

	return sum_error

def plot_and_save_sld_fit(param, dbin, zexp, interpolated_experimental_sld, vfprot, sldprot, bulk_sld, plotflag, graph2, sum_sld_fit, this_frame, number_of_frames, runtype, fit_error, results, sld_output_files, sldpath):
	z, A= param
	SLDfit=A*(vfprot*sldprot)+A*(1-vfprot)*bulk_sld+(1-A)*bulk_sld
	SLDexp=interpolated_experimental_sld(zexp)

	if(plotflag == 1):

		plt.figure(3)
		plt.plot(zexp,SLDexp,'-',numpy.arange(len(SLDfit))*dbin+z,SLDfit)
		raw_input('Please press return to continue...\n')

	exp_data = []	
	for i in xrange(len(zexp)):
		exp_data.append([zexp[i],SLDexp[i]])	
	
	fit_data = []
	avg_sld_fit = []
	for i in xrange(len(SLDfit)):
		this_z = i*dbin+z
		fit_data.append([this_z,SLDfit[i]])
		sum_sld_fit[i] += SLDfit[i]
		avg_sld_fit.append([this_z,(sum_sld_fit[i]/this_frame)])
	
	if(plotflag == 2):
		if(runtype == 0):	
			graph2.plot(Gnuplot.Data(exp_data,using='1:2 w p ps 1',title='experimental SLD'),Gnuplot.Data(fit_data,using='1:2 w p ps 1',title='fit SLD')) 
		elif(runtype == 1):
			graph2.plot(Gnuplot.Data(exp_data,using='1:2 w p ps 1',title='experimental SLD'),Gnuplot.Data(fit_data,using='1:2 w p ps 1',title='fit SLD'),Gnuplot.Data(avg_sld_fit,using='1:2 w p ps 1',title='average fit SLD')) 

	outputfile = open(sld_output_files[this_frame-1],'w')
	average_outputfile = open(sldpath+'average_sld_fit.txt','w')

	st = "#  FRAME = "+str(this_frame)+": OPTIMUM Z POSITION = "+str(results[0][0])+" : OPTIMUM SURFACE COVERAGE = "+str(results[0][1])+" : FIT ERROR = "+str(fit_error)+"\n"

	outputfile.write(st)
	average_outputfile.write(st)

	st = "# z (angstroms) : SLD\n"

	outputfile.write(st)
	average_outputfile.write(st)
	
	for i in xrange(len(SLDfit)):
		outputfile.write("%f\t%e\n" % (fit_data[i][0],fit_data[i][1]))
		average_outputfile.write("%f\t%e\n" % (avg_sld_fit[i][0],avg_sld_fit[i][1]))			

	outputfile.close()
	average_outputfile.close()

	return 

def run_file_utilities(runname,number_of_frames):

	direxist=os.path.exists(runname)
	if(direxist==0):
		os.system('mkdir -p '+runname+'/')
#
#       global run administration
#
	sldpath=os.path.join(runname,'sld_mol')
	direxist=os.path.exists(sldpath)
	if(direxist==0):
		os.system('mkdir -p '+sldpath)


	sld_output_files = []

	for i in range(number_of_frames):
        	nst = str(i+1).zfill(5) 
        	sld_output_files.append(sldpath+'sldfile_'+nst+'.txt')
	
	return sldpath,sld_output_files

def get_deuterated_atoms(num_deut_regions,deut_low_res,deut_high_res):

	deuterated_residues = []
	for i in xrange((num_deut_regions)):
		for k in xrange(deut_low_res[i],deut_high_res[i]+1):
			deuterated_residues.append(k)	
	return deuterated_residues

def update_gui(fraction_done,txtOutput):
	progress_string='\nCOMPLETED '+str(0+1)+' of '+str(1000)+' : '+str(fraction_done*100.0)+' % done'
	report_string='STATUS\t'+str(fraction_done)
	txtOutput.put(report_string)

	return

def sld_mol_main(variables,txtOutput):

	runname,path,pdbfile,dcdfile,expdatafile,outputfile,plotflag,runtype,bulk_sld,xon,num_deut_regions,deut_low_res,deut_high_res,dbin,width,sldfit,sldoffset,zfit0,zfitmin,zfitmax,zevalmin,zevalmax,A0,Amin,Amax = unpack_variables(variables)

	m1 = sasmol.SasMol(0)
	m1.read_pdb(path+pdbfile)

	print 'calculating minmax for all frames'
	
	minmax = m1.calc_minmax_all_steps(path+dcdfile)	

	minmaxz = [minmax[0][2],minmax[1][2]]

	binnum = int((minmax[1][2]-minmax[0][2])/dbin)+1

	tdcdinputfile = m1.open_dcd_read(path+dcdfile)
	print 'reading a single dcd frame'
	m1.read_dcd_step(tdcdinputfile,0)
	result = m1.close_dcd_read(tdcdinputfile[0])

	print 're-opening dcdfile'

	dcdinputfile = m1.open_dcd_read(path+dcdfile)
	
	number_of_frames = dcdinputfile[2]

	print 'number_of_frames = ',number_of_frames
	print 'dcdinputfile = ',dcdinputfile

	sldpath,sld_output_files = run_file_utilities(runname,number_of_frames)

	lineintxtOutput=''.join(['=' for x in xrange(60)])
	#ttxt=time.ctime()
	ttxt=time.asctime( time.gmtime( time.time() ) ) 
	txtOutput.put("\n%s \n" %(lineintxtOutput))
	txtOutput.put("DATA FROM RUN: %s \n\n" %(ttxt))

	graph = 0
	graph2 = 0

	if(plotflag == 2):
		plot_type = 'gnuplot'
		graph = Gnuplot.Gnuplot(debug=1)
		graph.clear()
		graph('set title "Initial SLD Calculation"')
		graph.xlabel('Z position (Angstroms)')
		graph.ylabel('Scattering Length Density')
	
		graph2 = Gnuplot.Gnuplot(debug=1)
		graph2.clear()
		graph2('set title "Final Fitted SLD Calculation"')
		graph2.xlabel('Z position (Angstroms)')
		graph2.ylabel('Scattering Length Density')

	iwidth = int(width/dbin)

	SLDh2o = -0.567E-6
	SLDd2o = 6.33E-6

	deut_ratio = (bulk_sld - SLDh2o) / (SLDd2o - SLDh2o)

	deuterated_residues = get_deuterated_atoms(num_deut_regions,deut_low_res,deut_high_res)

	fitparam0 = [zfit0,A0]
	parambounds = [[zfitmin, zfitmax], [Amin, Amax]]
	zeval = [zevalmin, zevalmax]
	
	zexp,interpolated_experimental_sld,interpolated_variance_data = read_experimental_sld_file(path,expdatafile,plotflag,graph,sldoffset)

	outfile=open(sldpath+outputfile,'w')
	st = "#  FRAME : OPTIMUM Z POSITION : OPTIMUM SURFACE COVERAGE : FIT ERROR\n"
	outfile.write(st)

	radius,volumes,residue_sld_hydrogen,residue_sld_deuterium,residue_sld_electron,residues = setup_sld_calculation(m1,xon,deuterated_residues)

	if(runtype == 0):
		
		txtOutput.put("PHASE 1:\n")
		txtOutput.put("\tCALCULATING SLD FOR ENSEMBLE\n")
		if(xon==0):
			volprof,slprofH,slprofD = calc_sld(m1,dcdinputfile,0,dbin,xon,runtype,minmaxz,radius,volumes,residue_sld_hydrogen,residue_sld_deuterium,residue_sld_electron,residues,txtOutput)
			slprof = (1.0-deut_ratio)*slprofH + deut_ratio*slprofD
		else:
			volprof,slprof,slprofD = calc_sld(m1,dcdinputfile,0,dbin,xon,runtype,minmaxz,radius,volumes,residue_sld_hydrogen,residue_sld_deuterium,residue_sld_electron,residues,txtOutput)
		txtOutput.put("\nPHASE 2:\n")
		if(sldfit == 1):
			txtOutput.put("\tSMOOTHING AND INTERPOLATING SLD FIT FOR ENSEMBLE\n")
			print 'SETTING UP FIT FOR SLD'
		else:
			txtOutput.put("\tSMOOTHING AND INTERPOLATING SLD FOR ENSEMBLE\n")
			print 'SETTING UP FOR SLD'

               	update_gui(0.0/5.0,txtOutput)
		
		print 'SMOOTING DATA 1 SLD'
               	update_gui(1.0/5.0,txtOutput)
		
		volprofsmooth = smooth_spectra(width,volprof,dbin)
		slprofsmooth = smooth_spectra(width,slprof,dbin)
		
		print 'SMOOTING DATA 2 SLD'
               	update_gui(2.0/5.0,txtOutput)
	
		sldprofsmooth = numpy.zeros(len(volprofsmooth),numpy.float)
		for ii in xrange(len(volprofsmooth)):
			if(volprofsmooth[ii] > 0):
				sldprofsmooth[ii] = slprofsmooth[ii]/volprofsmooth[ii]
	
		#sldprofsmooth = slprofsmooth/volprofsmooth
		vfprofsmooth = volprofsmooth/numpy.max(volprofsmooth)
		zprofsmooth = numpy.arange(len(sldprofsmooth))*dbin
		zmaxprofsmooth = len(sldprofsmooth)*dbin
		
		print 'INTERPOLATING DATA'
               	update_gui(3.0/5.0,txtOutput)

		fnsldprofsmooth = interpolate.interp1d(zprofsmooth,sldprofsmooth,'cubic')
		fnvfprofsmooth = interpolate.interp1d(zprofsmooth,vfprofsmooth,'cubic')

		if(sldfit == 1):
			print 'OPTIMIZING SLD FIT'
               	update_gui(4.0/5.0,txtOutput)

		if(sldfit == 1):	
			txtOutput.put("\tOPTIMIZING SLD FIT FOR ENSEMBLE\n")
			results = scipy.optimize.fmin_tnc(calc_error, fitparam0, fprime=None, args=(dbin, zmaxprofsmooth, zeval, fnvfprofsmooth, fnsldprofsmooth, bulk_sld, interpolated_experimental_sld, interpolated_variance_data, txtOutput), approx_grad=1, bounds=parambounds)
		
		else:
			results=[]
			results.append(fitparam0)

		fit_error = calc_error(results[0],dbin, zmaxprofsmooth, zeval, fnvfprofsmooth, fnsldprofsmooth,bulk_sld, interpolated_experimental_sld, interpolated_variance_data, txtOutput)
		
		print results[0][0], results[0][1], fit_error
		outfile.write('1\t'+`results[0][0]`+'\t'+`results[0][1]`+'\t'+`fit_error`+'\n')
		
		print 'PLOTTING SLD FIT'
               	update_gui(5.0/5.0,txtOutput)

		this_frame = 1

		new_dbin = float(len(volprof))*dbin/float(len(volprofsmooth))
	
		sum_sld_fit = numpy.zeros(len(vfprofsmooth),numpy.float)
		plot_and_save_sld_fit(results[0], new_dbin, zexp, interpolated_experimental_sld, vfprofsmooth, sldprofsmooth, bulk_sld, plotflag, graph2, sum_sld_fit, this_frame, number_of_frames, runtype, fit_error, results, sld_output_files, sldpath)

	elif(runtype == 1):
	
		sum_flag = 0
		for i in xrange(number_of_frames):
			
			if(xon==0):
				volprof,slprofH,slprofD = calc_sld(m1,dcdinputfile,0,dbin,xon,runtype,minmaxz,radius,volumes,residue_sld_hydrogen,residue_sld_deuterium,residue_sld_electron,residues,txtOutput)
				slprof = (1.0-deut_ratio)*slprofH + deut_ratio*slprofD
			else:
				volprof,slprof,slprofD = calc_sld(m1,dcdinputfile,0,dbin,xon,runtype,minmaxz,radius,volumes,residue_sld_hydrogen,residue_sld_deuterium,residue_sld_electron,residues,txtOutput)
		
			volprofsmooth = smooth_spectra(width,volprof,dbin)
			slprofsmooth = smooth_spectra(width,slprof,dbin)
		
			sldprofsmooth = numpy.zeros(len(volprofsmooth),numpy.float)
			for ii in xrange(len(volprofsmooth)):
				if(volprofsmooth[ii] > 0):
					sldprofsmooth[ii] = slprofsmooth[ii]/volprofsmooth[ii]
	
			#sldprofsmooth = slprofsmooth/volprofsmooth
			vfprofsmooth = volprofsmooth/numpy.max(volprofsmooth)
			zprofsmooth = numpy.arange(len(sldprofsmooth))*dbin
			zmaxprofsmooth = len(sldprofsmooth)*dbin
		
			fnsldprofsmooth = interpolate.interp1d(zprofsmooth,sldprofsmooth,'cubic')
			fnvfprofsmooth = interpolate.interp1d(zprofsmooth,vfprofsmooth,'cubic')
		
			if(sldfit == 1):	
				results = scipy.optimize.fmin_tnc(calc_error, fitparam0, fprime=None, args=(dbin, zmaxprofsmooth, zeval, fnvfprofsmooth, fnsldprofsmooth, bulk_sld, interpolated_experimental_sld, interpolated_variance_data, txtOutput), approx_grad=1, bounds=parambounds)
			else:
				results=[]
				results.append(fitparam0)
		
			fit_error = calc_error(results[0],dbin, zmaxprofsmooth, zeval, fnvfprofsmooth, fnsldprofsmooth,bulk_sld, interpolated_experimental_sld, interpolated_variance_data, txtOutput)
			
			this_frame = i+1
			print results[0][0], results[0][1], fit_error
			outfile.write(`this_frame`+'\t'+`results[0][0]`+'\t'+`results[0][1]`+'\t'+`fit_error`+'\n')
	
			if(sum_flag == 0):
				sum_sld_fit = numpy.zeros(len(vfprofsmooth),numpy.float)
				sum_flag = 1
		
			new_dbin = float(len(volprof))*dbin/float(len(volprofsmooth))
		
			plot_and_save_sld_fit(results[0], new_dbin, zexp, interpolated_experimental_sld, vfprofsmooth, sldprofsmooth, bulk_sld, plotflag, graph2, sum_sld_fit, this_frame, number_of_frames, runtype, fit_error, results, sld_output_files, sldpath)
	
                	if(((i+1)%(float(number_of_frames)/100.0)==0 or (i<10))):
				print_status(i,number_of_frames,txtOutput)
	
	outfile.close()   

	txtOutput.put("\nscattering length profile(s) and statistics are saved in %s directory\n\n" % ('./'+sldpath))

	txtOutput.put("\n%s \n" %(lineintxtOutput))
	time.sleep(2)
	print 'SLD IS DONE'

	return

if __name__=='__main__':

	### start user input 
	### start user input 
	### start user input 

	runname='run_0'
	path='./'
	pdbfile = 'nef_nohis.pdb'
	dcdfile = 'aligned_rnz2.pdb_run_allrg_nohis.dcd'
	dcdfile = 'a5.dcd'
	expdatafile = 'SNS_dDPPG_myrdNef_nSLD.txt'
	outputfile = 'tmp_smoothdist.txt'
	plotflag = '2'    # 0 == NO PLOT, 1 == matplotlib, 2 == Gnuplot

	runtype = '1'	# 0 === average sld over all structures, 1 == best fit sld for each individual structure
	bulk_sld = '-0.51e-6'
		
	dbin = '0.2'		# heavy atom profile

	num_deut_regions = '1'
	deut_low_res = '1'
	deut_high_res = '200'

	### deuterated_residues = numpy.arange(10,211)
	
	width = '2.5'		# Gaussian smooting

	sldfit = '1'		# 0 == no fit, 1 == fit
	sldoffset = '0.0'		# offset to experimental sld

	zfit0 ='24.0'			# SLD fitting
	zfitmin ='20.0'
	zfitmax ='28.0'
	A0 = '0.30'
	Amin = '0.250'
	Amax = '0.35'

	zevalmin = '28.0'
	zevalmax = '145.0' #z-range to evaluate the error function
	
	xon = '0'

	### end user input 
	### end user input 
	### end user input 

	svariables={}

	svariables['runname']           = (runname,'string')
	svariables['path']              = (path,'string')
	svariables['pdbfile']           = (pdbfile,'string')
	svariables['dcdfile']           = (dcdfile,'string')
	svariables['expdatafile']	= (expdatafile,'string')
	svariables['outputfile']	= (outputfile,'string')

	svariables['runtype']		= (runtype,'int')
	svariables['bulk_sld']		= (bulk_sld,'float')

	svariables['xon']		= (xon,'int')
	svariables['num_deut_regions']	= (num_deut_regions,'int')
	svariables['deut_low_res']	= (deut_low_res,'int_array')
	svariables['deut_high_res']	= (deut_high_res,'int_array')

	svariables['dbin']		= (dbin,'float')
	svariables['width']		= (width,'float')

	svariables['sldfit']		= (sldfit,'int')
	svariables['sldoffset']		= (sldoffset,'float')
	
	svariables['zfit0']		= (zfit0,'float')
	svariables['zfitmin']		= (zfitmin,'float')
	svariables['zfitmax']		= (zfitmax,'float')
	svariables['zevalmin']		= (zevalmin,'float')
	svariables['zevalmax']		= (zevalmax,'float')

	svariables['A0']		= (A0,'float')
	svariables['Amin']		= (Amin,'float')
	svariables['Amax']		= (Amax,'float')

	svariables['plotflag']		= (plotflag,'int')
	
	import sassie.interface.input_filter as input_filter
	import sassie.interface.sld_filter as sld_filter

	error,variables=input_filter.type_check_and_convert(svariables)

	if(len(error)>0):
		print 'error = ',error
		sys.exit()	

	error=sld_filter.check_sld(variables)

	if(len(error)>0):
		print 'error = ',error
		sys.exit()	

	txtQueue=multiprocessing.JoinableQueue()
	sld_main(variables,txtQueue)
	#process=multiprocessing.Process(target=sld_main,args=(variables,txtQueue))
	#process.start()

