'''
    SASSIE  Copyright (C) 2011 Joseph E. Curtis
    This program comes with ABSOLUTELY NO WARRANTY; 
    This is free software, and you are welcome to redistribute it under certain
    conditions; see http://www.gnu.org/licenses/gpl-3.0.html for details.
''' 
import os,sys,string,platform,time
import getpass,readline

#       UPGRADE 
#
#       09/05/2011      --      initial coding 			:    	jc
#       06/27/2012      --      courtesy installer		:	jc
#       10/27/2012      --      altered for simple updates	:	jc
#
#LC      1         2         3         4         5         6         7
#LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''
	upgrade.py is the script to install a newer version of sassie
	
	this version supports MacOS X:
		 10.5.8 -- Leopard 
		 10.6.8 -- Snow Leopard
		 10.7.4 -- Lion
		 10.8   -- Mountain Lion (experimental)

	Ubuntu (10.04 LTS & 11.04) and CentOS (5.6, 5.8 & 6.0).   

'''

def log(logfile,line):

	print line

	logfile.write("%s\n" % (line))
	
	return	

def print_error(logfile,error):

	log(logfile,'\n>>> ERROR detected : '+error)
	
	return

def determine_os(logfile,install_path,current_path):

	supported_flag = 'NO'

	log(logfile,'platform.system() = '+platform.system())
	log(logfile,'platform.platform() = '+platform.platform())
	log(logfile,'platform.processor() = '+platform.processor())
	log(logfile,'platform.architecture()[0] = '+platform.architecture()[0])
	log(logfile,'platform.architecture()[1] = '+platform.architecture()[1])

	bit = platform.architecture()[0]

	local_os_type = platform.system()

	if(local_os_type == "Darwin"):
		log(logfile,'Apple OS Detected')
		(release,versioninfo,machine) = platform.mac_ver()
		log(logfile,'release = '+release)
		log(logfile,'machine = '+machine)

		os_type = [local_os_type,release,machine,platform.architecture()[0]]
	
		if(release == "10.5.8"):
			supported_flag = 'YES'
			python = 'python'
			cpst = 'cp -v '+current_path+'/sassie/sasmol/extensions/matrix_math/matrix_math.py '+current_path+'/sassie/sasmol/extensions/matrix_math/old_matrix_math.py'	
        		result = os.popen(cpst).readlines()
        		for line in result: log(logfile,line)

			cpst = 'cp -f -v '+current_path+'/sassie/sasmol/extensions/matrix_math/matrix_math_osx_10.5.8.py '+current_path+'/sassie/sasmol/extensions/matrix_math/matrix_math.py'	
        		result = os.popen(cpst).readlines()
        		for line in result: log(logfile,line)

		elif(release == "10.6.8"):
			supported_flag = 'YES'
			python = 'python'
		elif(release[:4] == "10.7"):
			supported_flag = 'YES'
			python = 'python'
		elif(release[:4] == "10.8"):
			supported_flag = 'YES'
			python = 'python'
		else:
			error = 'this installer does not handle your OS'
			error += ' update your leopard (10.5.X to 10.5.8)\n'
			error += 'or snow leopard (10.6.X to 10.6.8)\n'
			print_error(logfile,error)
	
	elif(local_os_type == "Linux"):

		log(logfile,'Linux OS Detected')

		if(bit == '64bit'):
			error = '\n\nthis installer does not handle Linux 64-bit machines'
			error +='\n\nINSTALLATION STOPPING NOW\n\n'
			print_error(logfile,error)

		try:
			(distname,version,id) = platform.linux_distribution()
		except:
			(distname,version,id) = platform.dist()
			
		log(logfile,'distname = '+distname)
		log(logfile,'version = '+version)
		log(logfile,'id = '+id)

		try:
			(lib,gversion) = platform.libc_ver()
			log(logfile,'glib  = '+lib+' version = '+gversion)
		except:
			log(logfile,'could not determine glib used to compile current python')

		if((distname == "redhat" or distname == 'CentOS') and version == "5.6"):
			supported_flag = 'YES'
			python = install_path + '/bin/python2.7'
		elif((distname == "redhat" or distname == 'CentOS Linux') and version == "6.0"):
			supported_flag = 'YES'
			python = install_path + '/bin/python2.7'
		elif(distname == "Ubuntu" and version == "10.04"):
			supported_flag = 'YES'
			python = 'python'
		elif(distname == "Ubuntu" and version == "11.04"):
			supported_flag = 'YES'
			python = 'python'
		else:
			error = 'this installer does not handle your OS'
			log(logfile,error)
			print_error(logfile,error)

		os_type = [local_os_type,distname,version,platform.architecture()[0]]

	if(supported_flag == 'NO'):
		log(logfile,'\n >> your OS is not currently supported: but the installer may work\n')
		log(logfile,' >> if you have a "flavor" of either RedHat or Debian\n')

		try_anyway = ""

		st0 = '\n \t\t (0) QUIT NOW \n'
		st1 = '\n \t\t (1) REDHAT TYPE (CentOS, RHEL, Fedora,etc.)\n'
		st2 = '\n \t\t (2) DEBIAN TYPE (Debian, Ubuntu, Linux Mint, etc.)\n'
		st3 = '\n >> quit (0) or enter OS type do you want to mimic ? : (1) or (2) '
		log(logfile,st0) ; log(logfile,st1) ; log(logfile,st2)

		while(try_anyway != "0" and try_anyway != "1" and try_anyway != "2"):
			try_anyway = raw_input(st3)	

		st4 = '\n >> you entered '+try_anyway
	
		log(logfile,st4)

		if(try_anyway == "0"):
			log(logfile,'\n\n QUITING NOW \n\n')
			sys.exit()
		elif(try_anyway == "1"):
			local_os_type = "Linux"
			distname = "redhat"
			os_type = [local_os_type,distname,version,platform.architecture()[0]]
			supported_flag = 'MAYBE'
			python = install_path + '/bin/python2.7'
		elif(try_anyway == "2"):
			local_os_type = "Linux"
			distname = "Ubuntu"
			os_type = [local_os_type,distname,version,platform.architecture()[0]]
			supported_flag = 'MAYBE'
			python = 'python'

	log(logfile,'>OS_TYPE = '+os_type[0])
	log(logfile,'>IS OS SUPPORTED? = '+supported_flag)
	
	return os_type,python

def compile_extensions(logfile,current_path,python):

	spath = current_path+'/sassie/sasmol/'
	path = current_path+'/sassie/sasmol/extensions/'
	dpath = path+'dcdio/'
	mpath = path+'mask/'

	os.chdir(dpath)
	buildst1 = python +' setup_dcdio.py build'
	result = os.popen(buildst1).readlines() 
	for line in result: log(logfile,line)
	
	os.chdir(mpath)
	buildst2 = python +' setup_mask.py build'
	result = os.popen(buildst2).readlines() 
	for line in result: log(logfile,line)

	os.chdir(current_path)	
	
	cpst1 = 'cp '+dpath+'build/lib*/_dcdio.so '+spath
	cpst2 = 'cp '+dpath+'dcdio.py '+spath
	cpst3 = 'cp '+mpath+'build/lib*/_mask.so '+spath
	cpst4 = 'cp '+mpath+'mask.py '+spath

	result = os.popen(cpst1).readlines() 
	for line in result: log(logfile,line)
	result = os.popen(cpst2).readlines() 
	for line in result: log(logfile,line)
	result = os.popen(cpst3).readlines() 
	for line in result: log(logfile,line)
	result = os.popen(cpst4).readlines() 
	for line in result: log(logfile,line)
	
	if(not os.path.isfile(spath+'_dcdio.so')):
		print_error(logfile,"_dcdio.so did not install correctly!\n\nINSTALLATION STOPPING NOW\n\n")
		sys.exit()
	if(not os.path.isfile(spath+'dcdio.py')):
		print_error(logfile,"dcdio.py did not install correctly!\n\nINSTALLATION STOPPING NOW\n\n")
		sys.exit()
	if(not os.path.isfile(spath+'_mask.so')):
		print_error(logfile,"_mask.so did not install correctly!\n\nINSTALLATION STOPPING NOW\n\n")
		sys.exit()
	if(not os.path.isfile(spath+'mask.py')):
		print_error(logfile,"mask.py did not install correctly!\n\nINSTALLATION STOPPING NOW\n\n")
		sys.exit()

	return

def compile_sassie(logfile,current_path,python):

	buildst = python+' setup.py build'

	installst = python+' setup.py install'
	
	result = os.popen(buildst).readlines() 
	for line in result: log(logfile,line)
	result = os.popen(installst).readlines() 
	for line in result: log(logfile,line)

	return

def install_sassie(logfile,os_type,current_path,install_path,python):

	compile_extensions(logfile,current_path,python)
	compile_sassie(logfile,current_path,python)

	return

def module_installation(logfile,current_path,install_path,module_list):

	sasgui_file = current_path+'/testing/sassie_gui.py'
	sasgui = open(sasgui_file,'r').readlines()
	sasgui2 = open(sasgui_file,'w')

	aprogs = ['Contrast Calculator',
		'Align', 
		'Data Interpolation', 
		'Coordinate Tools',
		'Merge Utilities',
		'Monomer Monte Carlo', 
		'Complex Monte Carlo', 
		'Energy Minimization', 
		'Docking',
		'SasCalc', 
		'SLDMOL', 
		'EM to SAS', 
		'Chi-Square Filter', 
		'Density Plot',
		'APBS']


	log(logfile,'\n >> writing sassie_gui.py script for installed modules')
	check = 0	
	all_dict = {}
	for i in xrange(len(sasgui)):
		lin = string.split(sasgui[i])
		if(len(lin)>0 and lin[0] == '}' and check == 1):
			for prog in aprogs:
				if prog in module_list:
					tline = '\t\t"'+prog+'"\t:\t"YES",'
				else:
					tline = '\t\t"'+prog+'"\t:\t"NO",'
				log(logfile,'altering sassie_gui.py file: '+tline)
				if prog == aprogs[-1]:	
					sasgui2.write(tline[:-1]+'\n')
				else:
					sasgui2.write(tline+'\n')
		#	sasgui2.write('\t\t\t\t}\n')
			check = 0		
			
		elif(len(lin)>0):
			if(lin[0] == 'desired_modules'): 
				check = 1	
				sasgui2.write(sasgui[i])
		if(check == 0):			
			sasgui2.write(sasgui[i])

	sasgui2.close()
	
	return

def clean_up(os_type):

	return

def check_permissions(path):

	try:
		existvalue = os.access(path,os.F_OK)
	except:
		existvalue = False
 	try:
		readvalue = os.access(path,os.R_OK)
	except:
		readvalue = False
	try:
		writevalue = os.access(path,os.W_OK)
	except:
		writevalue = False

	return existvalue,readvalue,writevalue

def check_path(path):

	error=[]
	ev,rv,wv=check_permissions(path)
	if(not ev or not rv or not wv):
		error.append('permission error in input file path '+path+'  [PATH = '+str(ev)+': READ PERMISSION = '+str(rv)+': WRITE PERMISSION = '+str(wv)+']')
		if(ev==False):
			error.append('path does not exist')
		elif(rv==False):
			error.append('read permission not allowed')
		elif(wv==False):
			error.append('write permission not allowed')
		return error

	return error

def test_installation(logfile,os_type):

	log(logfile,'\nTESTING INSTALLATION\n')

	return

def preamble(logfile):

	log(logfile,'\nSASSIE UPDATER : version 0.99 : 10/27/2012\n')

	user = getpass.getuser()
	current_path = os.getcwd()
	
	log(logfile,'executed by user : '+user+' on '+time.ctime())
	log(logfile,'current directory : '+current_path)
	log(logfile,'hostname : '+platform.node()+'\n')

	#s if(user != "root"):
	#s 	print_error(logfile,'\n\nYOU MUST RUN THIS SCRIPT AS ROOT\n\n')
	#s 	log(logfile,'\n\nINSTALLATION STOPPING NOW\n\n')
	#s 	sys.exit()

	log(logfile,'\n >> This script will compile and install software')
	log(logfile,'\n >> NOTE: USE MUST USE THE SAME VERSION OF PYTHON YOU') 
	log(logfile,'\n >> HAVE USED TO RUN SASSIE BEFORE')
	log(logfile,'\n >> NOTE: this script assumes you have a previous version')
	log(logfile,'\n >> of SASSIE installed and you are running this script')
	log(logfile,'\n >> to update only.  If this is not the case then run')
	log(logfile,'\n >> the installer.py script as root')

	st7 = '\nDo you wish to continue? : (0)==no or (1)==yes '

	continue_installation = False
	while(continue_installation != "0" and continue_installation != "1"):
		continue_installation = raw_input(st7)

	if(continue_installation=="0"):
		log(logfile,'\n\n STOPPING HERE\n\n')
		sys.exit()

	print '\n'

	install_path = '/usr/local/bin'

	return current_path,install_path

def epilogue(logfile,os_type):

	log(logfile,'\n\nINSTALLATION COMPLETE \n\n')
	log(logfile,'\n >> you can start sassie by running the script')
	log(logfile,' >> "sassie_gui.py" which is in a directory called')
	log(logfile,' >> "testing" in this directory')
	log(logfile,'\n >> you should copy this file to your own working directory')
	log(logfile,' >> to work on your own projects')

	return

def set_modules(os_type, **kwargs):
 
	full_module_list = ['Contrast_Calculator','Align','Data Interpolation','Merge Utilities','Coordinate Tools','Monomer Monte Carlo','Complex Monte Carlo','Energy Minimization','Docking','SasCalc','SLDMOL','EM to SAS','Chi-Square Filter','Density Plot','APBS']

	excluded_modules = []

	module_list = []

	for module in full_module_list:
		if module not in excluded_modules:
			module_list.append(module)
	
	return module_list

def remove_old_sassie(**kwargs):

	save_log = False

	if 'log_file' in kwargs:
		save_log = True
		logfile = kwargs['log_file']

	from distutils.sysconfig import get_python_lib
	lib_path = get_python_lib()
	if save_log: 
		log(logfile,lib_path)
	else:
		print lib_path

	rmst = 'rm -Rf '+lib_path+'/sassie*'
	print 'rmst = ',rmst
	result = os.popen(rmst).readlines()
	for line in result: 
		if save_log: 
			log(logfile,line)
		else:
			print line

	rmst = 'rm -Rf ./build/'
	print 'rmst = ',rmst
	result = os.popen(rmst).readlines()
	for line in result: 
		if save_log: 
			log(logfile,line)
		else:
			print line
	
	return

def install_me():

	logfile = open('log_sassie_upgrade.txt','w')

	current_path,install_path = preamble(logfile)
	
	os_type,python = determine_os(logfile,install_path,current_path)

	print 'os_type = ',os_type

	remove_old_sassie(log_file=logfile)

	module_list = set_modules(os_type)

	module_installation(logfile,current_path,install_path,module_list)

	install_sassie(logfile,os_type,current_path,install_path,python)

#	test_installation(logfile,os_type)

	epilogue(logfile,os_type)

	logfile.close()

	return


if __name__=='__main__':
      install_me()

