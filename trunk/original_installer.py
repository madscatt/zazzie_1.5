'''
    SASSIE  Copyright (C) 2011-2015 Joseph E. Curtis
    This program comes with ABSOLUTELY NO WARRANTY;
    This is free software, and you are welcome to redistribute it under certain
    conditions; see http://www.gnu.org/licenses/gpl-3.0.html for details.
'''
import os,sys,string,platform,time
import getpass, readline

#       INSTALLER
#
#       09/05/2011      --      initial coding 			    :   jc
#       06/27/2012      --      courtesy installer		:	jc
#       12/12/2013      --      added sasview/tamd 		:	jc
#       02/19/2014      --      osx mavericks 10.9 		:	jc
#       11/22/2014      --      updated for 2.0 branch  :               jc
#
#LC      1         2         3         4         5         6         7
#LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''
	installer.py is the script to install sassie and required dependencies
	
	this version supports MacOS X:
		 10.5.8 -- Leopard 
		 10.6.8 -- Snow Leopard
		 10.7.4 -- Lion
		 10.8   -- Mountain Lion (experimental)
		 10.9   -- Mavericks (experimental)

	Ubuntu (10.04 LTS & 11.04) and CentOS (5.6, 5.8 & 6.0).   

	See the README file in this directory for information on how
	to use this installer script.
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
			cpst = 'cp -v '+current_path+'/sassie/sasmol/extensions/matrix_math/matrix_math.f '+current_path+'/sassie/sasmol/extensions/matrix_math/old_matrix_math.f'	
        		result = os.popen(cpst).readlines()
        		for line in result: log(logfile,line)

			cpst = 'cp -f -v '+current_path+'/sassie/sasmol/extensions/matrix_math/matrix_math_osx_10.5.8.f '+current_path+'/sassie/sasmol/extensions/matrix_math/matrix_math.f'	
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
		elif(release[:4] == "10.9"):
			print '>>> Mavericks installation'
			supported_flag = 'YES'
			python = 'python'

			xcodest = 'xcode-select --install'
        		result = os.popen(xcodest).readlines()
        		for line in result: log(logfile,line)
        		st3 = '>>>> press a key AFTER Xcode command line tools are installed <<<<'
        		dum_input = raw_input(st3)
		else:
			error = 'this installer does not handle your OS'
			error += ' update your leopard (10.5.X to 10.5.8)\n'
			error += 'or snow leopard (10.6.X to 10.6.8)\n'
			error += 'or to lion (10.7), Mountain Lion (10.8)\n'
			error += 'or to Mavericks (10.9)\n'
			print_error(logfile,error)
	
	elif(local_os_type == "Linux"):

		log(logfile,'Linux OS Detected')

		if(bit == '64bit'):
			error = '\n\nthis installer does not robustly handle Linux 64-bit machines'
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
		elif(distname == "Ubuntu" and version == "14.04"):
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

def find_files(logfile,files_to_check):

	programs_needed = []

	for i in xrange(len(files_to_check)):
		locatefile = os.popen('which '+files_to_check[i]).readlines()
		if(locatefile):
			log(logfile,'found '+locatefile[0])
		else:
			log(logfile,'file '+files_to_check[i]+' not found!')
			programs_needed.append(files_to_check[i])

	if(len(programs_needed) > 0):	
		log(logfile,'\n >> dependencies summary: programs needed = '+', '.join(programs_needed))
	else:
		log(logfile,'\n >> dependencies summary:\n >> it seems that you have the programs needed to continue the installation\n')

	return programs_needed

def apt_install(logfile,program):

	log(logfile,'>>installing '+program+' using apt')
	error=[]
	st1 = 'apt-get -y install '+program
	try:
		result = os.popen(st1).readlines()
		for line in result:
			log(logfile,line)
	except:
		error.append = '>>> failed to install '+program+'\nINSTALLATION STOPPING\n\n'
		print_error(error)
		sys.exit()

	return error		

def yum_install(logfile,program):

	log(logfile,'>>installing '+program+' using yum')
	error=[]
	st1 = 'yum -y install '+program
	try:
		result = os.popen(st1).readlines()
		for line in result:
			log(logfile,line)
	except:
		error.append = '>>> failed to install '+program+'\nINSTALLATION STOPPING\n\n'
		print_error(error)
		sys.exit()

	return error		

def install_gfortran(logfile,current_path,install_path):

     	(release,versioninfo,machine) = platform.mac_ver()  

	log(logfile,'> installing gfortran')
        locpath = current_path+'/sassie/distribution/macosx/'

	if(release[:4] != "10.9"):
		untarst = 'tar -xvf '+locpath+'gfortran-intel-leopard.tar -C /'
	elif(release[:4] == "10.9"):
		untarst = 'tar -xvf '+locpath+'gfortran-4.9-bin_mavericks.tar -C /'
        result = os.popen(untarst).readlines()
        for line in result: log(logfile,line)
	
	return 
	
def install_gnuplot(logfile,current_path,install_path):

        log(logfile,'> installing gnuplot')
        locpath = current_path+'/sassie/distribution/macosx/'
        untarst = 'tar -xvf '+locpath+'gnuplot-4.4.3.tar.gz -C '+locpath
        result = os.popen(untarst).readlines()
        for line in result: log(logfile,line)

        os.chdir(locpath+'gnuplot-4.4.3')

        configst = './configure --prefix=/usr/local --with-readline=builtin'
        result = os.popen(configst).readlines()
        for line in result: log(logfile,line)

        makest = 'make'
        result = os.popen(makest).readlines()
        for line in result: log(logfile,line)

        instst = 'make install'
        result = os.popen(instst).readlines()
        for line in result: log(logfile,line)

        os.chdir(current_path)

        if(not os.path.isfile('/usr/local/bin/gnuplot')):
                print_error(logfile,"gnuplot was not installed correctly!\n\nINSTALLATION STOPPING NOW\n\n")
                sys.exit()

	return
	
def install_pcre(logfile,current_path,install_path):
	
	locpath = current_path+'/sassie/distribution/general/'
	untarst = 'tar -xvjf '+locpath+'pcre-8.20.tar.bz2 -C '+locpath
	result = os.popen(untarst).readlines()
	for line in result: log(logfile,line)

	os.chdir(locpath+'pcre-8.20')

	configst = './configure --prefix=/usr/local/'
	result = os.popen(configst).readlines()
	for line in result: log(logfile,line)

	makest = 'make'
	result = os.popen(makest).readlines()
	for line in result: log(logfile,line)

	makest = 'make install'
	result = os.popen(makest).readlines()
	for line in result: log(logfile,line)

	os.chdir(current_path)

def install_swig(logfile,current_path,install_path):

	install_pcre(logfile,current_path,install_path)

	locpath = current_path+'/sassie/distribution/general/'
	untarst = 'tar -xvzf '+locpath+'swig-2.0.4.tar.gz -C '+locpath
	result = os.popen(untarst).readlines()
	for line in result: log(logfile,line)
	
	os.chdir(locpath+'swig-2.0.4')
	
	configst = './configure --prefix=/usr/local/'
	result = os.popen(configst).readlines()
	for line in result: log(logfile,line)

	makest = 'make'
	result = os.popen(makest).readlines()
	for line in result: log(logfile,line)

	makest = 'make install'
	result = os.popen(makest).readlines()
	for line in result: log(logfile,line)

	os.chdir(current_path)

	return

def install_gnuplot_py(logfile,current_path):

	log(logfile,'> installing gnuplot-py')
	locpath = current_path+'/sassie/distribution/general/'
	
	untarst = 'tar -xvf '+locpath+'gnuplot-py-1.8.tar -C '+locpath
	result = os.popen(untarst).readlines()
	for line in result: log(logfile,line)

	os.chdir(locpath+'gnuplot-py-1.8')
	
	buildst = 'python setup.py build'
	result = os.popen(buildst).readlines()
	for line in result: log(logfile,line)
		
	instst = 'python setup.py install'
	result = os.popen(instst).readlines()
	for line in result: log(logfile,line)
	
	os.chdir(current_path)
	
	return

def darwin_dependencies(logfile,current_path,install_path):

	os_type = "Darwin"
	
	log(logfile,' >> checking Mac OS installation\n')
	log(logfile,' >> NOTE: the installer assumes that you have installed')
	log(logfile,' >> the developer applications in their normal locations\n')	

	xcode="/Developer/Applications/Xcode.app/Contents/MacOS/Xcode"
	xcode_LION_NEW="/Applications/Xcode.app/Contents/MacOS/Xcode"
	xcode_ML="/Applications/Xcode44-DP5.app/Contents/MacOS/Xcode"
	xcode_Mavericks="/Applications/Xcode5-DP.app/Contents/MacOS/Xcode"

	if(not os.path.isfile(xcode) and not os.path.isfile(xcode_Mavericks) and not os.path.isfile(xcode_ML) and not os.path.isfile(xcode_LION_NEW)):
		error = '\n\nXcode is required to complete the installation\n'
		error += 'You can find this on your original Mac OS Install disk\n'
		error += 'Make sure to run "software update" as well\n'
		print_error(logfile,error)
		log(logfile,'XCODE WAS NOT FOUND: QUITING NOW\n\n')
		sys.exit()
	else:
		log(logfile,'found '+xcode+'\n')

	epd_python = "/Library/Frameworks/Python.framework/Versions/Current/bin/python"
	#epd_python = "/Library/Frameworks/EPD64.framework/Versions/Current/bin/python"

	if(not os.path.isfile(epd_python)):
		error = '\n\nEnthought "Free" Python is required to complete the installation\n'
		error += 'You can find the installer in '+current_path+'/sassie/distribution/macosx\n'
		print_error(logfile,error)
		log(logfile,'EPD Python WAS NOT FOUND: QUITING NOW\n\n')
		sys.exit()

	aquaterm = '/Applications/AquaTerm.app/Contents/MacOS/AquaTerm'

	if(not os.path.isfile(aquaterm)):
		error = '\n\nAquaTerm is required to complete the installation\n'
		error += 'You can find the installer in '+current_path+'/sassie/distribution/macosx\n'
		print_error(logfile,error)
		log(logfile,'AquaTerm WAS NOT FOUND: QUITING NOW\n\n')
		sys.exit()

	files_to_check = ['gcc','g++','gfortran','gnuplot','swig']

	programs_needed = find_files(logfile,files_to_check)

	if('gcc' in programs_needed or 'g++' in programs_needed):
		error = '\n >> Apple Developer Tools should be installed\n'
		error += '>> for some reason Xcode is installed yet gcc or g++ could not be found\n'
		error += '>> YOU MAY NEED TO INSTALL COMMAND LINE TOOLS\n'
		error += '>> YOU MAY NEED TO INSTALL COMMAND LINE TOOLS\n'
		error += '>> YOU MAY NEED TO INSTALL COMMAND LINE TOOLS\n'
		log(logfile,'\n\n >> INSTALLATION STOPPED: NO GCC/G++ FOUND: QUITING NOW\n\n')
 		sys.exit()

	if('gfortran' in programs_needed):
		install_gfortran(logfile,current_path,install_path)

	if('gnuplot' in programs_needed):
		install_gnuplot(logfile,current_path,install_path)
	
	if('swig' in programs_needed):
		error = '\n >> swig needs to be installed\n'
		#log(logfile,'\n\n >> INSTALLATION STOPPED: NO swig FOUND: QUITING NOW\n\n')
 		#sys.exit()
		install_swig(logfile,current_path,install_path)

	install_gnuplot_py(logfile,current_path)

	install_mdx(logfile,current_path,install_path,os_type)
	
	return


def install_python(logfile,current_path,install_path):

	locpath = current_path+'/sassie/distribution/general/'
	untarst = 'tar -xvf '+locpath+'Python-2.7.2.tar -C '+locpath
	
	result = os.popen(untarst).readlines()
	for line in result: log(logfile,line)
	
	os.chdir(locpath+'Python-2.7.2')
		
	configst = './configure --prefix='+install_path
	result = os.popen(configst).readlines()
	for line in result: log(logfile,line)
	
	makest = 'make'
	result = os.popen(makest).readlines()
	for line in result: log(logfile,line)
	
	instst = 'make altinstall'
	result = os.popen(instst).readlines()
	for line in result: log(logfile,line)
	
	untarst = 'tar -xvf '+locpath+'numpy-1.6.1.tar -C '+locpath
	result = os.popen(untarst).readlines()
	for line in result: log(logfile,line)
		
	os.chdir(locpath+'numpy-1.6.1')
		
	buildst = install_path+'/bin/python2.7 setup.py build'
	result = os.popen(buildst).readlines()
	for line in result: log(logfile,line)
	
	instst = install_path+'/bin/python2.7 setup.py install'
	result = os.popen(instst).readlines()
	for line in result: log(logfile,line)
	
	untarst = 'tar -xvf '+locpath+'scipy-0.9.0.tar -C '+locpath
	result = os.popen(untarst).readlines()
	for line in result: log(logfile,line)
	
	os.chdir(locpath+'scipy-0.9.0')
		
	result = os.popen(buildst).readlines()
	for line in result: log(logfile,line)
	
	result = os.popen(instst).readlines()
	for line in result: log(logfile,line)
		
	untarst = 'tar -xvf '+locpath+'gnuplot-py-1.8.tar -C '+locpath
	result = os.popen(untarst).readlines()
	for line in result: log(logfile,line)

	os.chdir(locpath+'gnuplot-py-1.8')
	
	result = os.popen(buildst).readlines()
	for line in result: log(logfile,line)
		
	result = os.popen(instst).readlines()
	for line in result: log(logfile,line)
	
	os.chdir(current_path)

	return

def redhat_dependencies(logfile,current_path,install_path):

	os_type = "redhat"
	
	log(logfile,' >> checking RedHat installation\n')
	files_to_check = ['gcc','g++','gfortran','gnuplot','tcsh','swig']
	programs_needed = find_files(logfile,files_to_check)

	files_to_append = ['tcl','tcl-devel','tk','tk-devel','compat-libstdc++-33','readline-devel','ncurses-devel','subversion','blas','blas-devel','lapack','lapack-devel','zlib','zlib-devel']
	
	for file in files_to_append:
		programs_needed.append(file)
	
	log(logfile,'>> updating system via yum')

	result = os.popen('yum -y update').readlines()
	for line in result: log(logfile,line)

	try:
		error = yum_install(logfile,'gcc-c++')
	except:
		pass

	for program in programs_needed:
		if(program == 'gcc'):
			error = yum_install(logfile,'gcc')
			error = yum_install(logfile,'compat-gcc-34')
			error = yum_install(logfile,'compat-gcc-34-c++')
		if(program == 'g++'):
			error = yum_install(logfile,'gcc-c++')
		if(program == 'gfortran'):
			error = yum_install(logfile,'gcc-gfortran')
		else:	
			error = yum_install(logfile,program)

	result = os.popen('touch '+current_path+'/README')
	for line in result: log(logfile,line)

	install_python(logfile,current_path,install_path)

	return

def ubuntu_dependencies(logfile,current_path,install_path,version):
	
	os_type = "Ubuntu"
	
	log(logfile,' >> checking Ubuntu installation\n')
        if(version == '14.04'):
	    files_to_check = ['gcc','g++','gfortran','gnuplot-X11','tcsh','swig']
	
        else:
            files_to_check = ['gcc','g++','gfortran','gnuplot','tcsh','swig']

	programs_needed = find_files(logfile,files_to_check)

	files_to_append = ['tcl8.4','tcl8.4-dev','tk8.4','tk8.4-dev','libfftw3-3','libfftw3-dev','libstdc++5','python-numpy','python-scipy','python-tk','python-all-dev','python-gnuplot']
	for file in files_to_append:
		programs_needed.append(file)
	
	result = os.popen('cp /etc/apt/sources.list /etc/apt/sources.list.backup').readlines()
	for line in result:
		log(logfile,line)

	sources=[]
	sources.append("# the following were added during sassie installation")
	sources.append("# by "+getpass.getuser()+" on "+time.ctime())
	sources.append("deb http://archive.ubuntu.com/ubuntu hardy main")
	sources.append("deb-src http://archive.ubuntu.com/ubuntu hardy main")
	sources.append("deb http://archive.ubuntu.com/ubuntu hardy universe")
	sources.append("deb-src http://archive.ubuntu.com/ubuntu hardy universe")

	sourcefile = open('/etc/apt/sources.list','a')
	for source in sources:
		sourcefile.write("%s\n" % (source))
		log(logfile,'appending line to sources.list : '+source)
	sourcefile.close()

        if(version == '10.04'):
		result = os.popen('cp /etc/apt/apt.conf /etc/apt/apt.conf.backup').readlines()
		for line in result:
			log(logfile,line)
	
		apt_string = 'APT::Cache-Limit "99999999";'
		apt_file = open('/etc/apt/apt.conf','a')
		apt_file.write("%s\n" % (apt_string))
		apt_file.close()


	log(logfile,'>> updating system via apt-get update')

	result = os.popen('apt-get update').readlines()
	for line in result: log(logfile,line)

	if(version == '10.04'):
     	    result = os.popen('apt-get -y install g++').readlines()
	    for line in result: log(logfile,line)
        elif(version == '14.04'):
#            programs_needed.append('gnuplot-x11')
            programs_needed.append('libx32z1')
            programs_needed.append('libx32ncurses5')
            programs_needed.append('libbz2-dev')
            bit = platform.architecture()[0]
            if(bit == '64bit'):
     	       
                result = os.popen('dpkg --add-architecture i386').readlines()
	        for line in result: log(logfile,line)
     	        
                result = os.popen('apt-get -y install libc6:i386').readlines()
	        for line in result: log(logfile,line)
               
                debfile = open('/etc/apt/sources.list.d/ia32-libs-raring.list','a')
                debfile.write("deb http://archive.ubuntu.com/ubuntu/ raring main restricted universe multiverse\n")
                debfile.close()
     	        
                result = os.popen('apt-get -y update').readlines()
	        for line in result: log(logfile,line)

                result = os.popen('apt-get -y install ia32-libs').readlines()
	        for line in result: log(logfile,line)

                result = os.popen('rm /etc/apt/sources.list.d/ia32-libs-raring.list').readlines()
	        for line in result: log(logfile,line)
                
                result = os.popen('apt-get -y update').readlines()
	        for line in result: log(logfile,line)

                result = os.popen('apt-get -y install gcc-multilib --fix-missing').readlines()
	        for line in result: log(logfile,line)

	for program in programs_needed:
		error = apt_install(logfile,program)
		if(len(error)>0): print 'error = ',error

	return

def check_and_install_dependices(logfile,os_type,current_path,install_path):

	log(logfile,'\nCHECKING AND INSTALLING DEPENDENCIES\n')
	if(os_type[0] == 'Darwin'):
		darwin_dependencies(logfile,current_path,install_path)
	elif(os_type[1] == 'redhat' or os_type[1] == 'CentOS' or os_type[1] == 'CentOS Linux'):
		redhat_dependencies(logfile,current_path,install_path)
	elif(os_type[1] == 'Ubuntu'):
		ubuntu_dependencies(logfile,current_path,install_path,os_type[2])

	return

def install_toppar(logfile,current_path,install_path):

	cpst = 'cp -R '+current_path+'/sassie/distribution/general/toppar '+install_path

	result = os.popen(cpst).readlines()
	for line in result:
		log(logfile,line)

	if(not os.path.isdir(install_path+'/toppar')):
		print_error(logfile,"toppar did not install correctly!\n\nINSTALLATION STOPPING NOW\n\n")
		sys.exit()

	return

def compile_extensions(logfile,current_path,python):

	spath = current_path+'/sassie/sasmol/'
	path = current_path+'/sassie/sasmol/extensions/'
	dpath = path+'dcdio/'
	mpath = path+'mask/'
	vpath = path+'sasview/'
	mmpath = path+'matrix_math/'

	os.chdir(dpath)
	buildst1 = python +' setup_dcdio.py build'
	result = os.popen(buildst1).readlines() 
	for line in result: log(logfile,line)
	
	os.chdir(mpath)
	buildst2 = python +' setup_mask.py build'
	result = os.popen(buildst2).readlines() 
	for line in result: log(logfile,line)

        os.chdir(vpath)
        buildst3 = python +' setup_sasview_vmd.py build'
        result = os.popen(buildst3).readlines()
        for line in result: log(logfile,line)

        os.chdir(mmpath)
        buildst4 = python +' setup_matrix_multiply.py build'
        result = os.popen(buildst4).readlines()
        for line in result: log(logfile,line)

	os.chdir(current_path)	
	
	cpst1 = 'cp '+dpath+'build/lib*/_dcdio.so '+spath
	cpst2 = 'cp '+dpath+'dcdio.py '+spath
	cpst3 = 'cp '+mpath+'build/lib*/_mask.so '+spath
	cpst4 = 'cp '+mpath+'mask.py '+spath
	cpst5 = 'cp '+vpath+'build/lib*/_sasview_vmd.so '+spath
	cpst6 = 'cp '+vpath+'sasview_vmd.py '+spath
	cpst7 = 'cp '+mmpath+'build/lib*/matrix_math.so '+spath

	result = os.popen(cpst1).readlines() 
	for line in result: log(logfile,line)
	result = os.popen(cpst2).readlines() 
	for line in result: log(logfile,line)
	result = os.popen(cpst3).readlines() 
	for line in result: log(logfile,line)
	result = os.popen(cpst4).readlines() 
	for line in result: log(logfile,line)
	result = os.popen(cpst5).readlines() 
	for line in result: log(logfile,line)
	result = os.popen(cpst6).readlines() 
	for line in result: log(logfile,line)
	result = os.popen(cpst7).readlines() 
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
	if(not os.path.isfile(spath+'_sasview_vmd.so')):
       		print "_sasview_vmd.so did not install correctly!\n\nINSTALLATION STOPPING NOW\n\n"
    		sys.exit()
	if(not os.path.isfile(spath+'sasview_vmd.py')):
		print "sasview_vmd.py did not install correctly!\n\nINSTALLATION STOPPING NOW\n\n"
		sys.exit()
	if(not os.path.isfile(spath+'matrix_math.so')):
		print "matrix_math.so did not install correctly!\n\nINSTALLATION STOPPING NOW\n\n"
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

	install_toppar(logfile,current_path,install_path)
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
		'TAMD',
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

	log(logfile,'\nSASSIE INSTALLER: version 1.99 : 11/22/2014\n')

	user = getpass.getuser()
	current_path = os.getcwd()
	
	log(logfile,'executed by user : '+user+' on '+time.ctime())
	log(logfile,'current directory : '+current_path)
	log(logfile,'hostname : '+platform.node()+'\n')

	#s if(user != "root"):
	#s 	print_error(logfile,'\n\nYOU MUST RUN THIS SCRIPT AS ROOT\n\n')
	#s 	log(logfile,'\n\nINSTALLATION STOPPING NOW\n\n')
	#s 	sys.exit()

	log(logfile,'\n >> NOTE: This script may compile and install software')
	log(logfile,' >> including Python > 2.6.X, NumPy 1.6.1, SciPy 0.9.0')
	log(logfile,' >> Swig 2.0.4, gfortan 4.6.1, Gnuplot 4.4.3')
	log(logfile,'\n >> Installation on Linux machines may require an internet connection to download dependencies\n')
	log(logfile,' >> Apple installations will include Xcode 3.X or Xcode 4.X, EPD-free-7.X, & AquaTerm 1.1.0\n')
	log(logfile,' >> Apple installations may required you to install COMMAND LINE TOOLS\n')
	log(logfile,' >> The re-distributed software is open-source and existing license terms remain in effect\n\n')
	log(logfile,' >> The software will be installed in')
	log(logfile,'\n /usr/local/bin \n\n')
	log(logfile,' >> If you wish to install this software in your own Python distribution')
	log(logfile,' >> you probably can construct a method of your own from reading this installation script.\n')

	install_path = '/usr/local/bin'

	st1 = '\n >> press (0) to quit or (1) to continue or (2) to change installation path\n'
	decide=""
	while(decide != "0" and decide != "1" and decide != "2"):
                  decide = raw_input(st1)

	if(decide == "0"):
		log(logfile,'\n\nSTOPPING INSTALLATION NOW\n\n')
		sys.exit()
	elif(decide == "1"):
			error = check_path(install_path)
			if(len(error)>0):
				err = 'Default installation in /usr/local/bin is not valid : '
				err += error[0]
				print_error(logfile,err)
                                createpath=""
                                cst1 = '\nenter (0) to continue or (1) to try to create this path\n'
                                while(createpath != "0" and createpath != "1"):
                                        createpath = raw_input(cst1)

                                if(createpath == "1"):
                                        mkst = 'mkdir -p '+install_path
                                        try:
                                                result = os.popen(mkst).readlines()
                                                for line in result:
                                                        log(logfile,line)
                                        except:
                                                log(logfile,'\n > could not make that directory')

                                                log(logfile,'\n\n >> trying original option (2)\n\n')
                                                decide = "2"
                                else:
                                        log(logfile,'\n\n >> trying original option (2)\n\n')
                                        decide = "2"

	if(decide == "2"):
		pathinput=""
		st2 = '\n >> enter path for installation: '
		while(pathinput != "1"):
			newpath = raw_input(st2)
			error = check_path(newpath)
			if(len(error)==0):	
				pathinput = "1"
			else:
				print_error(logfile,error[0])
		log(logfile,'\n >> alternate installation path chosen : '+newpath)
		install_path = newpath
	
	return current_path,install_path

def epilogue(logfile,os_type):

	log(logfile,'\n\nINSTALLATION COMPLETE \n\n')
	log(logfile,'\n >> you can start sassie by running the script')
	log(logfile,' >> "sassie_gui.py" which is in a directory called')
	log(logfile,' >> "testing" in this directory')
	log(logfile,'\n >> you should copy this file to your own working directory')
	log(logfile,' >> to work on your own projects')

	if(os_type[1] ==  "redhat" or os_type[1] == "CentOS" or os_type[1] == "CentOS Linux"):
		log(logfile,'\n\n >> to run sassie on RedHat/CentOS machines')
		log(logfile,'>> you can use the new Python installation directly by typing\n')
		log(logfile,'\t > /usr/local/bin/bin/python2.7 sassie_gui.py\n\n')
		log(logfile,'>> or you can add a line to your .cshrc or .bashrc files\n')
		log(logfile,'  for csh/tcsh ')
		log(logfile,' set path = ( /usr/local/bin/bin $path )')
		log(logfile,' of for bash ')
		log(logfile,' export PATH=$PATH:/usr/local/bin/bin ')
		log(logfile,'>> make sure you either open a new terminal to access these changes')
		log(logfile,'>> or source the file from the command line')
		log(logfile,'\n\nthis will allow you to run sassie by typing\n')
		log(logfile,' python2.7 sassie_gui.py')
		log(logfile,'\n enjoy,\njec nist\n')


	return

def set_modules(os_type, **kwargs):

	full_module_list = ['Contrast Calculator','Align','Data Interpolation','Coordinate Tools','Merge Utilities','Monomer Monte Carlo','Complex Monte Carlo','TAMD','Energy Minimization','Docking','SasCalc','SLDMOL','EM to SAS','Chi-Square Filter','Density Plot','APBS'] 

	excluded_modules = []

	module_list = []

	for module in full_module_list:
		if module not in excluded_modules:
			module_list.append(module)
	
	return module_list

def install_me():

	logfile = open('log_sassie_install.txt','w')

	current_path,install_path = preamble(logfile)
	
	os_type,python = determine_os(logfile,install_path,current_path)

	module_list = set_modules(os_type)

	module_installation(logfile,current_path,install_path,module_list)

	#sys.exit()

	check_and_install_dependices(logfile,os_type,current_path,install_path)

	install_sassie(logfile,os_type,current_path,install_path,python)

	test_installation(logfile,os_type)

	epilogue(logfile,os_type)

	logfile.close()

	return


if __name__=='__main__':
      install_me()

