'''
    SASSIE  Copyright (C) 2011-2017 Joseph E. Curtis
    This program comes with ABSOLUTELY NO WARRANTY;
    This is free software, and you are welcome to redistribute it under certain
    conditions; see http://www.gnu.org/licenses/gpl-3.0.html for details.
'''
import os
import sys
import string
import platform
import time
import getpass
import readline

sys.path.append('./sassie/util/')
import sasconfig

#       INSTALLER
#
#       09/05/2011      --      initial coding              :   jc
#       06/27/2012      --      courtesy installer          :   jc
#       12/12/2013      --      added sasview/tamd          :   jc
#       02/19/2014      --      osx mavericks 10.9          :   jc
#       11/22/2014      --      updated for 2.0 branch      :   jc
#       05/24/2017      --      updated for 1.5 HPC branch  :   jc
#
# LC      1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''
	installer.py is the script to install sassie and required dependencies
	
	Ubuntu (14.04 LTS)

	See the README file in this directory for information on how
	to use this installer script.
'''


def log(logfile, line):

    print line

    logfile.write("%s\n" % (line))

    return


def print_error(logfile, error):

    log(logfile, '\n>>> ERROR detected : ' + error)

    return


def determine_os(logfile, install_path, current_path):

    supported_flag = 'NO'

    log(logfile, 'platform.system() = ' + platform.system())
    log(logfile, 'platform.platform() = ' + platform.platform())
    log(logfile, 'platform.processor() = ' + platform.processor())
    log(logfile, 'platform.architecture()[0] = ' + platform.architecture()[0])
    log(logfile, 'platform.architecture()[1] = ' + platform.architecture()[1])

    bit = platform.architecture()[0]

    local_os_type = platform.system()

    if(local_os_type == "Darwin"):
        log(logfile, 'Apple OS Detected')
        (release, versioninfo, machine) = platform.mac_ver()
        log(logfile, 'release = ' + release)
        log(logfile, 'machine = ' + machine)

        os_type = [local_os_type, release, machine, platform.architecture()[0]]

        error = 'this installer does not handle your OS'
        print_error(logfile, error)

    elif(local_os_type == "Linux"):

        log(logfile, 'Linux OS Detected')

        if(bit == '32bit'):
            error = '\n\nthis installer has not been tested for Linux 64-bit machines'
            print_error(logfile, error)

        try:
            (distname, version, id) = platform.linux_distribution()
        except:
            (distname, version, id) = platform.dist()

        log(logfile, 'distname = ' + distname)
        log(logfile, 'version = ' + version)
        log(logfile, 'id = ' + id)

        try:
            (lib, gversion) = platform.libc_ver()
            log(logfile, 'glib  = ' + lib + ' version = ' + gversion)
        except:
            log(logfile, 'could not determine glib used to compile current python')

        if(distname == "debian" and version == "jessie/sid"):
            supported_flag = 'YES'
            python = sys.executable
            log(logfile, 'installing using PYTHON = '+sys.executable)
        else:
            error = 'this installer does not handle your OS'
            log(logfile, error)
            log(logfile, 'distname = ' + distname)
            log(logfile, 'version = ' + version)
            print_error(logfile, error)

        os_type = [local_os_type, distname, version, platform.architecture()[
            0]]

    if(supported_flag == 'NO'):
        log(logfile, '\n >> your OS is not currently supported: but the installer may work\n')
        log(logfile, ' >> if you have a "flavor" of Debian\n')

        try_anyway = ""

        st0 = '\n \t\t (0) QUIT NOW \n'
        st2 = '\n \t\t (1) DEBIAN TYPE (Debian, Ubuntu, Linux Mint, etc.)\n'
        st3 = '\n >> quit (0) or enter OS type do you want to mimic ? : (1) '
        log(logfile, st0)
        log(logfile, st2)

        while(try_anyway != "0" and try_anyway != "1"):
            try_anyway = raw_input(st3)

        st4 = '\n >> you entered ' + try_anyway

        log(logfile, st4)

        if(try_anyway == "0"):
            log(logfile, '\n\n QUITING NOW \n\n')
            sys.exit()
        elif(try_anyway == "1"):
            local_os_type = "Linux"
            distname = "Ubuntu"
            (distname, version, id) = platform.dist()
            os_type = [local_os_type, distname, version, platform.architecture()[
                0]]
            supported_flag = 'MAYBE'
            python = sys.executable

    log(logfile, '>OS_TYPE = ' + os_type[0])
    log(logfile, '>IS OS SUPPORTED? = ' + supported_flag)

    return os_type, python


def find_files(logfile, files_to_check):

    programs_needed = []

    for i in xrange(len(files_to_check)):
        locatefile = os.popen('which ' + files_to_check[i]).readlines()
        if(locatefile):
            log(logfile, 'found ' + locatefile[0])
        else:
            log(logfile, 'file ' + files_to_check[i] + ' not found!')
            programs_needed.append(files_to_check[i])

    if(len(programs_needed) > 0):
        log(logfile, '\n >> dependencies summary: programs needed = ' +
            ', '.join(programs_needed))
    else:
        log(logfile, '\n >> dependencies summary:\n >> it seems that you have the programs needed to continue the installation\n')

    return programs_needed


def apt_install(logfile, program):

    log(logfile, '>>installing ' + program + ' using apt')
    error = []
    st1 = 'apt-get -y install ' + program
    try:
        result = os.popen(st1).readlines()
        for line in result:
            log(logfile, line)
    except:
        error.append = '>>> failed to install ' + \
            program + '\nINSTALLATION STOPPING\n\n'
        print_error(error)
        sys.exit()

    return error

def ubuntu_dependencies(logfile, current_path, install_path, version):

    os_type = "Ubuntu"
    os_type = "Linux"

    log(logfile, ' >> checking Ubuntu installation\n')
    if(version == 'jessie/sid'):
        files_to_check = ['gcc', 'g++', 'gfortran',
                          'swig']

    else:
        files_to_check = ['gcc', 'g++', 'gfortran', 'swig']

    programs_needed = find_files(logfile, files_to_check)

    files_to_append = ['libfftw3-3', 'libfftw3-dev']
                       
    for file in files_to_append:
        programs_needed.append(file)

    #result = os.popen(
    #    'cp /etc/apt/sources.list /etc/apt/sources.list.backup').readlines()
    #for line in result:
    #    log(logfile, line)

    sources = []
    #sources.append("# the following were added during sassie installation")
    #sources.append("# by " + getpass.getuser() + " on " + time.ctime())
    #sources.append("deb http://archive.ubuntu.com/ubuntu hardy main")
    #sources.append("deb-src http://archive.ubuntu.com/ubuntu hardy main")
    #sources.append("deb http://archive.ubuntu.com/ubuntu hardy universe")
    #sources.append("deb-src http://archive.ubuntu.com/ubuntu hardy universe")

    #sourcefile = open('/etc/apt/sources.list', 'a')
    #for source in sources:
    #    sourcefile.write("%s\n" % (source))
    #    log(logfile, 'appending line to sources.list : ' + source)
    #sourcefile.close()

    if(version == 'jessie/sid'):
        ##            programs_needed.append('gnuplot-x11')
        #programs_needed.append('libx32z1')
        #programs_needed.append('libx32ncurses5')
        #programs_needed.append('libbz2-dev')
        bit = platform.architecture()[0]

    print 'APT NOT ENABLED'
    print 'programs needed = ', programs_needed

    #for program in programs_needed:
    #    error = apt_install(logfile, program)
    #    if(len(error) > 0):
    #        print 'error = ', error

    return


def check_and_install_dependices(logfile, os_type, current_path, install_path):

    log(logfile, '\nCHECKING AND INSTALLING DEPENDENCIES\n')
    if(os_type[1] == 'debian'):
        ubuntu_dependencies(logfile, current_path, install_path, os_type[2])
    else:
        log(logfile, 'unable to install dependencies for os_type[1] = ' + os_type[1])

    return


def install_toppar(logfile, current_path, install_path):

    cpst = 'cp -R ' + current_path + \
        '/sassie/distribution/general/toppar ' + install_path

    result = os.popen(cpst).readlines()
    for line in result:
        log(logfile, line)

    if(not os.path.isdir(install_path + '/toppar')):
        print_error(
            logfile, "toppar did not install correctly!\n\nINSTALLATION STOPPING NOW\n\n")
        sys.exit()

    return


def modify_sascalc_library_makefile(logfile):

    log(logfile, 'writing new Makefile with data in sasconfig.py')

    original_makefile = 'sassie/calculate/sascalc/sascalc_library/cpp_and_cuda_source/Makefile'
    new_makefile = 'sassie/calculate/sascalc/sascalc_library/cpp_and_cuda_source/Makefile_new'

    os.system('cp '+ original_makefile + ' ' + original_makefile + '_backup')

    lines = [] 

    lines.append('INCDIR=-I' + sasconfig.__cuda_include_path__ + '\n')
    lines.append('CC=' + sasconfig.__cuda_gpp__ + '\n')
    lines.append('NVCC=' + sasconfig.__cuda_nvcc__ + '\n')

    with open(original_makefile, 'r') as original_file:
        with open(new_makefile, 'w') as new_file:
            for line in lines:
                new_file.write(line)
            new_file.write(original_file.read()) 

    os.rename(new_makefile, original_makefile)

    return

def compile_sascalc_library(logfile, current_path):

    sascalc_library_path = current_path + '/sassie/calculate/sascalc/sascalc_library'
    sascalc_cpp_and_cuda_source = sascalc_library_path + '/cpp_and_cuda_source'

    os.chdir(sascalc_cpp_and_cuda_source)

    result = os.popen('mkdir lib').readlines()
    for line in result:
        log(logfile, line)

    result = os.popen('make').readlines()
    for line in result:
        log(logfile, line)

    return

def compile_extensions(logfile, current_path, python):

    
    modify_sascalc_library_makefile(logfile)
    compile_sascalc_library(logfile, current_path)

    return


def compile_sassie(logfile, current_path, python):

    buildst = python + ' setup.py build'

    installst = python + ' setup.py install'

    os.chdir(current_path)
    result = os.popen(buildst).readlines()
    for line in result:
        log(logfile, line)
    result = os.popen(installst).readlines()
    for line in result:
        log(logfile, line)

    return

def install_sascalc_library(logfile, current_path, python):

    os.chdir(current_path)
    os.chdir('sassie/calculate/sascalc/sascalc_library/cpp_extension')
    installst = python + ' setup_sascalc_api.py install'
    result = os.popen(installst).readlines()
    for line in result:
        log(logfile, line)

    os.chdir(current_path)

    return

def install_sassie(logfile, os_type, current_path, install_path, python):

    install_toppar(logfile, current_path, install_path)
    compile_extensions(logfile, current_path, python)
    compile_sassie(logfile, current_path, python)
    install_sascalc_library(logfile, current_path, python)

    return

def clean_up(os_type):

    return


def check_permissions(path):

    try:
        existvalue = os.access(path, os.F_OK)
    except:
        existvalue = False
    try:
        readvalue = os.access(path, os.R_OK)
    except:
        readvalue = False
    try:
        writevalue = os.access(path, os.W_OK)
    except:
        writevalue = False

    return existvalue, readvalue, writevalue


def check_path(path):

    error = []
    ev, rv, wv = check_permissions(path)
    if(not ev or not rv or not wv):
        error.append('permission error in input file path ' + path + '  [PATH = ' + str(
            ev) + ': READ PERMISSION = ' + str(rv) + ': WRITE PERMISSION = ' + str(wv) + ']')
        if(ev == False):
            error.append('path does not exist')
        elif(rv == False):
            error.append('read permission not allowed')
        elif(wv == False):
            error.append('write permission not allowed')
        return error

    return error


def test_installation(logfile, os_type):

    log(logfile, '\nTESTING INSTALLATION\n')

    return


def preamble(logfile):

    log(logfile, '\nSASSIE INSTALLER: 1.5 HPC : 05/24/2017\n')

    user = getpass.getuser()
    current_path = os.getcwd()

    log(logfile, 'executed by user : ' + user + ' on ' + time.ctime())
    log(logfile, 'current directory : ' + current_path)
    log(logfile, 'hostname : ' + platform.node() + '\n')

    if(user != "root"):
        print_error(logfile,'\n\nYOU MUST RUN THIS SCRIPT AS ROOT\n\n')
        log(logfile,'\n\nINSTALLATION STOPPING NOW\n\n')
        sys.exit()

    install_path = sasconfig.__installation_bin_path__

    log(logfile, '\n >> NOTE: This script will compile and install software')
    log(logfile, '\n >> Several dependencies are required to be installed\n \
      before using this script: see https://github.com/madscatt/zazzie_1.5\n \
      for instructions\n')
    log(logfile, ' >> Redistributed software is open-source and existing license terms remain in effect\n\n')
    log(logfile, ' >> The software will be installed in')
    log(logfile, '\n ' + install_path + '\n\n')

    st1 = '\n >> press (0) to quit or (1) to continue \n'
    decide = ""
    while(decide != "0" and decide != "1"):
        decide = raw_input(st1)

    if(decide == "0"):
        log(logfile, '\n\nSTOPPING INSTALLATION NOW\n\n')
        sys.exit()
    elif(decide == "1"):
        error = check_path(install_path)
        if(len(error) > 0):
            err = 'Default installation in ' + install_path + ' is not valid : '
            err += error[0]
            print_error(logfile, err)
            createpath = ""
            cst1 = '\nenter (0) to continue or (1) to try to create this path\n'
            while(createpath != "0" and createpath != "1"):
                createpath = raw_input(cst1)

            if(createpath == "1"):
                mkst = 'mkdir -p ' + install_path
                try:
                    result = os.popen(mkst).readlines()
                    for line in result:
                        log(logfile, line)
                except:
                    log(logfile, '\n > could not make that directory')
                    log(logfile, '\n\nSTOPPING INSTALLATION NOW\n\n')
                    sys.exit()
                    
    return current_path, install_path

def epilogue(logfile, os_type):

    log(logfile, '\n\nINSTALLATION COMPLETE \n\n')

    return

def install_me():

    logfile = open('log_sassie_install.txt', 'w')

    current_path, install_path = preamble(logfile)

    os_type, python = determine_os(logfile, install_path, current_path)

    check_and_install_dependices(logfile, os_type, current_path, install_path)

    install_sassie(logfile, os_type, current_path, install_path, python)

    sys.exit()

    test_installation(logfile, os_type)

    epilogue(logfile, os_type)

    logfile.close()

    return


if __name__ == '__main__':
    install_me()
