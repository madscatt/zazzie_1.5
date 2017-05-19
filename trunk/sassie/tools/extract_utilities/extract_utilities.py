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

#
#	EXTRACT UTILITIES
#
#	08/27/2012	--	initial coding			                        :	jc
#	09/06/2012	--	fixed for large DCD files	                    : 	hz
#	04/19/2015	--	added SAS options                               :   jc
#	10/09/2016	--	added sascalc and existing folder timestamps    :   jc
#
#	 1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#								       *      **
'''
	EXTRACT_UTILITIES is a module that allows one to extract coordinates and/or
    SAS profiles from PDB/DCD files and/or a directory containing SAS profiles.
	The multi-frame files trajectory is saved into new a PDB/DCD file.  SAS profiles
    are saved to a new directory.   Options are provided to extract
	single structures and/or SAS, structures and/or SAS over a range, by reading a list of frame 
	numbers from a simple text file, via a 'weights' file, or by a user-supplied
    sampling frequency.

	INPUT: 
	
		PDB or DCD file with molecular structure data (multi-frames)
		Extraction option
		Option value
		Name of output file		
        SAS profile type
        Folder containing input SAS profiles

	OUTPUT:
        	
		PDB or DCD file with requested coordinates
        and/or
        Folder containing SAS profiles

###
###	NOTE: should be straightforward to add basis filtering logic to pull out
###	      subsets.  Although this will take some work and may be better to
###           include with "builder" and/or "pdbrx" type projects
###
###	NOTE2: You could re-write this so that the single frame and range options
###	       work just like the text_file and weight_file versions (i.e. using masks)
###
###
'''

import os
import sys
import locale
import numpy
import string
import time
import glob
import sassie.interface.input_filter as input_filter
import sasmol.sasmol as sasmol
import sassie.util.folder_management as folder_management

import datetime


def print_failure(message, txtOutput):

    txtOutput.put("\n\n>>>> RUN FAILURE <<<<\n")
    txtOutput.put(">>>> RUN FAILURE <<<<\n")
    txtOutput.put(">>>> RUN FAILURE <<<<\n\n")
    txtOutput.put(message)

    return


def unpack_variables(variables):
    '''
    This method returns the variables to be used in the program.

    '''

    runname = variables['runname'][0]
    path = variables['path'][0]
    pdb_filename = variables['pdb_filename'][0]
    trajectory_filename = variables['trajectory_filename'][0]
    option = variables['option'][0]
    local_value = variables['local_value'][0]
    output_filename = variables['output_filename'][0]
    extract_trajectory = variables['extract_trajectory'][0]
    extract_sas = variables['extract_sas'][0]

    sas_type = variables['sas_type'][0]
    sas_paths = variables['sas_paths'][0]

    return runname, path, pdb_filename, trajectory_filename, option, local_value, output_filename, extract_trajectory, extract_sas, sas_type, sas_paths


def copy_spec_files(sas_files, runname, sas_output_path, suffix, mask):

    string_fill = 5

    number = 1
    number_of_saved_files = 0

    for name in sas_files:
        if str(number) in mask:
            number_of_saved_files += 1
            numst = str(number_of_saved_files).zfill(string_fill)
            runst = runname + '_' + numst + suffix[1:]
#            print 'runst: ', runst
            this_file = os.path.join(sas_output_path, runst)
            cpst = 'cp  ' + name + ' ' + this_file
#            print 'cpst: ', cpst
            os.system(cpst)

        number += 1

    return number_of_saved_files


def get_range_sas(local_value):

    local_values = string.split(local_value, '-')
    first = locale.atoi(local_values[0])
    last = locale.atoi(local_values[1])
    return range(first, last + 1)


def get_text_file_sas(local_value):

    infile = open(local_value, 'r').readlines()

    mask = []
    for i in xrange(len(infile)):
        lin = string.split(infile[i])
        if(len(lin) > 0):
            mask.append(lin[0])

    return mask


def get_weight_file_sas(local_value):
    weights = numpy.loadtxt(local_value)
    fweights = weights[:, 2]
    file_numbers = weights[:, 0]
    x2 = weights[:, 1]
    for i in xrange(len(x2)):
        if(fweights[i] == 1):
            #             print 'st = ',int(weights[i][0]),' : x2 = ',x2[i]
            pass
    frame_list = ((numpy.nonzero(fweights))[0]).tolist()
    #mask = [str(i) for i in frame_list]
    mask = [str(int(file_numbers[i])) for i in frame_list]
    return mask


def get_frequency(local_value, number_of_spec_files, **kwargs):

    mask = []

#    print 'nspecfil = ',number_of_spec_files

    try:
        coordinate_flag = kwargs['coordinate_flag']
#        print 'extracting coordinates: local_value = ',local_value

    except:
        coordinate_flag = False

    if not coordinate_flag:
        for number in xrange(1, number_of_spec_files + 1, int(local_value)):
            mask.append(str(number))
    elif coordinate_flag:
        for number in xrange(0, number_of_spec_files, int(local_value)):
            mask.append(number)

    return mask


def get_sas_mask(extract_option, value, **kwargs):

    mask = []

    if extract_option == 'single_frame':
        mask.append(value)
    elif extract_option == 'range':
        rangelist = get_range_sas(value)
        mask = [str(i) for i in rangelist]
    elif extract_option == 'text_file':
        mask = get_text_file_sas(value)
    elif extract_option == 'weight_file':
        mask = get_weight_file_sas(value)
    elif extract_option == 'sampling_frequency':
        number_of_spec_files = kwargs['number_of_spec_files']
        mask = get_frequency(value, number_of_spec_files)
    elif extract_option == 'all':
        number_of_spec_files = kwargs['number_of_spec_files']
        value = '1'
        mask = get_frequency(value, number_of_spec_files)

    # print mask

    return mask


def extract_sas_files(runname, sas_paths, sas_type, extract_option, output_path, local_value, txtOutput):

#added support for subdirectories under sascalc

    sas_paths = [x.strip() for x in sas_paths.split(',')]
    if len(sas_paths) > 1:
        base_paths = []
        print 'sas_paths = ',sas_paths
        for sas_path in sas_paths:
            base_paths.append(os.path.basename(os.path.normpath(sas_path)))

        print 'base_paths = ', base_paths

        print 'sas_paths = ',sas_paths
        sas_input_paths = sas_paths
    else:
        for sas_path in sas_paths:
            base_paths = []
            print 'sas_path: ',sas_path
            new_paths = []
            for root, dirs, files in os.walk(sas_path, topdown=False):
#               for name in files:
#                   print(os.path.join(root, name))

                for name in dirs:
                    print(os.path.join(root, name))
                    new_paths.append(os.path.join(root, name))
                    base_paths.append(os.path.basename(os.path.normpath(os.path.join(root, name))))

            print 'new_paths: ', new_paths
            print 'base_paths: ', base_paths
            if base_paths == []:
                base_paths.append(os.path.basename(os.path.normpath(sas_path)))
                new_paths.append(sas_path)
            print 'base_paths = ', base_paths
            print 'new_paths = ', new_paths
            sas_input_paths = new_paths

    sas_output_paths = []

    if(sas_type == 0):
        for base_path in base_paths:
            sas_output_paths.append(os.path.join(
                output_path, 'sascalc', base_path))
#        print 'sas_output_paths = ', sas_output_paths
        suffix = ['*.iq', '*.log']
        extra = []
        num_ex_files = [0, 0, 0, 0]
        # for root, dirs, files in os.walk(sas_path):
        #    all_sas_paths.append(dirs)
        # ['neutron_D2Op_100']
        #sas_path = os.path.join(sas_path, all_sas_paths[0][0])

    elif(sas_type == 1):
        sas_output_paths.append(os.path.join(output_path, 'xtal2sas'))
        suffix = ['*.iq', '*.log']
        extra = ['*.inf', '*.crd', '*.ans', '*.pr']
        num_ex_files = [1, 1, 1, 1]
    elif(sas_type == 2):
        sas_output_paths.append(os.path.join(output_path, 'cryson'))
        suffix = ['*.int', '*.log']
        extra = ['*.sav', '*.flm', '*.alm', '*.ans']
        num_ex_files = [1, 1, 1, 1]
    elif(sas_type == 3):
        sas_output_paths.append(os.path.join(output_path, 'crysol'))
        suffix = ['*.int', '*.log']
        extra = ['*.sav', '*.flm', '*.alm', '*.ans']
        num_ex_files = [1, 1, 1, 1]

#    utc_datetime = datetime.datetime.utcnow().strftime("%Y_%b_%d_%H:%M:%S")
#    new_directory_name = sas_output_path + '_' + utc_datetime + '_UTC'

    print 'sas_output_paths = ', sas_output_paths
    print len(sas_output_paths)

    print 'sas_input_paths = ', sas_input_paths 
    print len(sas_input_paths)


    for sas_output_path in sas_output_paths:
        direxist = os.path.exists(sas_output_path)
        if(direxist):
            try:
                modification_date = folder_management.modification_date(
                    sas_output_path)
                new_directory_name = sas_output_path + '_' + modification_date + '_UTC'
                result = os.system('mv ' + sas_output_path +
                                   ' ' + new_directory_name)
                result = os.system('mkdir -p ' + sas_output_path)
            except:
                message = 'can not create project directory: ' + sas_output_path
                message = '\n or can not create backup project directory: ' + \
                    new_directory_name + '_UTC'
                message += '\nstopping here\n'
                print_failure(message, txtOutput)
                sys.exit()
        else:
            result = os.system('mkdir -p ' + sas_output_path)

    copy_extras = False

    num_iq_files = 1
    num_log_files = 1

    count = 0

    for sas_input_path in sas_input_paths:
        sas_output_path = sas_output_paths[count]
        spec_files = []
        log_files = []
        extra_files = []

        for name in glob.glob(os.path.join(sas_input_path, suffix[0])):
            spec_files.append(name)
        for name in glob.glob(os.path.join(sas_input_path, suffix[1])):
            log_files.append(name)
        for ex in extra:
            this_extra = []
            for name in glob.glob(os.path.join(sas_input_path, ex)):
                this_extra.append(name)
            if(copy_extras == False and len(this_extra) > 0):
                copy_extras = True
                print 'copying extra sas files'
            this_extra.sort()
            extra_files.append(this_extra)
#            print 'extra files: ',extra_files

        spec_files.sort()
        log_files.sort()

        if(len(spec_files) == 0):
            print 'sas_input_path,suffix = ', os.path.join(sas_input_path, suffix[0])
            sys.exit()

        if(extract_option in ['sampling_frequency', 'all']):
            mask = get_sas_mask(extract_option, local_value,
                                number_of_spec_files=len(spec_files))

        else:
            mask = get_sas_mask(extract_option, local_value)

        num_iq_files = copy_spec_files(
            spec_files, runname, sas_output_path, suffix[0], mask)
        num_log_files = copy_spec_files(
            log_files, runname, sas_output_path, suffix[1], mask)


#        print 'num_iq_files = ', num_iq_files
#        print 'num_log_files = ', num_log_files

        if(copy_extras == True):
            for j in xrange(len(extra)):
                num_ex_files[j] = copy_spec_files(
                    extra_files[j], runname, sas_output_path, extra[j], mask)
#            print 'num_ex_files = ', num_ex_files

        txtOutput.put('wrote %i sas files to %s\n' %
                      (num_iq_files, sas_output_path))
#        print 'sas_input_path = ',sas_input_path
        print 'wrote %i sas files to %s\n' % (num_iq_files, sas_output_path)
        #output_log_file.write('wrote %i sas files to %s\n' % (num_iq_files,sas_output_path))
        count += 1

    return


def get_single_frame(local_value):
    frame = locale.atoi(local_value) - 1
    return [frame]


def get_range(local_value):

    local_values = string.split(local_value, '-')
    first = locale.atoi(local_values[0]) - 1
    last = locale.atoi(local_values[1]) - 1
    return range(first, last + 1)


def get_text_file(local_value):

    infile = open(local_value, 'r').readlines()

    frame_list = []
    for i in xrange(len(infile)):
        lin = string.split(infile[i])
        if(len(lin) > 0):
            this_value = locale.atoi(lin[0]) - 1
            frame_list.append(this_value)

    # print 'frame_list = ',frame_list
    return frame_list


def get_weight_file(local_value):
    weights = numpy.loadtxt(local_value)
    fweights = weights[:, 2]
    x2 = weights[:, 1]
    for i in xrange(len(x2)):
        if(fweights[i] == 1):
            #           print 'st = ',int(weights[i][0]),' : x2 = ',x2[i]
            pass
    frame_list = ((numpy.nonzero(fweights))[0]).tolist()
    return frame_list


def extract_coords(path, output_path, pdb_filename, trajectory_filename, infile_type, output_filename, txtOutput, rangelist):

    m1 = sasmol.SasMol(0)
    m1.read_pdb(pdb_filename)
    natoms = m1.natoms()

    output_filename = output_path + output_filename

    coor = numpy.zeros((1, natoms, 3), numpy.float32)

    txtOutput.put("writing frames to %s \n" % (output_filename))
    if(infile_type == 'dcd'):
        print '> input file is a dcd file'
        if(output_filename[-3:] == 'dcd'):
            print '> output file is a dcd file'
            j = 0
            m2 = sasmol.SasMol(1)
            m2.read_pdb(pdb_filename, fastread=True)
            coor = numpy.zeros((1, natoms, 3), numpy.float32)
            dcdfile = m1.open_dcd_read(trajectory_filename)
            number_of_frames = dcdfile[2]
            print 'number_of_frames = ', number_of_frames

#			print 'rangelist = ',rangelist
#			print 'max(rangelist) = ',max(rangelist)
            dcdoutfile = m2.open_dcd_write(output_filename)
            for i in range(number_of_frames):
                print '.',
                sys.stdout.flush()
                m1.read_dcd_step(dcdfile, i)
                if i in rangelist:
                    print '\nextracting coordinates from frame = ', i
                    coor[0, :, :] = m1.coor()[0]
                    m2.setCoor(coor)
                    m2.write_dcd_step(dcdoutfile, 0, j + 1)
                    j += 1
                if(i > max(rangelist) + 1):
                    break

            m2.close_dcd_write(dcdoutfile)
            m1.close_dcd_read(dcdfile[0])

        elif(output_filename[-3:] == 'pdb'):
            print '> output file is a pdb file'
            m2 = sasmol.SasMol(1)
            m2.read_pdb(pdb_filename)  # ,fastread = True)
            j = 0
            dcdfile = m1.open_dcd_read(trajectory_filename)
            number_of_frames = dcdfile[2]
            print 'number_of_frames = ', number_of_frames
            for i in range(number_of_frames):
                print '.',
                sys.stdout.flush()
                m1.read_dcd_step(dcdfile, i)
                coor[0, :, :] = m1.coor()[0]
                m2.setCoor(coor)
                if i in rangelist:
                    print '\nextracting coordinates from frame = ', i
                    if(j == 0):
                        m2.write_pdb(output_filename, 0, 'w', model=j + 1)
                    else:
                        m2.write_pdb(output_filename, 0, 'a', model=j + 1)
                    j += 1
                if(i > max(rangelist) + 1):
                    break

            with open(output_filename, "a") as myfile:
                myfile.write("END\n")

    elif(infile_type == 'pdb'):
        m1.read_pdb(trajectory_filename)
        natoms = m1.natoms()
    #	coor = numpy.zeros(((last-first+1),natoms,3),numpy.float32)
        coor = numpy.zeros((1, natoms, 3), numpy.float32)
        number_of_frames = m1.number_of_frames()

        print '> input file is a pdb file'
        if(output_filename[-3:] == 'dcd'):
            print '> output file is a dcd file'
            j = 0
            m2 = sasmol.SasMol(1)
            m2.read_pdb(trajectory_filename, fastread=True)
            dcdoutfile = m2.open_dcd_write(output_filename)
            for i in range(number_of_frames):
                print '.',
                sys.stdout.flush()
                if i in rangelist:
                    print '\nextracting coordinates from frame = ', i
                    coor[0, :, :] = m1.coor()[i]
                    m2.setCoor(coor)
                    m2.write_dcd_step(dcdoutfile, 0, j + 1)
                    j += 1
                if(i > max(rangelist) + 1):
                    break

        elif(output_filename[-3:] == 'pdb'):
            print '> output file is a pdb file'
            m2 = sasmol.SasMol(1)
            m2.read_pdb(pdb_filename, fastread=True)
            j = 0
            for i in range(number_of_frames):
                print '.',
                sys.stdout.flush()
                if i in rangelist:
                    print '\nextracting coordinates from frame = ', i

                    coor[0, :, :] = m1.coor()[i]
                    m2.setCoor(coor)
                    if(j == 0):
                        m2.write_pdb(output_filename, 0, 'w', model=j + 1)
                    else:
                        m2.write_pdb(output_filename, 0, 'a', model=j + 1)
                    j += 1
                if(i > max(rangelist) + 1):
                    break

            with open(output_filename, "a") as myfile:
                myfile.write("END\n")

    txtOutput.put("wrote %i frames to %s \n" %
                  (len(rangelist), output_filename))
    return


def get_number_of_frames(infile_type, trajectory_filename):

    m = sasmol.SasMol(0)

    if infile_type == 'pdb':
        m.read_pdb(trajectory_filename)
        number_of_frames = m.number_of_frames()
    elif infile_type == 'dcd':
        dcdfile = m.open_dcd_read(trajectory_filename)
        number_of_frames = dcdfile[2]

    return number_of_frames


def extract_data(variables, txtOutput):

    runname, path, pdb_filename, trajectory_filename, option, local_value, output_filename, extract_trajectory, extract_sas, sas_type, sas_path = unpack_variables(
        variables)

    pdb_filename = path + pdb_filename
    trajectory_filename = path + trajectory_filename

    output_path = runname + '/extract_utilities/'
#    print 'output_path = ',output_path

    vers = 'version 0.2 : 04/23/15 : jc'
    direxist = os.path.exists(output_path)
    if(direxist == 0):
        try:
            result = os.system('mkdir -p ' + output_path)
        except:
            message = 'can not create project directory: ' + output_path
            message += '\nstopping here\n'
            print_failure(message, txtOutput)
            sys.exit()

    # ttxt=time.ctime()
    ttxt = time.asctime(time.gmtime(time.time()))
    st = ''.join(['=' for x in xrange(60)])

    txtOutput.put("\n%s \n" % (st))
    txtOutput.put("DATA FROM RUN: %s \n\n" % (ttxt))

    if extract_trajectory:

        infile_type = trajectory_filename[-3:]

        txtOutput.put("reading frames from %s \n" % (trajectory_filename))

        if extract_sas:
            fraction_done = 0.25
        else:
            fraction_done = 0.50

        report_string = 'STATUS\t' + str(fraction_done)
        txtOutput.put(report_string)

        if(option == 'single_frame'):

            rangelist = get_single_frame(local_value)

        elif(option == 'range'):

            rangelist = get_range(local_value)

        elif(option == 'text_file'):

            rangelist = get_text_file(local_value)

        elif(option == 'weight_file'):

            rangelist = get_weight_file(local_value)

        elif(option == 'sampling_frequency'):
            number_of_frames = get_number_of_frames(
                infile_type, trajectory_filename)
#            print 'local_value = ',local_value
            rangelist = get_frequency(
                local_value, number_of_frames, coordinate_flag=True)

        elif(option == 'all'):
            local_value = 1
            number_of_frames = get_number_of_frames(
                infile_type, trajectory_filename)
            rangelist = get_frequency(
                local_value, number_of_frames, coordinate_flag=True)

        extract_coords(path, output_path, pdb_filename, trajectory_filename,
                       infile_type, output_filename, txtOutput, rangelist)

    if extract_sas:

        txtOutput.put("\nextracting SAS profiles\n")
        #txtOutput.put("sas_type = %s\n" % (str(sas_type)))
        txtOutput.put("sas_path = %s\n" % (sas_path))

        extract_sas_files(runname, sas_path, sas_type, option,
                          output_path, local_value, txtOutput)

        fraction_done = 0.5
        report_string = 'STATUS\t' + str(fraction_done)
        txtOutput.put(report_string)

    fraction_done = 1.0
    report_string = 'STATUS\t' + str(fraction_done)
    txtOutput.put(report_string)

    txtOutput.put("\n%s \n" % (st))
    print '\nEXTRACT UTILITIES IS DONE'
    time.sleep(2.5)

    return

if __name__ == '__main__':

    # BEGIN USER EDIT
    # BEGIN USER EDIT
    # BEGIN USER EDIT

    runname = 'run_0'
    #pdb_filename ='min3.pdb'
    pdb_filename = 'hiv1_gag.pdb'
    #trajectory_filename = 'c7.dcd'
    #trajectory_filename = 'min3.pdb'
    trajectory_filename = 'hiv1_gag_200_frames.dcd'
    trajectory_filename = 'demo/run_0/monomer_monte_carlo/run_0.dcd'
    #option ='single_frame'
    option = 'range'
    #option ='weight_file'
    #local_value = '1'
    local_value = '3-22'
    #local_value = 'weights_file.txt'
    output_filename = 'output.dcd'
    #output_filename = 'output.pdb'
    extract_trajectory = 'True'
    extract_trajectory = 'False'
    extract_sas = 'True'
    path = './'
    sas_type = '0'  # sascalc
#    sas_type = '1' #xtal2sas
#    sas_type = '2' #cryson
#    sas_type = '3' #crysol
    sas_paths = 'run_0/sascalc/neutron_D2Op_100'
    sas_paths = 'demo/run_0/sascalc/neutron_D2Op_100,demo/run_0/sascalc/neutron_D2Op_80, demo/run_0/sascalc/xray'
#    sas_paths = 'demo/run_0/xtal2sas'
#    sas_paths = 'demo/run_0/cryson'
#    sas_paths = 'demo/run_0/crysol'

    # END USER EDIT
    # END USER EDIT
    # END USER EDIT

    svariables = {}

    svariables['runname'] = (str(runname), 'string')
    svariables['pdb_filename'] = (str(pdb_filename), 'string')
    svariables['trajectory_filename'] = (str(trajectory_filename), 'string')
    svariables['option'] = (str(option), 'string')
    svariables['local_value'] = (str(local_value), 'string')
    svariables['output_filename'] = (str(output_filename), 'string')
    svariables['extract_trajectory'] = (str(extract_trajectory), 'boolean')
    svariables['extract_sas'] = (str(extract_sas), 'boolean')
    svariables['path'] = (str(path), 'string')
    svariables['sas_type'] = (str(sas_type), 'int')
    svariables['sas_paths'] = (str(sas_paths), 'string')

    error, variables = input_filter.type_check_and_convert(svariables)
    if len(error) > 0:
        print 'error = ', error
        sys.exit()

    import sassie.interface.extract_utilities_filter as extract_utilities_filter

    error = extract_utilities_filter.check_extract_utilities(variables)
    if len(error) > 0:
        print 'error = ', error
        sys.exit()

    import multiprocessing
    txtQueue = multiprocessing.JoinableQueue()

    # process=multiprocessing.Process(target=extract_utilities.extract_data,args=(variables,txtQueue))
    # process.start()
    extract_data(variables, txtQueue)
