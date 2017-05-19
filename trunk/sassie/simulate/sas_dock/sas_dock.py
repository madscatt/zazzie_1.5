import os
import sys
import string
import time
import numpy
import subprocess
import multiprocessing
import sasmol.sasmol as sasmol
import sassie.util.sasconfig as sasconfig
import sassie.util.module_utilities as module_utilities
import sassie.interface.input_filter as input_filter
from sassie.simulate.sas_dock.convergence_test import calc_spatial_convergence_all
from sassie.simulate.sas_dock.write_flags import *

app = 'sas_dock'

'''

    docking is a module that performs RosettaDock full protocol docking , and generates
    a dcd file with docking decoy sructures, in which protein partner 2 is globally moved around
    protein partner 1.

'''

class module_variables():

    def __init__(self, parent=None):
            self.app = app

bin_path = sasconfig.__bin_path__

rosetta_path = os.path.join(bin_path,'ROSETTA_3.4','rosetta_bin_linux_2015.19.57819_bundle','main','source','bin')
prepack_exe = os.path.join(rosetta_path, 'docking_prepack_protocol.linuxgccrelease')
docking_exe = os.path.join(rosetta_path, 'docking_protocol.default.linuxgccrelease')
fixbb_exe = os.path.join(rosetta_path, 'fixbb.default.linuxgccrelease')
database_path = ' -database ' + os.path.join(bin_path , 'ROSETTA_3.4', 'rosetta_bin_linux_2015.19.57819_bundle', 'main', 'database') + ' -nodelay  2>&1 ' 

class simulation():

    def __init__(self, parent=None):
        pass

    def main(self, input_variables, txtOutput):

        self.mvars = module_variables()

        self.run_utils = module_utilities.run_utils(app, txtOutput)

        self.run_utils.setup_logging(self)

        self.log.debug('in main')

        self.unpack_variables(input_variables)

        self.run_utils.general_setup(self)

        self.docking_lr(input_variables, txtOutput)

        self.sascalc_x2(input_variables, txtOutput)

        self.epilogue()

        return


    def unpack_variables(self, variables):
    
            mvars = self.mvars
            mvars.runname = variables['runname'][0]
            mvars.pdbfile = variables['pdbfile'][0]
            mvars.number_processors = variables['number_processors'][0]
            mvars.maximum_dock_cycles = variables['maximum_dock_cycles'][0]
            mvars.partner1_partner2 = variables['partner1_partner2'][0]
            mvars.sas_option = variables['sas_option'][0]
            mvars.interpolated_data_file = variables['interpolated_data_file'][0]
            mvars.spatial_voxel_convergence_tolerance = variables['spatial_voxel_convergence_tolerance'][0]
            mvars.number_of_q_values = variables['number_of_q_values'][0]
            mvars.maximum_q_value = variables['maximum_q_value'][0]
            mvars.i0 = variables['i0'][0]
    
            return 
    
    
    def print_failure(message, txtOutput):
    
            txtOutput.put("\n\n>>>> RUN FAILURE <<<<\n")
            txtOutput.put(">>>> RUN FAILURE <<<<\n")
            txtOutput.put(">>>> RUN FAILURE <<<<\n\n")
            txtOutput.put(message)
    
            return
    
    
    
    def docking_lr(self, variables, txtOutput):
        '''
        DOCKING DRIVER is the function to read in variables from GUI input and
        used to run docking with Rosetta program for two protein partners
    
        INPUT:  variable descriptions:
    
        runname:                                        directory path to write all output
        pdbfile:                                        input pdb file (reference)
        number_processors:                              number of core processors requested
        partner1_partner2:                              protein chain (partner2) to move around psuedo-fixed
                                                        protein chain(s) (partner1) (partner2)
    
        OUTPUT:
                                            runname.dcd 
                                            docking_scores
    
        txtOutput:        TK handler for output to GUI textbox
    
        files stored in ~/runname/sas_dock directory:
        outfile:          runname.dcd
                          fullatom_reference.pdb
                          docking_scores
    
        '''
    
        log = self.log
        pgui = self.run_utils.print_gui
    
        # start gui output
        pgui("\n%s \n" % ('=' * 60))
        pgui("DATA FROM RUN: %s \n\n" % time.asctime( time.gmtime( time.time() ) ))

        mvars = self.mvars
    
        path = os.path.join(mvars.runname, 'work')
        direxist = os.path.exists(path)
        if(direxist == 0):
            try:
                result = os.system('mkdir -p ' + path)
            except:
                message = 'can not create project directory: ' + path
                message += '\nstopping here\n'
                print_failure(message, txtOutput)
            if(result != 0):
                message = 'can not create project directory: ' + path
                message += '\nstopping here\n'
                print_failure(message, txtOutput)
    
    
        write_flags(mvars.pdbfile, mvars.partner1_partner2, path)
    
        def p(run_cmd, log_write):
            
            dock=subprocess.Popen(run_cmd, shell=True, executable='/bin/bash')
            dock.wait()
            if log_write:
                log.info('Finished Rosetta Docking %s' % run_cmd)
    
            return
    
        # prepack sidechains of initial structure--generates pdbfile_prepack.pdb
        # display progress
        fraction_done = (0 + 1) * 1.0 / 10.0
        report_string = 'STATUS\t%f' % fraction_done
        pgui(report_string)
    
        run_cmd=prepack_exe + ' -no_nstruct_label -s ' + mvars.pdbfile + ' @' + path + os.sep + 'flags_pre ' \
        '-no_nstruct_label -no_scores_in_pdb -out:path:all ' + path + \
            database_path + ' > ' + path + os.sep + 'prepack.log'
        log.info('Starting Rosetta Prepack %s' % run_cmd)
        prepack=subprocess.Popen(run_cmd, shell=True, executable='/bin/bash')
        prepack.wait()
    
        log.info('Finished Rosetta Prepack %s' % run_cmd)
        cmd='mv ' + path + os.sep + 'prepack__' + \
            mvars.pdbfile + ' ' + path + os.sep + 'dock.pdb'
        os.system(cmd)
    
        #set reference pdb for merging full-atom pdb files
        mvars.refpdb = path + os.sep + 'dock.pdb'

        cmd='rm -f ' + path + os.sep + 'away*.pdb'
        os.system(cmd)
        cmd='rm -f ' + path + os.sep + 'initial*.pdb'
        os.system(cmd)

        # display progress
        fraction_done = (1 + 1) * 1.0 / 10.0
        report_string = 'STATUS\t%f' % fraction_done
        pgui(report_string)
        #number_docking_structures_cycle = 10000
        number_docking_structures_percycle = 8
        number_docking_structures_processor = int(number_docking_structures_percycle/mvars.number_processors)

        #Start cycles of Approx. 10000 docking decoy structures
        #mvars.maximum_dock_cycles = 0
        for dock_cycle in xrange(mvars.maximum_dock_cycles):
            #Docking--generate number of decoy structures = number_processors x number_docking_structures_processor
            log.info('Begin Docking Cycle: %s' % int(dock_cycle + 1))
            
            processes = []
            for i in xrange(mvars.number_processors):
                istr='%04d' % (i + 1)
                seed=numpy.random.randint(1000000, high=10000000, size=None)
                run_cmd=docking_exe + ' -s ' + path + os.sep + 'dock.pdb -suffix _' + str(i + 1) +' -run:jran ' + str(seed) + ' @' + path + os.sep + \
                'flags -nstruct ' + str(number_docking_structures_processor) + ' -out:file:scorefile score_' + istr + ' -out:path:all ' \
                + path + database_path + '| egrep -vf ' + path + os.sep + 'ignore_list >' + path + os.sep + 'dock_' + str(i + 1) + '.log'
        
                log.info('Starting Rosetta Docking %s' % run_cmd)
        
                log_write=True
                dock = multiprocessing.Process(target=p, args=(run_cmd, log_write,))
                processes.append(dock)
        
            for d in processes:
                d.start()
        
            for d in processes:
                d.join(timeout=54000)
        
        
            # display progress
            fraction_done = (5 + 1) * 1.0 / 10.0
            report_string = 'STATUS\t%f' % fraction_done
            pgui(report_string)
        
        
            #collect scores files
            cmd = 'cat ' + path + os.sep + 'score_* >> ' + path + os.sep + 'docking_scores.txt'
            os.system(cmd)
            cmd = 'grep -v SEQUENCE ' + path + os.sep + 'docking_scores.txt | grep -v interchain > '  + path + os.sep + 'scores'
            os.system(cmd)
            cmd = 'mv ' + path + os.sep + 'scores ' + path + os.sep + 'docking_scores.txt'
            os.system(cmd)
            mvars.scores_file = path + os.sep + 'docking_scores.txt'


            trajectory_names = ''
            for i in xrange(mvars.number_processors):
                for j in xrange(number_docking_structures_processor):
                    if(number_docking_structures_processor < 100000):
                      jstr='%05d' % (j + 1)
                    if(number_docking_structures_processor < 10000):
                      jstr='%04d' % (j + 1)
                    if(i == 0 and j == 0):
                        trajectory_names=trajectory_names + \
                        path + os.sep + 'dock_' + str(i + 1) + '_' + jstr + '.pdb'
                    else:
                        trajectory_names=trajectory_names + ',' + \
                        path + os.sep + 'dock_' + str(i + 1) + '_' + jstr + '.pdb'
        
            #merge  full-atom pdb files from current docking cycle
            output_log_file=open(os.path.join(path, 'merge.log'), 'w')
            trajectory_names=string.split(trajectory_names, ',')
            m1 = sasmol.SasMol(0)
            m2 = sasmol.SasMol(0)
            m2.read_pdb(mvars.refpdb,fastread=True)
            m1.read_pdb(mvars.refpdb)
            mvars.dcdfile = mvars.runname + '_cycle.dcd'
            log.info('opening new dcd file to store trajectory: %s' %
                        os.path.join(self.runpath, mvars.dcdfile))
        
            outfile_name = str(os.path.join(path, mvars.dcdfile))
            dcdoutfile = m2.open_dcd_write(outfile_name)
            count = 0
            coor = numpy.zeros((1,m2.natoms(),3),numpy.float32)
            for this_trajectory_name in trajectory_names:
        
                    pdbfile = m1.read_pdb(this_trajectory_name)
                    number_of_frames = m1.number_of_frames()
        
                    for j in xrange(number_of_frames):
                        coor[0,:,:] = m1.coor()[j]
                        m2.setCoor(coor)
                        m2.write_dcd_step(dcdoutfile,0, count + 1)
                        count += 1
        
            m2.close_dcd_write(dcdoutfile)
        
            if (dock_cycle == 0):
                cmd = 'mv ' + outfile_name + ' ' + path + os.sep + mvars.runname + '_dock.dcd'
                os.system(cmd)
            else:
                #merge cummulative dcd file with current cycle dcd file
                trajectory_cycle_names = path + os.sep + mvars.runname + '_dock.dcd,' + outfile_name

                output_log_file=open(os.path.join(path, 'merge_cycles.log'), 'w')
                trajectory_cycle_names=string.split(trajectory_cycle_names, ',')
                m1 = sasmol.SasMol(0)
                m2 = sasmol.SasMol(0)
                m2.read_pdb(mvars.refpdb,fastread=True)
                m1.read_pdb(mvars.refpdb)
        
                mvars.dcdfile = mvars.runname + '_merge.dcd'
                log.info('opening new dcd file to store trajectory: %s' %
                         os.path.join(self.runpath, mvars.dcdfile))
        
                outfile_name = str(os.path.join(path, mvars.dcdfile))
                dcdoutfile = m2.open_dcd_write(outfile_name)
                count = 0
                coor = numpy.zeros((1,m2.natoms(),3),numpy.float32)
                for this_trajectory_name in trajectory_cycle_names:
        
                    dcdfile = m1.open_dcd_read(this_trajectory_name)
                    number_of_frames = dcdfile[2]
        
                    for j in xrange(number_of_frames):
                        m1.read_dcd_step(dcdfile,j)
                        coor[0,:,:] = m1.coor()[0]
                        m2.setCoor(coor)
                        m2.write_dcd_step(dcdoutfile,0, count + 1)
                        count += 1
        
                m2.close_dcd_write(dcdoutfile)
                
                cmd = 'mv ' + outfile_name + ' ' + path + os.sep + mvars.runname + '_dock.dcd'
                os.system(cmd)
                
                
                cmd='rm -f ' + path + os.sep + '*voxels*'
                os.system(cmd)
                cmd='rm -f ' + path + os.sep + '*.eps'
                os.system(cmd)
                cmd='rm -f ' + path + os.sep + '*.png'
                os.system(cmd)

            # display progress
            fraction_done = (9 + 1) * 1.0 / 10.0
            report_string = 'STATUS\t%f' % fraction_done
            pgui(report_string)
        
        
            cmd='rm -f ' + path + os.sep + 'dock_*.pdb'
            os.system(cmd)
            cmd='rm -f ' + path + os.sep + 'score*'
            os.system(cmd)
            cmd='rm -f ' + path + os.sep + '*cycle.dcd'
            os.system(cmd)
            cmd='rm -f ' + path + os.sep + '*dock*.log'
            os.system(cmd)
            # test for spactial voxel convergence
            outfile_name = path + os.sep + mvars.runname + '_dock.dcd'
            calc_spatial_convergence_all(mvars.refpdb, [outfile_name])
            ref = path + os.sep + 'dock'
            voxels = numpy.loadtxt(ref + '_occupied_voxels.npy')
            n_row, _, = voxels.shape
            sum = 0.0
            for i in xrange(n_row):
                #if i > n_row-2000:
                if i > n_row-8:
                     sum = sum + voxels[i,1] - voxels[i-1,1]
    
            #avg_voxel_diff = sum/2000.0
            #if mvars.maximum_dock_cycles == 20:
            avg_voxel_diff = sum/8.0
            log.info('avg_voxel_diff: %s' % avg_voxel_diff)
            if mvars.maximum_dock_cycles == 5:
                if avg_voxel_diff < mvars.spatial_voxel_convergence_tolerance:
                    log.info('End Docking Cycles--spatial_voxel_convergence_tolerance is satisfied')
                    break   # End docking cycles

        
        mvars.dcdfile_fullatom = path + os.sep + mvars.runname + '_dock.dcd'
        docking_results = os.path.join(mvars.runname, app)
        direxist = os.path.exists(docking_results)
        if(direxist == 0):
            try:
                result = os.system('mkdir -p ' + docking_results)
            except:
                message = 'can not create project directory: ' + docking_results
                message += '\nstopping here\n'
                print_failure(message, txtOutput)
            if(result != 0):
                message = 'can not create project directory: ' + docking_results
                message += '\nstopping here\n'
                print_failure(message, txtOutput)

        cmd = 'mv ' + path + os.sep + mvars.runname + '_dock.dcd ' + docking_results + os.sep + mvars.runname + '_dock.dcd'
        os.system(cmd)
        cmd = 'mv ' + mvars.refpdb + ' ' + docking_results + os.sep +  mvars.runname + '_dock.pdb'
        os.system(cmd)
        cmd = 'mv ' + ref + '_occupied_voxels.npy ' + docking_results + os.sep + mvars.runname + '_occupied_voxels.npy'
        os.system(cmd)
        return
    
    def sascalc_x2(self, variables, txtOutput):
        log = self.log
        from convergence_test import calc_chisq
        #USE SASSIE TO CALCULATE SAS PROFILES AND CHISQ FILTER
        mvars = self.mvars

        if (mvars.sas_option == 'neutron'):
            sasfolders = os.path.join(mvars.runname, 'neutron_D2Op_100')
        else:
            sasfolers = os.path.join(mvars.runname, 'xray')
        direxist = os.path.exists(sasfolders)
        log.info('sasfolders: ' + sasfolders)
        if(direxist == 0):
            try:
                result = os.system('mkdir -p ' + sasfolders)
            except:
                message = 'can not create project directory: ' + sasfolders
                message += '\nstopping here\n'
                print_failure(message, txtOutput)
            if(result != 0):
                message = 'can not create project directory: ' + sasfolders
                message += '\nstopping here\n'
                print_failure(message, txtOutput)

        chisq_scores = os.path.join(mvars.runname, app)
        log.info('chi_square scores: ' + chisq_scores)
        direxist = os.path.exists(chisq_scores)
        if(direxist == 0):
            try:
                result = os.system('mkdir -p ' + chisq_scores)
            except:
                message = 'can not create project directory: ' + chisq_scores
                message += '\nstopping here\n'
                print_failure(message, txtOutput)
            if(result != 0):
                message = 'can not create project directory: ' + chisq_scores
                message += '\nstopping here\n'
                print_failure(message, txtOutput)

        #print mvars.dcdfile_fullatom
        #print mvars.refpdb
        
        #SPACE FOR CALL SASCALC (mvars.dcdfile_fullatom, mvars.refpdb, mvars.i0, ........
        #SPACE FOR CALL CHI_SQUARE (mvars.dcdfile_fullatom, mvars.refpdb, interpolated_data_file, mvars.sas_option \
        #     mvars.i0, mvars.number_of_q_values,  mvars.maximum_q_value, ........

        #CAll calc_chisq to calculate chi-sq values
        x2outfile = chisq_scores + os.sep + 'x2file.txt'
        log.info('x2outfile : '+x2outfile)
        calc_chisq(sasfolders, mvars.interpolated_data_file, x2outfile)

        path = os.path.join(mvars.runname, 'work')
        mvars.scores_file = path + os.sep + 'docking_scores.txt'
        log.info('scores_file : '+mvars.scores_file)
        file = open(mvars.scores_file,"r").readlines()

        # Read LR docking scores
        score_list=[]
        for line in file:
            score = string.split(line)
            score_list.append(float(score[1]))

        #Read ChiSQ (X2) VALUES
        file = open(x2outfile,"r").readlines()
        x2_list=[]
        i = 0
        for line in file:
            #skip first 2 (header) lines
            i += 1
            x2 = string.split(line)
            if i > 2:
                x2_list.append(float(x2[1]))

        log.info('outfile : ' + chisq_scores + os.sep +  mvars.runname + '_scores_vs_x2.txt')
        outfile = open(chisq_scores + os.sep +  mvars.runname + '_scores_vs_x2.txt','w')
        #print('x2_list = ', x2_list)
        #print('score_list = ', score_list)
        #print('len(x2_list) = ',len(score_list))
        #print('len(score_list) = ',len(score_list))

        for i in range(len(score_list)):
        #    print('i = %i, x2_list[i] = %s, score_list[i] = %s\n' % (i, x2_list[i], score_list[i]))
            outfile.write('%s %s\n' %  ( x2_list[i], score_list[i]))

        return

    def epilogue(self):
        '''
        method to print out simulation results and to move results
        to appropriate places.
        '''

        log = self.log
        log.debug('in epilogue')
        pgui = self.run_utils.print_gui

        self.run_utils.clean_up(log)

        pgui('Configurations and statistics saved in %s directory\n\n' % (self.runpath+os.path.sep))

        pgui('DCD data were written to %s\n\n' % os.path.join(self.runpath, self.mvars.runname + '_dock.dcd'))

        pgui('Reference pdb file written to %s\n\n' % os.path.join(self.runpath, self.mvars.runname + '_dock.pdb'))

        pgui('Docking scores were written to %s\n\n' % os.path.join(self.runpath, self.mvars.runname + '_scores_vs_x2.txt'))

        pgui('Spatial voxels data were written to %s\n\n' % os.path.join(self.runpath, self.mvars.runname + '_occupied_voxels.npy'))

        pgui("\n" + "=" * 60 + " \n")
        time.sleep(0.1)

        return

