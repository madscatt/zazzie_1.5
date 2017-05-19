import sys
import sassie.simulate.sas_dock.sas_dock as sas_dock
import sassie.interface.input_filter as input_filter
import sassie.interface.sas_dock_filter as sas_dock_filter
import multiprocessing

svariables = {}

# BEGIN USER EDIT
# BEGIN USER EDIT
# BEGIN USER EDIT

runname = 'run_0'
pdbfile = '1ahw.pdb'
number_processors = '4'
maximum_dock_cycles = '1'
partner1_partner2 = 'AB_C'
sas_option = 'neutron'
#spatial_voxel_convergence_tolerance = '0.005'
interpolated_data_file = '1ahw_iq.dat'
spatial_voxel_convergence_tolerance = '100.0'
number_of_q_values = '13'
maximum_q_value = '0.12'
i0 = '1.000'

# END USER EDIT
# END USER EDIT
# END USER EDIT

svariables['runname'] = (runname, 'string')
svariables['pdbfile'] = (pdbfile, 'string')
svariables['number_processors'] = (number_processors, 'int')
svariables['maximum_dock_cycles'] = (maximum_dock_cycles, 'int')
svariables['partner1_partner2'] = (partner1_partner2, 'string')
svariables['sas_option'] = (sas_option, 'string')
svariables['interpolated_data_file'] = (interpolated_data_file, 'string')
svariables['spatial_voxel_convergence_tolerance'] = (spatial_voxel_convergence_tolerance, 'float')
svariables['number_of_q_values'] = (number_of_q_values, 'int')
svariables['maximum_q_value'] = (maximum_q_value, 'float')
svariables['i0'] = (i0, 'float')


error, variables = input_filter.type_check_and_convert(svariables)
if(len(error) > 0):
    print 'error = ', error
    sys.exit()

error = sas_dock_filter.check_docking(variables)
if(len(error) > 0):
    print 'error = ', error
    sys.exit()

txtQueue = multiprocessing.JoinableQueue()

simulation = sas_dock.simulation()
simulation.main(variables, txtQueue)

# the following line is for the 1.0 execution only
this_text = txtQueue.get(True, timeout=0.1)

