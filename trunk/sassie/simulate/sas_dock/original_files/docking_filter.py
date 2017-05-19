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
import os,sys
import sassie.interface.input_filter as input_filter
import sasmol.sasmol as sasmol

def check_docking(variables,**kwargs):

    runname = variables['runname'][0]
    pdbfile = variables['pdbfile'][0]
    number_processors = variables['number_processors'][0]
    maximum_dock_cycles = variables['maximum_dock_cycles'][0]
    partner1_partner2 = variables['partner1_partner2'][0]
    sas_option = variables['sas_option'][0]
    interpolated_data_file = variables['interpolated_data_file'][0]
    spatial_voxel_convergence_tolerance = variables['spatial_voxel_convergence_tolerance'][0]
    number_of_q_values = variables['number_of_q_values'][0]
    maximum_q_value = variables['maximum_q_value'][0]
    i0 = variables['i0'][0]




#    xon = variables['xon'][0]
#    if xon not in ['neutron','xray','neutron_and_xray']:
#        error.append("Either neutron or xray input need to be checked")
#        return error


    error=[]
    error = input_filter.check_name(runname)

    if(error!=[]):
        return error

    error=input_filter.check_file_exists(interpolated_data_file)
    if(len(error) != 0):
        error.append('input file, '+interpolated_data_file+', does not exist')
        return error

    error=input_filter.check_file_exists(pdbfile)
    if(len(error) != 0):
        error.append('input pdb file, '+pdbfile+', does not exist')
        return error
    ev,value=input_filter.check_pdb_dcd(pdbfile,'pdb')
    if(ev == 0):
        error.append('check input pdb file: '+pdbfile)
        return error
    if(value == 0):
        error.append( 'input pdb file, '+pdbfile+', is not a valid pdb file')
        return error
    try:
        m1 = sasmol.SasMol(0)
        m1.read_pdb(pdbfile)
        number_of_frames = m1.number_of_frames()
        print '> found '+str(number_of_frames)+' frames in PDB file'
        chains = m1.chains()
        underscore = 0
        partners = partner1_partner2.strip()
        for i in xrange(len(partners)):
            if(partners[i] == '_'):
                underscore = underscore + 1
        if(underscore != 1 or partners[0] == '_' or partners[len(partners)-1] == '_'):
            error.append('incorrect label usage in partner1_partner2 variable--incorrect or misplaced "_" symbol')
        partner_chains = []     
        for i in xrange(len(partners)):
            if(partners[i] != '_'): 
                   partner_chains.append(partners[i])
        if(partner_chains != chains):
            error.append('incorrect order of chains and "_" symbol in partner1_partner2 variable relative to chain order in input pdb file')
    except:
        error.append('could not open PDB file '+pdbfile+' to check number of frames')
        return error
    if(number_of_frames < 1):
        error.append('PDB file has no frames : '+pdbfile)
        return error

    if number_processors <= 0:
        error.append('number of processors needs to be greater than zero')
        return error
    elif spatial_voxel_convergence_tolerance <= 0:        
        error.append('spatial voxels convergence tolerance needs to be greater than zero')
        return error
    elif maximum_dock_cycles <= 0:        
        error.append('maximum number of docking cycles needs to be greater than zero')
        return error
    elif number_of_q_values <= 0:        
        error.append('maximum number of q values needs to be greater than zero')
        return error
    elif maximum_q_value <= 0:        
        error.append('maximum q value needs to be greater than zero')
        return error
    elif i0 <= 0:        
        error.append('I(0) needs to be greater than zero')
        return error

    # advanced options

    # check ProDy command

    # check Rg cutoffs

    '''
    elif(lowrg > highrg):
        error.append( 'low Rg cutoff is larger than high Rg cutoff, lowrg = '+str(lowrg)+' highrg = '+str(highrg))
        return error
    elif(lowrg < 0 or highrg < 0):
        error.append( 'Rg cutoffs need to be >= zero, lowrg = '+str(lowrg)+' highrg = '+str(highrg))
        return error

    elif(zflag != 0 and zflag != 1):
        error.append( 'ERROR in Z coordinate filter selection: zflag == 0 for "no" and 1 for "yes", zflag = '+str(zflag))
        return error
    elif(cflag != 0 and cflag != 1):
        error.append( 'ERROR in atomic constraints selection: cflag == 0 for "no" and 1 for "yes", cflag = '+str(cflag))
        return error
    elif(cflag == 1):
        err = input_filter.check_file_exists(confile)
        if(err != []):
            lerr=['ERROR in constraint filename selection: ']
            lerr.append(err)
            error.append(lerr[0]+err[0])
            return error
        filter_flag = 1
        m1=sasmol.SasMol(0)
        m1.read_pdb(pdbfile)
        err = constraints.read_constraints(m1,confile,filter_flag)
        if(err != []):
            error.append(err[0])
            return error

    '''


    return error
