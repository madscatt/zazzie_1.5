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
import os,sys,locale,string
import input_filter

def check_minimize(variables,**kwargs):
      
    runname     = variables['runname'][0]
    infile		= variables['infile'][0]
    pdbfile		= variables['pdbfile'][0]
    outfile		= variables['outfile'][0]
    nsteps		= variables['nsteps'][0]
    parmfile	= variables['parmfile'][0]
    psffile		= variables['psffile'][0]
    ncpu		= variables['ncpu'][0]
    keepout		= variables['keepout'][0]
    dcdfreq		= variables['dcdfreq'][0]
    infiletype	= variables['infiletype'][0]

    md	    	= variables['md'][0]
    mdsteps		= variables['mdsteps'][0]
    dielect		= variables['dielect'][0]
    temperature	= variables['temperature'][0]

    error=[]
    error = input_filter.check_name(runname)
    if(error!=[]):
        return error

    if "no_file_check" not in kwargs:
        path='./'

        ev,rv,wv=input_filter.check_permissions(path)
        if(not ev or not rv or not wv):
            error.append('permission error in input file path '+path+'  [code = '+str(ev)+str(rv)+str(wv)+']')
            if(ev==False):
                error.append('path does not exist')
            elif(rv==False):
                error.append('read permission not allowed')
            elif(wv==False):
                error.append('write permission not allowed')
            return error

    if(ncpu < 1):
        error.append( 'ncpu needs to be >= 1, ncpu = '+str(ncpu))
        return error
    if(nsteps < dcdfreq):
        error.append( 'max steps < dcd write frequency, nsteps = '+str(nsteps))
        return error
    elif(keepout != 0 and keepout != 1):
        error.append( 'keepout == 0 for "no" and 1 for "yes", keepout = '+str(keepout))
        return error
    elif(nsteps < 20):
        error.append( 'max steps < 20 : max steps = '+str(nsteps))
        return error

    
    parmfiles = string.split(parmfile,',')
    for this_parmfile in parmfiles:
        error=input_filter.check_file_exists(this_parmfile)
        if(len(error) != 0):
            error.append('check parameter path and filename: '+str(this_parmfile)+' JEC')
            error.append('check parameter path and filename: '+str(parmfile)+' JEC')
            return error

    error=input_filter.check_file_exists(psffile)

    if(len(error) != 0):
        error.append('\ncan not find psf file in "./"\n')
        return error

          ### now check and see if input file is a correct pdb or dcd and compare to psf file    
	
    error=input_filter.check_file_exists(infile)
    if(len(error) != 0):
        error.append('check input filename (dcd)  path + filename : '+infile)
        return error

    error=input_filter.check_file_exists(pdbfile)
    if(len(error) != 0):
        error.append('check input filename (pdb)  path + filename : '+pdbfile)
        return error

    ev,value=input_filter.check_pdb_dcd(infile,'pdb')
    if(ev == 1): # file exists
        if(value == 0):         # it is not a pdb file
            ev,value=input_filter.check_pdb_dcd(infile,'dcd')
            if(value == 1):
                cvalue=input_filter.certify_dcd_psf(infile,psffile)
                if(cvalue == 0):
                    error.append('input file '+infile+' and psf file '+psffile+' are not compatible')
                    return error
            else:
                error.append('input file '+infile+' is not a valid pdb or dcd file')
                return error
            infiletype='dcd'
        else:                   # is a pdb file
            ev,cvalue=input_filter.certify_pdb_psf(infile,psffile)
            if(cvalue == 0):
                error.append('input file '+infile+' and psf file '+psffile+' are not compatible')
                return error
            infiletype='pdb' 
    else:
        error.append('input (pdb or dcd) file '+infile+' does not exist')
        return error

    ### now check and see if input pdb file (reference) is correct and matches pdb/dcd file 

    ev,value=input_filter.check_pdb_dcd(pdbfile,'pdb')
    if(ev == 1): # file exists
        if(value == 0):         # it is not a pdb file
            error.append('input pdb file'+pdbfile+' is not valid')
            return error
        elif(infiletype=='dcd'):
            cvalue=input_filter.certify_pdb_dcd(pdbfile,infile)
            if(cvalue == 0):
                error.append('input pdb file '+pdbfile+' and input file '+infile+' are not compatible')
                return error
                        
        elif(infiletype=='pdb'):
            locvariables=['name']
            value,result1=input_filter.get_pdb_stats(infile,locvariables)
            value,result2=input_filter.get_pdb_stats(pdbfile,locvariables)
            if(result1 != result2):
                error.append('reference pdb file '+pdbfile+' and input pdb file '+infile+' are not compatible')
                return error
    else:
        error.append('input (pdb or dcd) file '+pdbfile+' does not exist')
        return error

    variables['infiletype']=(infiletype,'string')

    if(md==1 or md==2):
        if(mdsteps <= 1):
            error.append('number of md steps needs to be greater than zero : '+str(mdsteps))
            return error
        elif(mdsteps < 20):
            error.append('number of md steps needs to be multiple of twenty : '+str(mdsteps))
            return error
        elif(mdsteps%2 != 0):
            error.append('number of md steps needs to be divisible by two : '+str(mdsteps))
            return error
        elif(dielect < 0.0):
            error.append('dielectric needs to be greater than zero : '+str(dielect))
            return error
        elif(temperature < 0.0):
            error.append('temperature needs to be greater than zero : '+str(temperature))
            return error

    return error



