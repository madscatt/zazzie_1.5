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
import input_filter

def check_em_to_sans(variables,**kwargs):

    runname		= variables['runname'][0]           
    emfiletype	= variables['emfiletype'][0]        
    inputpath	= variables['inputpath'][0]         
    emdensityfile	= variables['emdensityfile'][0]     
    angpix		= variables['angpix'][0]            
    pdbfile		= variables['pdbfile'][0]           
    threshold	= variables['threshold'][0]         
    sansfile	= variables['sansfile'][0]          
    npoints		= variables['npoints'][0]           
    qmax		= variables['qmax'][0]              

    error=[]
    error = input_filter.check_name(runname)
    if(error!=[]):
        return error

    if 'no_file_check' not in kwargs:
        emdensityfile = inputpath+'/'+emdensityfile
        ev,rv,wv=input_filter.check_permissions(inputpath)

        if(not ev or not rv or not wv):
            error.append('permission error in input file path '+inputpath+'  [code = '+str(ev)+str(rv)+str(wv)+']')
            if(ev==False):
                error.append('path does not exist')
            elif(rv==False):
                error.append('read permission not allowed')
            elif(wv==False):
                error.append('write permission not allowed')
            return error
    
    if(emfiletype != 0 and emfiletype != 1):
        error.append('emfiletype == 0 for "cube file" and 1 for "mrc file", emfiletype = '+str(emfiletype))
        return error
    elif(npoints < 1):
        error.append('npoints needs to be > 1, npoints = '+str(npoints))
        return error
    elif(qmax < 0.0):
        error.append('qmax needs to be > 0, qmax = '+str(qmax))
        return error
    elif(len(angpix)!=3):
        error.append('angpix needs to have three values, npoints = '+str(angpix))
        return error
    for i in xrange(len(angpix)):
        tval=angpix[i]
        if tval < 0.0 :
            error.append('angpix needs to be > 0, angpix = '+str(angpix))
            return error
        
    error=input_filter.check_file_exists(emdensityfile)
    if(len(error) != 0):
        error.append('check input emdensity file path + filename')
        return error

    return error

