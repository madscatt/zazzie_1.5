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
import os
import string
import sys

'''

This module is called by DOCKING.PY

Writes input files for Rosseta Prepack and Docking.

'''


def write_flags(pdbfile, partner1_partner2, path):
    '''
    WRITE_FLAGS is the function to write a generic input
    for Rosetta Docking Program.


    INPUT:  variable descriptions:

    pdbfile:                 name of pdbfile
    partner1_partner2:       protein chain(s) comprising 1st partner (fixed) and second partner (moved)
    OUTPUT:
    flage_pre:               file with input options for Rosetta "prepacking"
    flage:                   file with input options for Rosetta "full docking protocol"
    ignore_list              file for use with Rosetta programs to reduce output by excluding lines from output

    '''
    # flags_prepack

    outfile = open(path+ os.sep +'flags_pre', 'w')

    outfile.write('%s\n' % ('-s ' + pdbfile))
    outfile.write('%s\n' % ('-ex1'))
    outfile.write('%s\n' % ('-ex2aro'))
    outfile.write('%s\n' % ('-docking:partners ' + partner1_partner2))
    outfile.write('%s\n' % ('-mute core.util.prof ## dont show timing info'))
    outfile.write('%s\n' % ('-out:overwrite'))
    outfile.write('%s\n' % ('-out:file:fullatom #output in fullatom scorefile'))
    outfile.write('%s\n' % ('-out:file:scorefile prepack_score'))
    outfile.write('%s\n' % ('-mute core.io.database'))
    outfile.write('%s\n' % ('-unboundrot ' + pdbfile))
    outfile.write('%s\n' % ('-ignore_zero_occupancy false'))
    
    outfile.close()
    
    # flags_docking
    
    outfile = open(path+ os.sep +'flags', 'w')
    
    outfile.write('%s\n' % ('-ex1'))
    outfile.write('%s\n' % ('-ex2aro'))
    outfile.write('%s\n' % ('-native ' + pdbfile))
    outfile.write('%s\n' % ('-docking:partners ' + partner1_partner2))
    outfile.write('%s\n' % ('-randomize1'))
    outfile.write('%s\n' % ('-randomize2'))
    outfile.write('%s\n' % ('-dock_pert 3 8'))
    outfile.write('%s\n' % ('-mute core.util.prof ## dont show timing info'))
    outfile.write('%s\n' % ('-mute core.io.database'))
    outfile.write('%s\n' % ('-run:constant_seed'))
    outfile.write('%s\n' % ('-use_input_sc'))
    outfile.write('%s\n' % ('-unboundrot ' + pdbfile))
    outfile.write('%s\n' % ('-ignore_zero_occupancy false'))
    outfile.write('%s\n' % ('-spin'))
    
    outfile.close()

    # ignore_list

    outfile = open(path+ os.sep +'ignore_list', 'w')

    outfile.write('%s\n' % ('\/source\/bin\/'))
    outfile.write('%s\n' % ('\/database\/'))
    outfile.write('%s\n' % ('read.*from'))
    outfile.write('%s\n' % ('Reading\ in.*lib'))
    outfile.write('%s\n' % ('Basic\ usage\:'))
    outfile.write('%s\n' % ('Reading.*library'))
    outfile.write('%s\n' % ('Read.*lines.*from.*file'))
    outfile.write('%s\n' % ('Reading.*file'))
    outfile.write('%s\n' % ('INSERT\ INTO.*database'))
    outfile.write('%s\n' % ('Finished.+in [0-9]+ seconds.'))
    outfile.write('%s\n' % ('time'))
    outfile.write('%s\n' % ('Time'))
    outfile.write('%s\n' % ('seconds'))
    outfile.write('%s\n' % ('library took .+ seconds to load'))
    outfile.write('%s\n' % ('TIMING'))
    outfile.write('%s\n' % ('Rosetta version'))
    outfile.write('%s\n' % ('svn_version'))
    outfile.write('%s\n' % ('TIME_STAMP'))
    outfile.write('%s\n' % ('Warning: Unable to locate database file .*Dunbrack[0-9]+.lib.bin'))
    outfile.write('%s\n' % ('current_directory='))
    outfile.write('%s\n' % ('reating tricubic splines'))
    #outfile.write('%s\n' % ('core'))
    #outfile.write('%s\n' % ('protocols'))
    #outfile.write('%s\n' % ('basic'))

    outfile.close()

    return

if __name__ == '__main__':
    print 'running as a main process'

else:
    print 'running as a spawned process'
