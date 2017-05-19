import sassie.sasmol.sasmol as sasmol
import numpy,sys,random

a = sasmol.SasMol(0)

a.read_pdb('simple.pdb')

segname = a.segname()

print 'segname = ',segname

print '\n> THIS IS THE FULL PDB (A)'
print 'a.coor = ',a.coor()[0]
print 'a.coor()[0] = ',a.coor()[0]
print 'a.coor()[0][0] = ',a.coor()[0][0]
print 'a.coor()[0][-1] = ',a.coor()[0][-1]

print 'a.coor().shape = ',a.coor().shape

basis1 = 'segname[i] == "ENDA"'
basis2 = 'segname[i] == "ENDB"'

error,mask1 = a.get_subset_mask(basis1)
error,mask2 = a.get_subset_mask(basis2)

b = sasmol.SasMol(1)
error = a.copy_molecule_using_mask(b,mask2,0)
if(error != []): print 'error = ',error
print '\n> THIS IS THE A COPY OF ENDB (B)'
print 'bcoor = ',b.coor()[0]
print 'b.coor().shape = ',b.coor().shape
ran=random.random()
print 'random pertubation = ',ran
b.coor()[0] = b.coor()[0]+ran

print '\n\n'

error = a.set_coor_using_mask(b,0,mask1)
if(error != []): print 'error = ',error

print '\n> THIS IS THE NEW FULL PDB (A) w/ ENDB + random'
print 'new acoor = ',a.coor()[0]
print 'a.coor().shape = ',a.coor().shape

print '\n\n'

sys.exit()




natoms = a.natoms()
indicies = indicies=numpy.nonzero(mask1*numpy.arange(1,natoms+1))[0]
print 'local_indicies = ',indicies

ncoords = a.coor()[0]
print 'ncoords = ',ncoords
print 'ncoords.shape = ',ncoords.shape
print 'ncoords.flatten() = ',ncoords.flatten()

three_indicies = []
for i in xrange(len(indicies)):
	this_index = indicies[i]*3
	three_indicies.append([this_index,this_index+1,this_index+2])
print '3I = ',numpy.array(three_indicies).flatten()
ti = numpy.array(three_indicies).flatten()

c = numpy.take(ncoords,ti)
c.shape = (-1,3)
print 'c = ',c

numpy.put(b.coor()[0],ti,c)

print 'new bcoor = ',b.coor()[0]
print 'b.segname() = ',b.segname()
print 'b.coor().shape = ',b.coor().shape


'''

#temp=numpy.take(ncoords,(0,1,2),-1)
xcoor = numpy.take(ncoords,(0,0),1)
print 'xcoor = ',xcoor  
ycoor = numpy.take(ncoords,(1,0),1)
print 'ycoor = ',ycoor 
zcoor = numpy.take(ncoords,(2,0),1)
print 'zcoor = ',zcoor
#temp=numpy.take(ncoords,([0,1,2],[0,1]))
#print 'temp = ',temp

qsel = numpy.arange(0,natoms)
totsnz1=numpy.zeros(len(temp)*3,numpy.int32)
lts=0
print 'len(temp) = ',len(temp)
print 'len(qsel) = ',len(qsel)
for ts in xrange(len(temp)):
	totsnz1[lts]=(qsel[ts]*3)
	totsnz1[lts+1]=(qsel[ts]*3)+1
	totsnz1[lts+2]=(qsel[ts]*3)+2
	lts=lts+3

numpy.put(b.coor()[0],totsnz1,temp)


print 'new bcoor = ',b.coor()[0]

'''
