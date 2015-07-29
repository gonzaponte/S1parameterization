from ROOT import *
import os
import sys

gROOT.Load('libirene')

def PickData( evt, hit_map ):
    for hit in evt.GetSensorHits():
        if hit.GetDetectorName() == 'PmtR11410':
            hit_map[ hit.GetID() ] += hit.GetAmplitude()

def check_files( files ):
    missing = filter( lambda f: not os.path.exists( f ), files )
    if not missing: return files

    print 'The following files are missing:'
    for i,m in enumerate(missing): print i,'->',m
    ask = raw_input( 'Do you want to carry on with the remaining files? (y/n) ' )
    if 'y' in ask:
        return filter( os.path.exists, files )
    else:
        print 'Leaving the program...'
        sys.exit()

PMT_IDs = range(12)
inputfolder  = '/data4/NEXT/users/gmartinez/topo/jobs/nexus/S1/data/'
outputfolder = '/data4/NEXT/users/gmartinez/topo/jobs/nexus/S1/data/'

outputstr = 'PointID ' + ' '.join( 'PMT' + str(id) for id in PMT_IDs ) + '\n'

npoints  = 118096
nphotons = 1e6
nphotons = 1.0/nphotons

filenames = [ inputfolder + 'S1_{0}.next'.format(i) for i in range(npoints) ]
filenames = check_files( filenames )

evt = irene.Event()
for pointID, filename in enumerate(filenames):
    if not pointID % 10000: print 'point number', pointID
    datafile = TFile( filename )
    datatree = datafile.Get('EVENT')
    datatree.SetBranchAddress('EventBranch',evt)
    hit_map = { id : 0 for id in PMT_IDs }
    for j in range( datatree.GetEntries() ):
        datatree.GetEntry(j)
        PickData( evt, hit_map )
    for id in PMT_IDs: hit_map[id] *= nphotons
    outputstr += ' '.join( map( str, (pointID,) + zip(*sorted(hit_map.items()))[1] ) ) + '\n'

open( outputfolder + 'S1table.dat', 'w' ).write( outputstr )
