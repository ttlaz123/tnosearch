import sys
import os
LINKERDIR=os.environ["TNO_PATH"]
sys.path.insert(0, LINKERDIR)
import time

from LinkerLib import printPercentage
from LinkerLib import Triplet
from LinkerLib import Detection
from LinkerLib import writeTriplets
from LinkerLib import pickleTriplets
import LinkerLib as LL

import numpy as np
import pickle
import argparse
import time
import subprocess

from astropy.io import fits
from astropy.table import Table
'''
input: --a list of triplets with objids as dets
       --a dictionary from objid to Detection objects
       --chunk number of the chunk
output: --name of the output file
'''
def writeProcessingFile(triplets, detDict, chunkname, savename):
    trackIDs = {}
    tripList = []
    for x in range(len(triplets)):
        trip = triplets[x]
        #print(trip)
        if(isinstance(trip, Triplet)):
            trip.trackid = x
            tripdets = trip.dets
            tripList.append(trip)
            nex = tripList[x].trackid
            #print('trackid: ' + str(nex))
            try:
                print(tripList[469])
            except:
                pass
        else:
            tripdets = trip[1]
 
        trackID = x
        dets = []
        try:
            print(tripList[469])
            print('test2^^')
        except:
            pass
        for objid in tripdets:
            if isinstance(objid, int):
                dets.append(detDict[objid])
            elif isinstance(objid, Detection):
                dets.append(objid)
            else:
                print(str(objid) + ' is not an int or Detection')
        trackIDs[trackID] = dets
        try:
            print(tripList[469])
            print('test3^^')
            if(tripList[469].trackid != 469):
                exit()
        except:
            pass
        
    for x in range(len(tripList)):
        trip = tripList[x]
        if(x != trip.trackid):
            print(trip.trackid)
            print('should be1: ' + str(x))
    
    outName = LL.writeDetToOrb(trackIDs, 
            'siftRequest+' + chunkname + '+' + savename + '.fits')
    
    return outName, trackIDs, tripList 

'''
input: --name of fits file with a list of triplets to be fit
output: --name of file to open for fits file

'''
#TODO turn this into an actual wrapper rather than just calling it with subprocess
def callBulkFit(tripletFile):
    cmd = 'BulkFit -observationFile=' + str(tripletFile) + ' -orbitFile='
    paramName = tripletFile.split('+',1)[-1].split('.')[0]+'.orbit'
    cmd = cmd + paramName
    cmd = cmd +' -bindingFactor=4'
    print('running BulkFit...')
    time0 = time.time()
    os.system(cmd)
    print('done running BulkFit after ' + str(time.time()-time0) + ' seconds')

    return paramName


'''
input: --name of orbitfile
output: writes orbital parameters in aei to orbtifile
'''

def callBulkEle(orbitfile):
    cmd = 'BulkElements -orbitFile=' + str(orbitfile)
    print('running BulkElements...')
    time0 = time.time()
    os.system(cmd)
    print('done running BulkElements after ' + str(time.time()-time0) + ' seconds')
    return


'''
input: --a fits filename with trackIDs and orbital parameters
       --a dictionary from trackIDs to objids
       --a chisq threshold
output: a list of triplets with their objids
'''
def siftTrips(fitsName, trackDict, chiThresh=50):
    print('only keeping orbtis with chisq under ' + str(chiThresh))
    orbits = Table.read(fitsName, format='fits')
    
    # version 1
    time1 = time.time()
    mask_chisq = (orbits['CHISQ']<chiThresh) & (orbits['FLAGS']==0)
    orbitIDs = orbits[mask_chisq]['ORBITID']
    tripList = []
    for row in orbitIDs:
        trackid = row 
        objids = trackDict[trackid]
#        print(trackid)
        del trackDict[trackid]
        trip = Triplet(objids)
        trip.trackid = trackid
        tripList.append(trip)
    
    time2 = time.time()
    '''
    #version 2
    tripList2 = []
    for row in orbits:
        if row['CHISQ'] < chiThresh:
            trackid = '{0:06d}'.format(row['ORBITID'])
            objids = trackDict[trackid]
            tripList2.append(Triplet(objids))
    time3 = time.time()
    '''
    print('done after ' + str(time2-time1) + ' seconds')
    
    #print('v2: ' + str(time3-time2))
   
    return tripList

def main():
    args = argparse.ArgumentParser()
    args.add_argument('triplets',
                        help='path to pickle file; file has format ' + 
                        'chunk###+SNOBS_SEASON###_ML0#.pickle')
    args.add_argument('detFile', help='path to csv file with detections')
    args.add_argument('-o', '--orbit', action='store_true', help='produce orbitfile only')
    args = args.parse_args()
    
    detDict = LL.objidDictionary(args.detFile)
    print('\nopen pickle file ' + args.triplets) 
    with open(args.triplets, 'rb') as pickleFile:
        print('load pickle file')
        triplets = pickle.load(pickleFile)
        print('done loading')
    saveName = args.triplets.split('+')[-1].split('.')[0]
    chunkName = args.triplets.split('/')[-1].split('+')[0]
    
    triplets = list(set(triplets))
    outname, trackIDs, lets= writeProcessingFile(triplets, detDict, chunkName, saveName)
 #   print(trackIDs.keys())
    fitsFile = callBulkFit(outname)
    if(args.orbit):
        callBulkEle(fitsFile)
    else:
        lets = siftTrips(fitsFile, trackIDs)
    #saving triplets
    writeTriplets(lets, chunkName + '+goodtriplets+' + saveName + '.txt', args.orbit)
    pickleTriplets(lets, chunkName + '+goodtriplets+' + saveName + '.pickle')

if __name__=='__main__':
    main()
