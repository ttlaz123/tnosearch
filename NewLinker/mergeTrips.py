import sys
import os
tnopath = os.environ['TNO_PATH']
sys.path.insert(0, tnopath)
import LinkerLib as LL
from LinkerLib import Triplet
from LinkerLib import Detection
from LinkerLib import writeTriplets
from LinkerLib import pickleTriplets

from astropy.table import Table
import pickle
import time
import argparse

def newMergeTrips(triplets, savename, thresh=10):
    finalList = []
    totalSize = len(triplets)
    print('\nnumber to merge: '+ str(totalSize))
    checkList = {}
    rmList = {}
    trackid = 0
    time0 = time.time()
    trackList = {}
    objidTrack = {}
    tempTrack = 0
    for trip in triplets:
        trip.trackid = tempTrack
        trackList[tempTrack] = trip
        tempTrack += 1
        for det in trip.dets:
            if(det.objid in objidTrack):
                objidTrack[det.objid].append(trip.trackid)
            else:
                objidTrack[det.objid] = [trip.trackid]

    while(len(trackList)>0):
        LL.printPercentage(totalSize - len(trackList), totalSize, time.time()-time0)
        track, thisTrip = trackList.popitem()
        update = False
        for det in thisTrip.dets:
            for tempid in objidTrack[det.objid]:
                if(tempid not in trackList):
                    continue
                trip = trackList[tempid]
                if(thisTrip.shareM(trip, max(len(thisTrip.dets), len(trip.dets))/2 + 1)
                    or thisTrip.isSubset(trip)
                    or trip.isSubset(thisTrip)):
                    newDets = thisTrip.merge(trip)
                    checkList[trackid] = newDets
                    rmList[trackid] = (thisTrip, trip)
                    trackList.pop(tempid)
                    trackid += 1
                    update = True
                    break
            if(update):
                break
        if(not update):
            finalList.append(thisTrip)
    if(len(checkList)==0):
        return finalList
     
    goodMergeIds, badMergeIds = siftChecks(checkList, savename, thresh)
    for i, chi in goodMergeIds:
        trip = Triplet(checkList[i])
        trip.chiSq = chi
        finalList.append(trip)
    for i, chi in badMergeIds:
        trip1, trip2 = rmList[i]
        if(trip1.chiSq > 0 and trip2.chiSq > 0):
            if(trip1.chiSq/len(trip1.dets) > trip2.chiSq/len(trip2.dets)):
                finalList.append(trip2)
            else:
                finalList.append(trip1)
        else:
            if(trip1.realLength() > trip2.realLength()):
                finalList.append(trip1)
            else:
                finalList.append(trip2)
    return newMergeTrips(finalList, savename, thresh)

def getInitTrips(triplets, savename):
    trackID = 1
    trackDict = {}
    print('assigning trackids')
    trips = []
    for trip in triplets:
        trips.append(trip)
        trip.trackid = trackID
        trackDict[trackID] = trip.dets
        trackID += 1
    triplets = trips
    finalList = []  
    goodIDs, badIDs = siftChecks(trackDict, savename)
    for i, chi in goodIDs:
        trip = Triplet(trackDict[i])
        trip.chiSq = chi
        trip.trackid = i
        finalList.append(trip)
    for i, chi in badIDs:
        trip = Triplet(trackDict[i])
        trip.trackid = i
        trip.chiSq = chi
        print('bad triplet:')
        print(trip)
        finalList.append(trip)
    return finalList 
    
def siftChecks(trackDict, savename, thresh=30):
    outName = LL.writeDetToOrb(trackDict, 'siftRequestMerge+' + savename + '.fits')
    cmd = 'BulkFit -observationFile=' + str(outName) + ' -orbitFile='
    paramName = outName.split('.')[0] + '.orbit'
    cmd = cmd + paramName
    print('running BulkFit...')
    time0 = time.time()
    os.system(cmd)
    print('done after ' + str(time.time()-time0) + ' seconds')
    orbits = Table.read(paramName, format='fits')
    chiList = orbits['CHISQ'].tolist()
    orbitIDs = orbits['ORBITID'].tolist()
    dofList = orbits['DOF'].tolist()
    flagList = orbits['FLAGS'].tolist()
    goodid = []
    badid = [] 
    for x in range(len(chiList)):
        if(chiList[x] < thresh*dofList[x] and 
                    flagList[x] == 0 
                    and chiList[x] > 0.0):
            goodid.append((orbitIDs[x], chiList[x]))
        else:
            badid.append((orbitIDs[x], chiList[x]))
    return goodid, badid

def mergeFakes(triplets):
    fakeDict = {}
    for trip in triplets:
        fakeid = trip.dets[0].fakeid
        if(fakeid in fakeDict):
            fakeDict[fakeid].extend(trip.dets)
        else:
            fakeDict[trip.dets[0].fakeid] = trip.dets
    trips = []
    for trip in fakeDict.values():
        trips.append(Triplet(trip))
    return trips
    
def main():
    args = argparse.ArgumentParser()
    args.add_argument('triplets', nargs='+', help='list of triplets to merge')
    args.add_argument('-f', '--fake', action='store_true', 
            help='whether to look at fakes')
    args = args.parse_args()
    
    savename1 = args.triplets[0].split('/')[-1].split('.')[0]
    savename = 'merged+' + savename1 
    triplets = []
    
    print('loading triplets...')
    time0 = time.time()
    for trips in args.triplets:
        triplets.extend(pickle.load(open(trips, 'rb')))
    print('done loading after ' + str(time.time()-time0) + ' seconds')
    if(args.fake):
        mergedTrips = mergeFakes(triplets)
    else:
        siftedTrips = getInitTrips(triplets, savename1)
        mergedTrips = newMergeTrips(siftedTrips, savename1)
    print('\nsize of final list = ' + str(len(mergedTrips)))    
    finalList =[]
    for trip in mergedTrips:
        if(trip.realLength()>4):
            finalList.append(trip)
    print('size after reducin: ' + str(len(finalList)))
    writeTriplets(finalList, savename + '.txt')
    pickleTriplets(finalList, savename + '.pickle')

if __name__ == '__main__':
    main()
