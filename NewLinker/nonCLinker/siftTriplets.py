import sys
import os
tnopath = os.environ['TNO_PATH']
sys.path.insert(0, tnopath)
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
import ephem

from multiprocessing import Pool, cpu_count, Manager, Queue

#find the potential orbits given a list of tripletsi
'''
input: a list of triplets
output: goodList, a list of triplets that satisfy orbits
        badList, a list of triplets that don't satisfy orbits
        missList, a list of triplets that don't satisfy, but are P9 like or
            are the same fake for fakes
'''
def potentialTriplets(triplets, detDict, chiThresh=50):#, queue):
   
    #queue is 0 if no multiprocessing
    size = len(triplets)
    counter = 0
    time0 = time.time()
    goodList = []
    missList = []
    # for each detection, loop through each possible combination of triplets
    for trip in triplets:
        if(not isinstance(trip,Triplet)):
            dets = [detDict[x] for x in trip] 
            trip = Triplet(dets)
        if(not isinstance(trip.dets[0], Detection)):
            dets = [detDict[x] for x in trip.dets]
            trip.dets = dets
        if(trip.trackid == -1):
            trip.trackid = counter
       # if(queue == 0):
        if(counter%1000==0):
            printPercentage(counter, size, time.time()-time0)
        if(trip.magCheck):
            elements, errs = trip.calcOrbit()
            chisq = trip.getChiSq()
            #orbit is a swigpy object that can't be saved for some reason
            trip.orbit = 0
            if(chisq < chiThresh): #and elements['a'] > 2 and elements['e'] < 1 and chisq < 10):
            #if chisq < threshold:
                goodList.append(trip)
            else:
                # already missing quite a lot before this stage
                #if(trip.p9Like()):
                #    p9List.append(trip)
                if(trip.sameFake()):
                    missList.append(trip)
        counter += 1
    '''
    badid = []
    for bad in badList:
        for det in bad.dets:
            badid.append(det.fakeid)

    missid = [x.dets[0].fakeid for x in missList]
    goodid = [x.dets[0].fakeid for x in goodList]
    p9id = [x.dets[0].fakeid for x in p9List]
    '''
    #print('\n Number missing: ' + str(len(set(missid))))
    #print('Number good: ' + str(len(set(goodid))))
    #print('Number p9 like: ' + str(len(set(p9id))))
    #print('Number bad: ' + str(len(set(badid))))
    #print('Total: ' + str(len(set(badid+missid+goodid))))i
    #if(queue != 0):
    #    queue.put(allTrips)
    return goodList, missList#, p9List

#check triplets with multiproccessing
def multiProcessTriplets(triplets, Ncpu):
    print('sifting triplets')
    pool = Pool(Ncpu)
    manager = Manager()
    queue = manager.Queue()
    splitList = [([triplet],queue) for triplet in triplets]
    result = pool.map_async(potentialTriplets, splitList)
    #monitor loop
    size_old = 0
    time0 = time.time()
    while not result.ready():
        print('triplet filtering progress: ' + str(queue.qsize()) +
            ' of ' + str(len(triplets)) + ' triplets checked ' +
                str(float(queue.qsize())/len(triplets)*100) + '% at ' +
                str(time.time()-time0) + ' sec')
        time.sleep(1)
    
    #consolidate reults into a few lists
    goodTrips = []
    missTrips = []
    badTrips = []
    p9Trips = []
    allTrips = []
    for r in result.get():
        goodTrip, missTrip, badTrip, p9Trip, allTrip = r
        goodTrips.append(goodTrip)
        missTrips.append(missTrip)
        badTrips.append(badTrip)
        p9Trips.append(p9Trip)
        allTrips.append(allTrip)
    pool.close()
    goodTrips = [trip for triplets in goodTrips for trip in triplets]
    missTrips = [trip for triplets in missTrips for trip in triplets]
    badTrips = [trip for triplets in badTrips for trip in triplets]
    p9Trips = [trip for triplets in p9Trips for trip in triplets]
    allTrips = [trip for triplets in allTrips for trip in triplets]
    return goodTrips, missTrips, badTrips, p9Trips, allTrips

def main():
    args = argparse.ArgumentParser()
    args.add_argument('triplets',
                        help='path to pickle file; file has format ' + 
                        'chunk###+SNOBS_SEASON###_ML0#.pickle')
    #args.add_argument('chisq', nargs='?', default=5,
    #                    help='chisq threshold to be considered a good triplet')
    args.add_argument('detFile', help='path to csv file with detections')
    args.add_argument('-c', '--ncpus', help='number of cpus')
    args.add_argument('-p', '--printP', help='1 to print out percentage, 0 to not')
    args = args.parse_args()

    detDict = LL.objidDictionary(args.detFile)
    print('\nopen pickle file ' + args.triplets) 
    pickleFile = open(args.triplets, 'rb')
    print('load pickle file')
    triplets = pickle.load(pickleFile)
    print('done loading')
    Ncpu = 1
    if(args.ncpus):
        Ncpu = int(args.ncpus)

    #checks orbits for good fit
    goodTrips, missTrips = potentialTriplets(triplets, detDict)
    #multiProcessTriplets(triplets, Ncpu)


    #saving triplets
    saveName = args.triplets.split('+')[-1].split('.')[0]
    chunkName = args.triplets.split('/')[-1].split('+')[0]
    writeTriplets(goodTrips, chunkName + '+goodtriplets+' + saveName + '.txt', False)
    writeTriplets(missTrips, chunkName + '+misstriplets+' + saveName + '.txt', False)
    #writeTriplets(p9Trips, chunkName + '+p9triplets+' + saveName + '.txt', False)

    pickleTriplets(goodTrips, chunkName + '+goodtriplets+' + saveName + '.pickle')
    pickleTriplets(missTrips, chunkName + '+misstriplets+' + saveName + '.pickle')
    #pickleTriplets(p9Trips, chunkName + '+p9triplets+' + saveName + '.pickle')

if __name__=='__main__':
    main()
