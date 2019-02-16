import sys
import os
tnopath = os.environ['TNO_PATH']
sys.path.insert(0, tnopath)
from LinkerLib import Triplet
from LinkerLib import Detection
from LinkerLib import writeTriplets
from LinkerLib import pickleTriplets
from LinkerLib import printPercentage
import LinkerLib as LL

import time
import pickle
import argparse

def mergeChunks(folder, prefix1='small', prefix2='large'):
    
    normalTrips = []
    largeTrips = []
    print('listing files')
    files = os.listdir(folder)
    savename = files[0].split('+',2)[-1].split('.')[0]
    print('savename: ' + str(savename))
    counter = 0
    size = len(files)
    time0 = time.time()
    for f in files:
        printPercentage(counter, size, time.time()-time0)
        counter += 1
        if(f.endswith('.pickle') and prefix2 in f):
            with open(folder+f, 'rb') as fil:
                trips = pickle.load(fil)
                largeTrips.extend(trips)
        elif(f.endswith('.pickle') and prefix1 in f):
            with open(folder+f, 'rb') as fil:
                trips = pickle.load(fil)
                normalTrips.extend(trips)

    LL.pickleTriplets(normalTrips, 'chunkMerged+' + prefix1 + '+' + savename +'.pickle')
    LL.writeTriplets(normalTrips, 'chunkMerged+' + prefix1 + '+' + savename +'.txt')
    if(len(largeTrips) > 0):
        LL.pickleTriplets(largeTrips, 'chunkMerged+' + prefix2 + '+' + savename+'.pickle')
        LL.writeTriplets(largeTrips, 'chunkMerged+' + prefix2 + '+' + savename + '.txt')

def mergeFiles(folder, detDict):
    files = os.listdir(folder)
    trips = []
    savename = files[0].split('+')[-1].split('.')[0]
    counter = 0
    size = len(files)
    time0 = time.time()
    for f in files:
        LL.printPercentage(counter, size, time.time()-time0)
        counter += 1
        if(f.endswith('.pickle')):
            with open(folder+f, 'rb') as name:
                try:
                    triplets = pickle.load(name)
                    for trip in triplets:
                        for ind, det in enumerate(trip.dets):
                            if(not isinstance(det, Detection)):
                                det = detDict[det.objid]
                                trip.dets[ind] = det
                    trips+=(triplets)
                except EOFError:
                    print("file error: " + str(name))
                
                '''
                for trip in triplets:
                    if(trip.realLength()>3):
                        trips.append(trip)
                '''
    print('pickling size:' + str(len(trips)))
    LL.pickleTriplets(trips, 'chunkCombined+'+savename+'.pickle')
    LL.writeTriplets(trips, 'chunkCombined+'+savename+'.txt')

def main():
    args = argparse.ArgumentParser()
    args.add_argument('-f', '--folder', help='folder of chunks')
    args.add_argument('-d', '--detections', help='file with detections')
    args = args.parse_args()
    
    detDict = LL.objidDictionary(args.detections)

    if(args.folder):
        print('running on ' + args.folder)
        mergeFiles(args.folder, detDict)

if __name__=='__main__':
    main()
