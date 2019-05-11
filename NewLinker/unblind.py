import sys
import os
tnopath = os.environ['TNO_PATH']
sys.path.insert(0, tnopath)

import argparse
import pickle
from astropy.table import Table
import pandas as pd
from LinkerLib import Triplet
from LinkerLib import writeTriplets
from LinkerLib import pickleTriplets
import LinkerLib as LL


def makeFakeDict(fakesFile):
    print('reading from '+ fakesFile)
    tab = Table.read(fakesFile)
    fakeDict = {}
    fakeids = tab['ORBITID'].tolist()
    objids = tab['CATALOG_ID'].tolist()
    print('making dictionary')
    for x in range(len(fakeids)):
        fakeDict[objids[x]] = fakeids[x]
    print('done making dictionary')
    return fakeDict
    
    


def unblindTrips(triplets, fakeDict):
    print('unblinding triplets')
    for trip in triplets:
        for det in trip.dets:
            det.fakeid = fakeDict[det.objid]
    return triplets


def unblindPatch(fakeDict, patchNumber):
    patch = "{0:03d}".format(patchNumber)
    df = pd.read_csv('ccdRM+SNFIT_SEASON' + patch + '_ML06.csv')
    df = unblindDetections(df, fakeDict)
    df.to_csv('unblinded+SNFIT_SEASON' + patch + '_ML06.csv')


def unblindDetections(df, fakeDict):
    objids = df['objid'].tolist()
    fakeids = []
    for objid in objids:
        fakeids.append(fakeDict[objid])
    df['fakeid'] = pd.Series(fakeids, index=df.index)
    return df

def makeFakes(csvFile):
    seasonInfo = csvFile.split('+')[-1].split('.')[0]
    detList = LL.wrapDets(csvFile) 
    tripDict = {}
    for det in detList:
        if(det.fakeid == -1):
            continue
        if(det.fakeid in tripDict):
            tripDict[det.fakeid].dets.append(det)
        else:
            tripDict[det.fakeid] = Triplet([det])
    tripList = []
    for k in tripDict:
        tripList.append(tripDict[k])
    
    writeTriplets(tripList, 'trueTriplets+' + seasonInfo + '.txt')
    pickleTriplets(tripList, 'trueTriplets+' + seasonInfo + '.pickle')



def main():
    args = argparse.ArgumentParser()
    args.add_argument('-t', '--triplets', help='pickle file of list of triplets')
    args.add_argument('-f', '--fakes', help='file of fakeids and objids')
    args.add_argument('-c', '--detections', help='unblind detections of a patch')
    args.add_argument('-u', '--truths', help='csvfile with fakeids')
    args = args.parse_args()
    
    if(args.truths):
        makeFakes(args.truths)
        return

    if(args.detections):
        fakeDict = makeFakeDict(args.fakes)
        for patch in range(32, 100):
            print('unblinding patch: ' + str(patch))
            unblindPatch(fakeDict, patch)
        return


    base = args.triplets.split('.')[0]
    if('/' in base):
        base = base.split('/')[-1]

    triplets = []
    
    print('loading from: ' + args.triplets)
    with open(args.triplets, 'rb') as f:
        triplets = pickle.load(f)

    print('number of triplets:' + str(len(triplets)))
    
    fakeDict = makeFakeDict(args.fakes)

    unblindedTrips = unblindTrips(triplets, fakeDict)

    seasonInfo = base.split('+')[-1]
    writeTriplets(unblindedTrips, 'unblindedTriplets+' + seasonInfo + '.txt')
    pickleTriplets(unblindedTrips, 'unblindedTriplets+' + seasonInfo + '.pickle')

if __name__== '__main__':
    main()
    

