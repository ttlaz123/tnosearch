import sys
import os
tnopath = os.environ['TNO_PATH']
sys.path.insert(0, tnopath)

import argparse
import pickle
from astropy.table import Table

from LinkerLib import writeTriplets
from LinkerLib import pickleTriplets
import LinkerLib


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

def main():
    args = argparse.ArgumentParser()
    args.add_argument('-t', '--triplets', help='pickle file of list of triplets')
    args.add_argument('-f', '--fakes', help='file of fakeids and objids')
    args = args.parse_args()

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
    

