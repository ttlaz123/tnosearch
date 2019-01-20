import sys
import os
tnopath = os.environ['TNO_PATH']
sys.path.insert(0, tnopath)
import pickle
import time
from LinkerLib import Detection
from LinkerLib import Triplet
import numpy as np
import argparse
from LinkerLib import writeTriplets
from LinkerLib import pickleTriplets
from LinkerLib import Fake

def countTriplets(triplets):
    fakeids = []
    falsePositives = []
    fourPlus = []
    counterN = 0
    counterL = 0
    counterO = 0
    for trip in triplets:
        num = len(trip.dets)
        dof = (num-3)*2
        if(num>=4 and trip.chiSq/dof < 100):
            fourPlus.append(trip)
        
        if(trip.sameFake()):
            fakeids.append(trip.dets[0].fakeid)
            if trip.elements['a'] < 0:
                counterN += 1
            elif trip.elements['a'] > 200:
                counterL += 1
            else:
                counterO += 1      
            
        else:
            falsePositives.append(trip)
    print('total positive triplets : ' + str(len(fakeids)))
    print('total found orbits:' + str(np.unique(fakeids).size))
    print('total false positives: ' + str(len(falsePositives)))
    return falsePositives, counterN, counterL,counterO, fourPlus

def main():
    args = argparse.ArgumentParser()
    args.add_argument('tripletFile', nargs=1,
                        help='path to pickle file of triplets')
    args = args.parse_args()
    triplets = pickle.load(open(args.tripletFile[0], 'rb'))
    saveName = args.tripletFile[0].split('+')[-1].split('.')[0]
    falsePositives, x,y,z, fourPlus = countTriplets(triplets)
    print('negative a: ' + str(x))
    print('large a: ' + str(y))
    print('neither: ' + str(z))
    #writeTriplets(fourPlus, 'fourPlus+' + saveName + '.txt', True)
    #pickleTriplets(fourPlus, 'fourPlus+' + saveName + '.pickle')
if __name__ == '__main__':
    main()
