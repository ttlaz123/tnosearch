import sys
import os
tnopath = os.environ['TNO_PATH']
sys.path.insert(0, tnopath)
import argparse
import pickle
import time

from LinkerLib import Detection
from LinkerLib import Triplet
from LinkerLib import writeTriplets
from LinkerLib import pickleTriplets
from LinkerLib import printPercentage
import LinkerLib as LL


def testCands(trip, tests):
    #print(trip)
    dets = trip.merge(tests)
    testTrip = Triplet(dets)
    testTrip.getChiSq()
    #print(testTrip)
    return testTrip


def testCand(trip, det):
    #print(trip)
    dets = trip.dets[:]
    dets.append(det)
    testTrip = Triplet(dets)
    testTrip.getChiSq()
    #print(testTrip)
    return testTrip.getChiSq()

def chooseMin1(triplet, cands, chiThresh):
    print('start choosemin1')
    minChi = 1000
    bestDets = []
    a_before = triplet.calcOrbit()[0]['a']
    for x in range(len(cands)):
        testTrip = testCands(triplet, [cands[x]])
        chi = testTrip.getChiSq()
        if(testTrip.elements['a'] < 0 and a_before >0):
            continue
        if(chi < minChi):
            minChi= chi
            bestDets = cands[x]
            triplet = testTrip
    if(minChi > chiThresh):
        return False, triplet
    cands.remove(bestDets)
    return True, triplet

def chooseMin2(triplet, cands, chiThresh):
    print('start choosemin2')
    minChi = 1000
    bestDets = []
    size = len(cands)
    total = size*(size-1)/2
    counter = 0
    time0 = time.time()
    for x in range(size):
        for y in range(x):
        #print(triplet)
            counter += 1
            printPercentage(counter, total, time.time()-time0)
            testTrip = testCands(triplet, (cands[x], cands[y]))
            #print('test:' + str(testTrip.realLength()))
            #print('trip:' + str(triplet.realLength()))
            if(testTrip.realLength() < triplet.realLength()+2):
                continue
            chi = testTrip.getChiSq()
            #print(chi)
            #print(testTrip)
            #print(chi)
            if(chi < minChi):
                print('changing to ' + str(chi))
                minChi = chi
                bestDets = (cands[x], cands[y])
                triplet = testTrip
            #print(minChi)
    if(minChi > chiThresh):
        print('smallest: ' + str(minChi))
        return False, triplet
    
    cands.remove(bestDets[0])
    cands.remove(bestDets[1])
    return True, triplet

def chooseMin4(cands):
    print('start choosemin')
    minChi = 10000000000
    bestDets = [cands[0], cands[1], cands[2]]
    triplet = []
    for x in range(len(cands)):
        for y in range(x):
            for z in range(y):
                for w in range(z):
               # print(triplet)
                    
                    testTrip = testCands(Triplet([]),[cands[x],cands[y],cands[z], cands[w]])
                    chi = testTrip.getChiSq()
                    #print(chi)
                    #print(testTrip)
                    if(chi < minChi):
               # print('changing to ' + str(chi))
                        minChi = chi
                        bestDets = (cands[x], cands[y], cands[z], cands[w])
                        triplet = testTrip
    cands.remove(bestDets[0])
    cands.remove(bestDets[1])
    cands.remove(bestDets[2])
    cands.remove(bestDets[3])
    print('********************chosen***********************')
    print(triplet)
    return triplet

def chooseMin5():
    return

def extractTrip(triplet, chiThresh=30):
    cands = []
    trip = []
    print('\nsize of triplet: ' + str(len(triplet.dets)))
    for det in triplet.dets:
        if(det.lookAhead == -1):
            cands.append(det)
        else:
            trip.append(det)

    if(len(trip)!=3 and len(cands)==0):
        cands = trip[:]
        trip = chooseMin4(cands) 
    else:
        trip = Triplet(trip)
    if(trip == -1):
        return -1
    canAdd = True

    #print(trip)
    while(canAdd):
        canAdd, temp = chooseMin2(trip, cands, chiThresh)
        if(canAdd):
            print('added')
            trip = temp
    if(len(trip.dets)>3):
        canAdd = True
        while(canAdd):
            canAdd, temp = chooseMin1(trip, cands, chiThresh)
            if(canAdd):
                trip = temp
    print('Final size: ' +str(len(trip.dets)))
    return trip
   
def iterateTrips(tripList, savename='default'):
    extList = []
    time0 = time.time()
    counter = 0
    chiThresh = 50
    for trip in tripList:
        printPercentage(counter, len(tripList), time.time()-time0)
        time1 = time.time()
        counter += 1
        trip.chiSq=-1
        if(trip.getChiSq() < chiThresh):
            extList.append(trip)
            continue
        newTrip = extractTrip(trip, chiThresh)
        if(newTrip.realLength()>3):
            extList.append(newTrip)
        print('runtime: ' + str(time.time()-time1))
    writeTriplets(extList, savename + '.txt', False, False)
    pickleTriplets(extList, savename+'.pickle')


def main():
    args = argparse.ArgumentParser()
    args.add_argument('tripList', help='path to pickle file of a list of triples')
    args.add_argument('-t', '--tripNum', help='only check one triplet')
    args = args.parse_args()
    
    print('loading file: ' + args.tripList)
    tripList = pickle.load(open(args.tripList, 'rb'))
    if(args.tripNum):
        tripList = [tripList[int(args.tripNum)]]
      
    splitName = args.tripList.split('+')
    chunkName = splitName[0].split('/')[-1]
    savename = splitName[-1].split('.')[0]
    savename = chunkName + '+extractedTriplets+' + savename
    iterateTrips(tripList, savename)
    
if __name__ == '__main__':
    main()
