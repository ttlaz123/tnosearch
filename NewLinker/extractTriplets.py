import sys
import os
Linker_Path = os.environ["TNO_PATH"]
sys.path.insert(0, Linker_Path)

import argparse
import pickle
import time
import pandas as pd 
from astropy.table import Table

from LinkerLib import Detection
from LinkerLib import Triplet
from LinkerLib import writeTriplets
from LinkerLib import pickleTriplets
from LinkerLib import printPercentage
import LinkerLib as LL


TRACK_ID_HEADER = 0
DET_DICT = 0

'''
DEPRECATED
Input: --a list of original detections
       --a list of candidates detections
       --chiThresh per degree of freedom
       --number n for combinations of how many detections
Output: bool - whether it was possible to add two cands
        the two cands that contribute the least chiSq

Writes all combinations of n new candidate(s) to a file, then calls the 
orbit fitter in C++, which returns the outputted chisq and params in a 
separate file, which is then used to determine the best fit.

'''
def chooseMin(dets, cands, chiThresh, num):
    if(num == 2):
        testDict = choose2(dets, cands)
    elif(num == 1):
        testDict = choose1(dets, cands)
    assert(len(testDict) > 0)
    inputname = LL.writeDetToOrb(testDict, verbose=False)
    ### call C function
    ### read in file with TRACKID and chisq and stuff
    trackParams = callOrbitFitter(inputname) 
    trackParams.CHISQ.astype(float).fillna(99999999.9)
    threshPerDof = chiThresh*(len(dets) + num)
    # get the two best cands
    minRow = trackParams['CHISQ'].idxmin()
    minChi = trackParams.at[minRow, 'CHISQ']
    if(minChi > threshPerDof):
        return False, []
    
    minTrack = trackParams.at[minRow, 'ORBITID']
    bestCands = testDict['{0:05d}'.format(minTrack)]
  
    newCands = [x for x in bestCands if x not in dets]
    return True, newCands


'''
DEPRECATED
Input: --List of triplets
       --chiThresh per number of degrees of freedom

Output: --List of triplets extracted from respective candidates with realLength()
            greater than 3 and chiSq per ndof less than chiThresh
'''
def iterateTrips(tripList, chiThresh=50):
    extList = []
    time0 = time.time()
    size = 0
    for trip in tripList:
        size+= len(trip.cands)*(len(trip.cands)-1)/2
    counter = 0
    for x in range(len(tripList)):
        trip = tripList[x]
        trip.cands = removeDup(trip)
        counter += len(trip.cands)*(len(trip.cands)-1)/2
        printPercentage(counter, size, time.time()-time0)
        newTrip = extractTrip(trip, chiThresh)
        if(len(newTrip.dets)>4):
            extList.append(newTrip)
    return extList

'''
DEPRECATED
Input: --a triplet containing a list of candidates
       --chiThresh per degree of freedom

Output: --the new Nlet with added detections from the list of candidates

'''
def extractTrip(trip, chiThresh):
    #print('\nnumber of candidates: ' + str(len(trip.cands)))
    dets = trip.dets[:]
    cands = trip.cands[:]
    canAdd = True

    # add detections 2 at a time
    if(len(trip.cands) < 2):
        canAdd = False
    while(canAdd):
        canAdd, min2 = chooseMin(dets, cands, chiThresh, 2)
        dets = dets + min2
        for mi in min2:
            cands.remove(mi.objid)
        if(len(cands) < 2):
            canAdd = False
    
    # add any extra one detection if possible
    if(len(dets) > 3 and len(cands) > 0):
        canAdd = True
        while canAdd:
            canAdd, min1 = chooseMin(dets, cands, chiThresh, 1)
            dets = dets + min1
            for mi in min1:
                cands.remove(mi.objid)
            if(len(cands) < 1):
                canAdd = False

    #print('Final Size: ' + str(len(dets)))
    return Triplet(dets)




'''
Input: --name of input file
Output: --a pandas dataFrame with headers TRACK_ID,CHI_SQ

Calls the C++ orbit fitter, then extracts relevant information

'''
def callOrbitFitter(inputname):
    outputname = inputname.split('.')[0] + '.orbit'
    cmd = 'BulkFit -observationFile=' + str(inputname)
    cmd = cmd + ' -orbitFile=' + str(outputname)
    print('calling C function...')
    os.system(cmd)
    print('reading Table')
    dat = Table.read(outputname, format='fits')
    print('removing columns')
    print(dat.info())
    dat.remove_columns(['ABG', 'ABGCOV'])
    
    masked_err = dat['FLAGS']==0 
    dat = dat[masked_err]
    masked_err = dat['CHISQ']!=0.0
    dat = dat[masked_err]
    print(dat.info())
    print('converting to pandas')
    df = dat.to_pandas()
    
    return df

'''
Input: --list of detections
       --list of candidate detections
Output: --a dictionary from trackID to lists of combinations of detections
'''
def choose1(dets, cands, countStart = 0, maxCount=999999):
    testDict = {}
    counter = 0
    #get every combination of 1
    for x in range(len(cands)):
        counter += 1
        if(counter > maxCount):
            print('exceeded max count: ' + cands)
            exit()
        trackID = counter+countStart
        detList = dets + [DET_DICT[cands[x]]]

        testDict[trackID] = detList
    return testDict

'''
Input: --list of detections
       --list of candidate detections
Output: --a dictionary from trackID to lists of combinations of detections
'''
def choose2(dets, cands, countStart=0, maxCount=999999):
    testDict = {}
    counter = 0
    closeCounter = 0
    #get every combination of 2
    for x in range(len(cands)):
        for y in range(x):
            det1 = DET_DICT[cands[x]]
            det2 = DET_DICT[cands[y]]
            if(abs(det1.mjd-det2.mjd) > 0.5):
                counter += 1
                if(counter > maxCount):
                    print('exceeded max count: ' + cands)
                    exit()
                trackID = (counter+countStart)
                #print(trackID) 
                detList = dets + [det1, det2]
           
                testDict[trackID] = detList
            else:
                closeCounter += 1
    return testDict,closeCounter

'''
input: --dataframe with headers ORBITID and CHISQ
       --threshold for total chisq
output: --orbitID of minimum chisq
        --minimum chisq

        returns -1,-1 if the minchi is too large or doesn't exist
'''
def getMinChisq(df, thresh=50, header='CHISQ', header2='FLAGS'):
    df.CHISQ.astype(float).fillna(999999999.9)
    try:
        minRow = df[header].idxmin()
    except ValueError:
        #doesn't have a zero flag orbit
        return -1, -1
    minchi = df.at[minRow, header]
    if(minchi > thresh):
        return -1, -1
    else:
        orbID = df.at[minRow, 'ORBITID']
        return orbID, minchi


'''
remove the candidates in the triplet that are already in the detections
'''
def removeDup(triplet):
    objidList = [det.objid for det in triplet.dets]
    remainList = []
    for cand in triplet.cands:
        if(not cand in objidList):
            remainList.append(cand)
    return remainList

'''
simultaneously writes every combination of 2 for every triplet into a single file
'''
def simulIterTrips(tripList, chiThresh, chunkname, savename):
    chooseDict = {}
    #newList = []
    for trip in tripList:
        #TODO
        #fake = trip.dets[0].fakeid
        #if(fake == 180126319):
        #    newList.append(trip)
            
        chooseDict[trip.trackid] = 2
    #TODO remove
    #tripList = newList
    finalList = []
    maxCount = 10000
    while len(tripList) > 0:
        print('number in final list: ' + str(len(finalList)))
        trackDict = {}
        print('number of triplets left to extend: ' + str(len(tripList)))
        contList = []
        #make dictionary of all combinations with cands of all triplets
        closeCounter = 0
        addList = []
        for trip in tripList:
            orbitIDh = int(trip.trackid)*maxCount
            #print(orbitIDh)
            if(len(trip.cands) > 1 and chooseDict[trip.trackid] == 2):
                tempDict, close = choose2(trip.dets, trip.cands, orbitIDh, maxCount)
                contList.append(trip)
                if(not tempDict):
                    addList.append(trip)
                    chooseDict[trip.trackid] = 1
            elif(len(trip.cands) > 0 and chooseDict[trip.trackid] == 1):
                tempDict = choose1(trip.dets, trip.cands, orbitIDh, maxCount)
                contList.append(trip)
            elif(len(trip.cands) == 1 and len(trip.dets)>3):
                tempDict = choose1(trip.dets, trip.cands, orbitIDh, maxCount)
                contList.append(trip)
            elif(len(trip.cands) == 0 and len(trip.dets)>5):
                # no more continue
                finalList.append(trip)
                tempDict = {}
            closeCounter += close
            trackDict.update(tempDict)
        print('Close detections not considered: ' + str(closeCounter))
        tempname = 'siftRequest+' + chunkname + '+' + savename + '.fits'
        inputname = LL.writeDetToOrb(trackDict, tempname)
        if(inputname == False):
            tripList = addList
            continue
        print('running orbit fitter...')
        time0 = time.time()
        trackParams = callOrbitFitter(inputname)
        print('done after ' + str(time.time()-time0) + ' seconds')
        print('extracting minimums')
        for trip in contList:
            combInds = ((trackParams['ORBITID'] >= int(trip.trackid)*maxCount) &
                        (trackParams['ORBITID'] < ((int(trip.trackid)+1)*maxCount)))
            #print(trip.trackid)
            #print(len(trackParams))
            combs = trackParams[combInds]
            #print(len(combs))
            num = len(trip.dets) + chooseDict[trip.trackid]
            
            orbID, minchi = getMinChisq(combs, chiThresh*num)
            '''
            if(trip.sameFake() == 180079561):
                print(trip)
                print(minchi)
                print(combs)
            '''
            if(int(orbID) == -1):
                #print(trip)
                # no minchi that satisfied threshold
                if(chooseDict[trip.trackid] == 2 and trip.realLength()>3):
                    chooseDict[trip.trackid] = 1
                    addList.append(trip)
                elif(chooseDict[trip.trackid] == 1 and trip.realLength()>5):
                    finalList.append(trip)
                continue
            dets = trackDict[orbID]
            newdets =[]
            newdets = [x for x in dets if x not in trip.dets]
            for det in newdets:
                det.lookAhead = -1
                trip.dets.append(det)
            trip.cands = [x for x in trip.cands if DET_DICT[x] not in trip.dets]
            #print(trip)
            if(len(trip.cands)>1):
                trip.chiSq = minchi
                addList.append(trip)
            elif(len(trip.cands)==1):
                addList.append(trip)
                chooseDict[trip.trackid] == 1
            elif(len(trip.cands)==0 and trip.realLength()>5):
                trip.chiSq = minchi
                finalList.append(trip)
            else:
                pass
                #print('what')
        tripList = addList
    return finalList

def numcombs(tripList):
    count = 0
    trackCount = 0
    newList = []
    for trip in tripList:
        trackCount += 1
        trip.trackid = trackCount%100000
        l = len(trip.cands)
        count += l*(l-1)/2
        objids = [x.objid for x in trip.dets]
        trip.cands = [x for x in trip.cands if x not in objids]
        #if(trip.sameFake() == 180079561):
        newList.append(trip)
    print('total combs: ' + str(count))
    return newList

def main():
    args = argparse.ArgumentParser()
    args.add_argument('tripList', help='path to pickle file of a list of triplets;' + 
              'format: splitchunk######+crossCampaignTriplets+SN???_SEASON###_ML##.pickle')
    args.add_argument('csvFile', help='csv File if triplets are objids')
    args.add_argument('-t', '--chiThres', help='threshold for chisq per degree of freedom')
    args = args.parse_args()

    #getting name
    splitName = args.tripList.split('+')
    if('/' in splitName[0]):
        chunkName = splitName[0].split('/')[-1]
    else:
        chunkName = splitName[0]
     
    savename = splitName[-1].split('.')[0]
    savename = chunkName + '+extractedTriplets+' + savename
    print('saving to: ' + savename)

    print('loading file: ' + args.tripList)
    with open(args.tripList, 'rb') as f:
        tripList = pickle.load(f)
    tripList = numcombs(tripList)
    
    global DET_DICT
    DET_DICT = LL.objidDictionary(args.csvFile)

    chiThresh = 10
    if(args.chiThres):
        chiThresh = int(args.chiThres)
    
    extList = simulIterTrips(tripList, chiThresh, chunkName, savename) 
    writeTriplets(extList, savename + '.txt', True)
    pickleTriplets(extList, savename + '.pickle')

if __name__ == '__main__':
    main()
