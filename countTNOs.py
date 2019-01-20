import sys
import os
tnopath = os.environ['TNO_PATH']
sys.path.insert(0, tnopath)
from LinkerLib import Triplet
from LinkerLib import writeTriplets
from LinkerLib import pickleTriplets
from LinkerLib import printPercentage
import LinkerLib as LL

import numpy as np
import time
import pickle
import argparse
from astropy.table import Table

def getFakes(triplets):
    realTrips = []
    counter = 0
    time0=time.time()
    for trip in triplets:
        counter += 1
        printPercentage(counter, len(triplets), time.time()-time0)
        if (not trip.isReal()) and trip.magFilter(0.7) and trip.realLength() > 3:
            realTrips.append(trip)
    return realTrips



def getReals(triplets, nites=5):
    realTrips = []
    counter = 0
    time0=time.time()
    for trip in triplets:
        counter += 1
        printPercentage(counter, len(triplets), time.time()-time0)
        if trip.isReal() and trip.magFilter(10) and trip.realLength() > nites:
            before = len(trip.dets)
            trip.removeDupExp()
            after = len(trip.dets)
            if(before > after):
                print('changed from ' + str(before) + ' to ' + str(after))
            realTrips.append(trip)
    return realTrips

def countTriplets(folder):
    numTrips = 0
    print('counting from ' + folder)
    counter = 0
    totalList = []
    for filename in os.listdir(folder):
        if filename.endswith(".pickle"): # and filename[12:16] == 'miss':
            counter += 1
            tripList = pickle.load(open(folder+filename))
            numTrips += len(tripList)
            printPercentage(numTrips, 100, counter)
            totalList.extend(tripList)
    #LL.writeTriplets(totalList, 'misstriplets+SNALL_SEASON240_ML06.txt')
    return numTrips

def countChunks(siftedFolder, growFolder):
   return

def moveChunks(folder):
    files = os.listdir(folder)
    mvDir = 'SNALL_SEASON240_ML06_siftedChunks/'
    size = len(files)
    counter = 0
    time0 = time.time()
    for f in files:
        counter += 1
        LL.printPercentage(counter, size, time.time()-time0)
        chunkNum = f.split('+')[0][11:17]
        fileType = f.split('+', 1)[-1]
        mvFile = 'chunk' + chunkNum + '+' + fileType
        os.system('mv ' + folder+f +' ' + mvDir+mvFile)

def rmChunks(folder, start, end, season=240):
    for x in range(start, end):
        fileG = 'chunk{0:06}+goodtriplets+SNALL_SEASON'.format(x) + str(season)
        fileG = fileG + '_ML06.pickle'
        pathG = TNO_PATH+folder
        cmd = 'rm ' + pathG + fileG
        print(cmd)
        os.system(cmd)




def rmBadGrows(folder):
    numRm = 0
    print('removing from ' + folder)

    for filename in os.listdir(folder):
        if filename.endswith(".pickle"):
            #print(filename)
            tripList = []
            try:
                tripList = pickle.load(open(folder+filename))
            except EOFError:
                os.system('rm ' + folder+filename)
                txtFile = filename.split('.')[0] + '.txt'
                os.system('rm ' + folder+txtFile)
                numRm += 1
                print('removed')

            rm = False
            if(len(tripList) == 0):
                rm = False
            for trip in tripList:
                #if(trip.getChiSq() > 10):
                 #   print(trip)
                if(len(trip.dets) <= 3):
                    print('length:' + str(len(trip.dets)))
                    rm = True
                    break
                pinned = 0
                for det in trip.dets:
                    if(det.lookAhead == -1):
                        pinned += 1
                if(pinned == 0):
                    print('numPinned:' + str(pinned))
                    rm = True
                    break
             #rm = False
            if(rm):
                os.system('rm ' + folder+filename)
                txtFile = filename.split('.')[0] + '.txt'
                os.system('rm ' + folder+txtFile)
                numRm += 1
                print('removed')
    return numRm

def countFakes(tripList):
    ids = []
    goldids = []
    for trip in tripList:
        if(trip.sameFake(0)):
            goldids.append(trip.dets[0].fakeid)
        for det in trip.dets:
            if(det.fakeid < 800000000):
                ids.append(det.fakeid)
    return ids, goldids
            
def countFakesCsv(csvFile):
    detlist = LL.wrapDets(csvFile)
    tripDict = {}
    size = len(detlist)
    counter = 0
    time0 = time.time()
    for det in detlist:
        i = det.fakeid
        if(i > 800000000):
            i = 0
        if(i == 0):
            continue
        if(i not in tripDict):
            tripDict[i] = [det]
        else:
            tripDict[i].append(det)
        counter += 1
        printPercentage(counter, size, time.time()-time0)
    return list(tripDict.values())

#returns list of discoverables that aren't discovered
def missingFakes(fakeTruth, fakeTrips, thresh=2):
   
    niteThresh = 5
    truthDict = {}
    trueCounter = 0
    print('setting up dictionary')
    for trip in fakeTruth:
        trip = Triplet(trip)
        if(trip.containsNlet(100,3) and trip.realLength() >3):
            trueCounter += 1
        truthDict[trip.dets[0].fakeid] = trip
    print('dictionary done with ' + str(len(truthDict)) + ' fakes')
    print('number of discoverables: ' + str(trueCounter))
    for trip in fakeTrips:
        fakeid = trip.sameFake(thresh)
        try:
            del truthDict[fakeid]
        except KeyError:
            pass
            #print('not in dict: ' + str(trip.dets[0].fakeid))
    print('number of missing fakes: ' + str(len(truthDict)))
    discoverable = []
    for key in truthDict:
        if(key == 0):
            continue
        trip = truthDict[key]
        if(trip.containsNlet(100, 3) and trip.realLength() > niteThresh):
            discoverable.append(trip)
    print('number of missing discoverable fakes: ' + str(len(discoverable)))
    
    return discoverable

def getRealLength(filename, nites=6):
    name = filename.split('.')[0]
    if('/' in name):
        name = name.split('/')[-1]
    print('loading file:'  + filename)
    tripList = pickle.load(open(filename,'rb'))
    print('done loading')
    finalList = [trip for trip in tripList if trip.realLength() >= nites]
    writeTriplets(finalList, name + '.txt', True)
    pickleTriplets(finalList, name + '.pickle')
    

def getCands(filename, csvFile):
    tripList = pickle.load(open(filename,'rb'))
    detList = LL.wrapDets(csvFile)
    tripDict = {}
    for det in detList:
        i = det.fakeid
        if(i > 800000000):
            i = 0
        if(i not in tripDict):
            tripDict[i] = [det]
        else:
            tripDict[i].append(det)
    finalList = []
    for trip in tripList:
        fakeid = trip.sameFake(0)
        if not fakeid == 0:
            trip.cands = tripDict[fakeid]
            finalList.append(trip)
    name = filename.split('/')[-1].split('.')[0]

    writeTriplets(finalList, name + '.txt', True)
    pickleTriplets(finalList, name + '.pickle')
    return

def wrapIds(fakes, csvFile):
    detDict = LL.objidDictionary(csvFile)
    newTrips = []
    for tup in fakes:
        detList = [detDict[det] for det in tup[1]]
        newTrips.append(Triplet(detList))
    return newTrips

def insertFakes(tripList, fakename):
    fakeTable = Table.read(fakename, format='fits')
    objidList = fakeTable['CATALOG_ID'].tolist()
    fakeList = fakeTable['ORBITID'].tolist()
    fakeDict={}
    for x in range(len(fakeList)):
        fakeDict[objidList[x]] = fakeList[x]
        if(fakeList[x]==-1):
            fakeDict[objidList[x]] = 0
    print('inserting fakeids')
    for trip in tripList:
        for det in trip.dets:
            det.fakeid = fakeDict[det.objid]
    return tripList
def getMags(pickleFile, csvFile):
    objidDict = LL.objidDictionary(csvFile)
    triplets = pickle.load(open(pickleFile, 'rb'))
    for trip in triplets:
        for det in trip.dets:
            d = objidDict[det.objid]
            det.mag = d.mag
    return triplets


def main():
    args = argparse.ArgumentParser()
    args.add_argument('-d', '--folder', help='path to folder of chunks')
    args.add_argument('-f', '--fileType', help='prefex in front of the plus sign')
    args.add_argument('-c', '--csvFile', help='csv detection file')
    args.add_argument('-t', '--thresh', help='thresh for finding tno')
    args.add_argument('-p', '--pickleFile', action='append',
                        nargs='+', help='pickle file of triplets')
    args.add_argument('-r', '--reals', action='store_true', help='reals or fakes')
    args.add_argument('-i', '--fakeids', help='fakeids to store')
    args.add_argument('-l', '--realLength',action='store_true',  help='getreallength')
    args.add_argument('-m', '--mag', action='store_true', help='get real mags')
    args = args.parse_args()
    #getCands(args.pickleFile, args.csvFile)
    if(args.realLength and args.pickleFile):
        getRealLength(args.pickleFile[0][0])
        return
    if(args.fakeids and args.pickleFile and not args.csvFile):
        pickleFile = args.pickleFile[0]
        savename = pickleFile[0].split('/')[-1].split('.')[0]
        with open(pickleFile[0], 'rb') as p:
            tripList = pickle.load(p)
        newTripList = insertFakes(tripList, args.fakeids)
        writeTriplets(newTripList, savename +'.txt', True)
        pickleTriplets(newTripList, savename +'.pickle')
        return
    #return
    if(args.csvFile and args.pickleFile):
        #print(args.pickleFile)

        pickleFile = args.pickleFile[0]
        saveName = ('missingFakes+' + 
            pickleFile[0].split('+',1)[-1].split('.')[0])
        if(args.mag):
            savename = pickleFile[0].split('.')[0]
            if '/' in savename:
                savename = savename.split('/')[-1]
            trips = getMags(pickleFile[0], args.csvFile)
            LL.writeTriplets(trips, savename +'.txt')
            LL.pickleTriplets(trips, savename + '.pickle')
            return

        print('counting fakes')
        tripList = countFakesCsv(args.csvFile)
        print('\nopening pickle file')
        fakesList = []
        counter = 0
        time0 = time.time()
        for f in pickleFile:
            counter += 1
            LL.printPercentage(counter, len(pickleFile), time.time()-time0)
            if('*' in f):
                folder = f.split('/')[0]
                files = os.listdir(folder)
                fileparts = f.split('*')
                for x in range(len(folder)):
                    name = fileparts[0] + '{0:06d}'.format(x) + fileparts[1]
                    with open(name, 'rb') as p:
                        fakesList.extend(pickle.load(p))
            else:
                with open(f, 'rb') as p:
                    fakesList.extend(pickle.load(p))
        if(args.fakeids):
            fakesList =insertFakes(fakesList, args.fakeids)
        print('matching fakes')
        if(isinstance(fakesList[0], tuple)):
            fakesList = wrapIds(fakesList, args.csvFile)
        missingTrips = missingFakes(tripList, fakesList)
        writeTriplets(missingTrips, saveName+'.txt', False)
        pickleTriplets(missingTrips, saveName+'.pickle')
        return

    if(args.csvFile):
        tripList = countFakesCsv(args.csvFile)
        discoverable = 0
        for t in tripList:
            trip = Triplet(t)
            thresh = int(args.thresh)
            if(trip.realLength() > thresh):
                discoverable += 1
        print('\nnumber of fakes above ' + str(thresh) + ': ' + str(discoverable))
        return

    if(args.pickleFile):
        print(args.pickleFile)
        pickleFile = args.pickleFile[0][0]
        saveName = ('realsOnly+' + 
            pickleFile.split('+',1)[-1].split('.')[0])
        print('loading pickle file:' + pickleFile)
        triplets = pickle.load(open(pickleFile, 'rb'))
        print('done loading')
        if(args.reals):
            reals = getReals(triplets)
        else:
            reals = getFakes(triplets)
            saveName = ('fakesOnly+' + 
                pickleFile.split('+',1)[-1].split('.')[0])
        
        writeTriplets(reals, saveName + '.txt', True)
        pickleTriplets(reals, saveName + '.pickle')
        return
    files = os.listdir(args.folder)
    tripList = []
    saveName = ''
    print('opening files in: ' + args.folder)
    time0 = time.time()
    counter = 0
    uniqueCount = []
    goldCount = []
    for f in files:
        counter+=1
        printPercentage(counter, len(files), time.time()-time0)
        if(f.split('.')[-1] == 'pickle'):
            print('\nopening: ' + args.folder + f)
            with open(args.folder+f, 'rb') as tripFile:
                trips = pickle.load(tripFile)
                unique, gold = countFakes(trips)
                uniqueCount.extend(unique)
                goldCount.extend(gold)
    
    print('unique size: ' + str(np.unique(uniqueCount).size)) 
    print('gold size: ' + str(np.unique(goldCount).size))

if __name__=='__main__':
    main()
