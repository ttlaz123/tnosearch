import pickle
import os
import argparse
import time
import sys
tnopath = os.environ['TNO_PATH']
sys.path.insert(0, tnopath)
from LinkerLib import Triplet
from LinkerLib import pickleTriplets
from LinkerLib import writeTriplets
from LinkerLib import printPercentage
#splits a list of triplets into 'chunks' number of chunks
#plus one extra if #chunks doesn't divide the size of the triplets
'''
input: list of triplets, number of chunks to split the triplets
output: a list of list of triplets
'''
def splitList(tripList, chunks):
    print('spliting list of size ' + str(len(tripList)) + 
        ' into ' + str(chunks) + ' sized chunks')
    size = len(tripList)
    #sorry i made chunks be the chunksize instead of 
    # number of chunks i promise i'll fix the variables
    # and documentation soon
    chunkSize = chunks
    splitTrips = []
    time0 = time.time()
    detCounter = 0
    chunkList = []
    x = 0
    while(len(tripList) > 0):
        x += 1
        printPercentage(x, size, time.time()-time0)
        trip = (tripList.pop())
        chunkList.append(trip)
        add = len(trip.dets)
        combs = add*(add-1)/2
        detCounter += combs
        if(detCounter > chunkSize):
            detCounter = 0
            splitTrips.append(chunkList)
            chunkList = []
    return splitTrips


'''
input: a folder containing chunks of triplets, number of chunks to split the triplets
output: writes out the triplets
'''

def splitFolder(foldername, prefix, chunks):
    files = os.listdir(foldername)
    tripChunk = []
    detCounter = 0
    filCounter = 0
    chunkCounter = 1
    time0 = time.time()
    savename = files[0].split('+')[-1].split('.')[0]
    saveName = prefix + '+' + savename

    errors = []
    for f in files:
        filCounter += 1
        printPercentage(filCounter, len(files), time.time()-time0)
        name = f.split('+', 1)[-1]
        if name.endswith('.pickle') and name.startswith(prefix):
            with open(foldername + f, 'rb') as pickleFile:
                try:
                    tripList = pickle.load(pickleFile)
                except:
                    er = 'error: ' + f
                    print('\n' + er)
                    errors.append(er)   
                    continue

                while(len(tripList) > 0):
                    trip = tripList.pop()
                    combs = (len(trip.cands)*(len(trip.cands)-1)/2)
                    detCounter += combs
                    if(detCounter > chunks):
                        writeTriplets(tripChunk, 'splitchunk' + 
                            '{0:06d}'.format(chunkCounter) + 
                            '+' + saveName + '.txt', False)
                        pickleTriplets(tripChunk, 'splitchunk' + 
                            '{0:06d}'.format(chunkCounter) + 
                            '+' + saveName + '.pickle')  
                        tripChunk = [trip]
                        chunkCounter += 1
                        detCounter = combs
                    else:
                        tripChunk.append(trip)

    chunkCounter = 0 
    writeTriplets(tripChunk, 'splitchunk' + 
                            '{0:06d}'.format(chunkCounter) + 
                            '+' + saveName + '.txt', False)
    pickleTriplets(tripChunk, 'splitchunk' + 
                            '{0:06d}'.format(chunkCounter) + 
                            '+' + saveName + '.pickle')  
    return errors

def main():
    args = argparse.ArgumentParser()
    args.add_argument('-t', '--triplets', help='path to pickle file; ' + 
             'file has format triplet+SNOBS_SEASON###_ML0#.pickle')
    args.add_argument('-n', '--numChunks',  
        help='number of detections in each chunk')
    args.add_argument('-f', '--folder', help='folder of triplets')
    args.add_argument('-p', '--prefix', help='prefix of files in folder')
    args = args.parse_args()
    
    #save triplets

    
    numChunks = 100000
    if(args.numChunks):
        numChunks = int(args.numChunks)
    if(args.triplets):
        saveName = args.triplets[0].split('+')[-1].split('.')[0]
        saveName2 = args.triplets[0].split('+')[-2]
        if '/' in saveName2:
            saveName2 = saveName2.split('/')[-1]
        print(saveName2)
        print('open pickle file')
        with open(args.triplets, 'rb') as pickleFile:
            print('loading pickle file')
            triplets = pickle.load(pickleFile)
        tripChunks = splitList(triplets, numChunks)
        for x in range(len(tripChunks)):
            writeTriplets(tripChunks[x], 'splitchunk' + '{0:06d}'.format(x) + 
                    '+' + saveName2 + '+' + saveName + '.txt', False)
            pickleTriplets(tripChunks[x], 'splitchunk' + '{0:06d}'.format(x) + 
                    '+' + saveName2 + '+' + saveName + '.pickle')  
    elif(args.folder):
        errors = splitFolder(args.folder, args.prefix, numChunks) 
        print(errors)

if __name__ == '__main__':
    main() 
    
