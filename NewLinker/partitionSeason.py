'''
Takes a csv file with a table of detections and splits them up into several regions
Needs to be run before linkDetections

'''


import sys
import os
tnopath = os.environ['TNO_PATH']
sys.path.insert(0, tnopath)
import numpy as np
import pandas as pd
import argparse
import math
import GammaTPlotwStatTNOExFaster as gts
import time
import pickle
from LinkerLib import Region
from LinkerLib import Detection
from LinkerLib import printPercentage

#finds coordinates that would bound every object
def raDecBounds(ras, decs):
    ras = [x-360 if x>180 else x for x in ras] 
    minra = min(ras)-0.01
    maxra = max(ras)+0.01
    mindec = min(decs)-0.01
    maxdec = max(decs)+0.01
    return ((minra, maxra), (mindec, maxdec))

# returns maximum ra
'''
input: csvfile containing the boundaries of each season, season number
output: topleft and bottomright corner of the season
'''
def raDecRng(seasonBounds, season):
    print 'determining range'
    csvfile = pd.read_csv(seasonBounds)
    seasons = csvfile['Season']    
    for s in range(len(seasons)):
        row = csvfile.loc[[s]]
        if(int(row['Season']) == int(season)):
            return ((float(row['minRa']), float(row['maxRa'])), 
			(float(row['minDec']), float(row['maxDec'])))
    #season not found
    return ((-1,-1), (-1,-1))    

#initializes the boundaries of each region
'''
input: boundaries of the season, number of rows and columns to partition
output: 2D array with a region in each cell
'''
def initializeArray(minRa, maxRa, minDec, maxDec, numrow, numcol):
    #determine step for each col and row
    stepRA = (maxRa - minRa) / numcol
    stepDec = (maxDec - minDec) / numrow
    print('stepsize: ' + str(stepRA) + " " + str(stepDec)) 
    # initialize row by col array with Regions object
    # [0,0] contains min RA and min Dec
    # [numrow - 1, numcol - 1] contains max RA and max Dec
    numrow = int(numrow)
    numcol = int(numcol)
    arr = np.empty((int(numrow), int(numcol)), dtype = object)
    for i in range(numrow):
        for j in range(numcol):
            arr[i,j] = Region((minRa + j * stepRA, minRa + (j + 1) * stepRA),
                            (minDec + i * stepDec, minDec + (i + 1) * stepDec))
    return arr

#determines where the detection will be inserted
'''
input: fakeid of object, position of object, boundary of season
output: index to insert the object
'''
def insertDetection(fakeid, ra, dec, minRa, maxRa, minDec, maxDec, stepRA, stepDec):
    i=0
    j=0
    #no mag 20 fakes
    if(not fakeid>800000000):
        # a bunch of math that determines where the detection goes
        
        if(ra > 180):
            ra = ra - 360
        if ra != minRa and dec != minDec:
            j = math.ceil((ra - minRa) / stepRA) - 1
            i = math.ceil((dec - minDec) / stepDec) - 1
        elif ra == minRa:
            j = 0
            if dec == minDec:
                i = 0
            else:
                i = math.ceil((dec - minDec) / stepDec) - 1
        else:
            j = math.ceil((ra - minRa) / stepRA) - 1
            i = 0
        return i,j
    return -1, -1

# store each TNO in respective regions
# param rfile: csvfile to read from
# param row: number of rows of region
# param col: number of columns of region
def storeTNO(rfile, numrow, numcol):#, seasonBounds, season):
    #read in the data file
    print 'reading in file: ' + rfile
    df = pd.read_csv(rfile)
    df.columns = df.columns.str.lower()

   #sometimes the names are different
    df = df.rename(columns={'snobjid': 'objid', 'snfake_id': 'fakeid',
                            'ccdnum': 'ccd'})
    
    size = len(df["objid"])
    
    print 'number of detections: ' + str(size)
    
    #speeds up the process
    raList = df['ra'].tolist()
    decList = df['dec'].tolist()
    #mjdList = df['mjd'].tolist()
    
    # temporary change to take care of flux/mag issue
    #if 'mag' in df.columns:
    #    fluxList = df['mag'].tolist()
    #    fluxList = [10 ** ((31.4 - x) / 2.5) for x in fluxList]
    #else:
    #    fluxList = df['flux'].tolist()
    objidList = df['objid'].tolist()
    #expnumList = df['expnum'].tolist()
    #ccdList = df['ccd'].tolist()
    #bandList = df['band'].tolist()
    try:
        fakeidList = df['fakeid'].tolist()
    except KeyError:
        fakeidList = [0]*len(objidList)
    #determine the range of the season
    #raDec = raDecRng(seasonBounds, season)
    print('finding ra dec range')
    time0 = time.time()
    raDec = raDecBounds(raList, decList)
    print('ra dec range: ' + str(raDec) + ' after ' + 
            str(time.time()-time0) + ' secs')
    minRA = raDec[0][0]
    maxRA = raDec[0][1]
    minDec = raDec[1][0]
    maxDec = raDec[1][1]
    
    arr = initializeArray(minRA, maxRA, minDec, maxDec, numrow, numcol)
    stepRA = (maxRA - minRA) / numcol
    stepDec = (maxDec - minDec) / numrow
    
    startT = time.time()
    #insert each detection into a partition
    for y in range(size):
        if(y%1000==0):
            printPercentage(y,size, time.time()-startT)
        #wrap the object
        i,j = insertDetection(fakeidList[y], raList[y], decList[y], 
                        minRA, maxRA, minDec, maxDec, stepRA, stepDec)
        if(i>-1 and j>-1):
            arr[int(i)][int(j)].add(objidList[y])
    
    return arr

#writes the regions to a file
def writeRegion(filename, regions, numrow, numcol):
    print('writing file')
    with open(filename, 'w+') as txtFile:
        total = numrow*numcol
        counter = total
        detCounter = 0
        startT = time.time()
        for row in regions:
            for reg in row:
                txtFile.write(reg.toStr() + '\n')
                printPercentage(total-counter, total, time.time() - startT)
                counter-=1
                for det in reg.detections:
                    txtFile.write(str(det) + '\n')    
                    detCounter += 1
        print('\ninserted detections: ' + str(detCounter))

def main():
    args = argparse.ArgumentParser()
    args.add_argument('seasonFile', nargs=1,
                        help=('path to csv file for objects in the season;' +  
                        ' name has format ccdRM+SNOBS_SEASON###_ML0#.csv'))
    args.add_argument('-c', '--numCols', help='number of columns to partition season')
    args.add_argument('-r', '--numRows', help='number of rows to partition season')
    args = args.parse_args()

    numCols = 60
    numRows = 60
    if(args.numCols):
        numCols = int(args.numCols)
    if(args.numRows):
        numRows = int(args.numRows)
    regionList = storeTNO(args.seasonFile[0], numRows, numCols)
    print('\n')
    #maintain a consistent naming convention
    csvFile = args.seasonFile[0].split('+')[-1].split('/')[-1].split('.')[0]
    writeRegion('regions+' + csvFile + '.txt', regionList, 
                numCols, numRows)
    print('\n pickling file to regions+' + csvFile + '.pickle')
    with open('regions+' + csvFile + '.pickle', 'wb') as f:
        pickle.dump(regionList, f)
    print('\n')

if __name__ == '__main__':
    main()
