import sys
import os
tnopath = os.environ['TNO_PATH']
sys.path.insert(0, tnopath)
import math
import numpy as np
import pandas as pd
import argparse
from LinkerLib import Region
from LinkerLib import Detection
from LinkerLib import printPercentage
import LinkerLib as LL
import pickle
import time
    
#checks if two rectangles overlap given min/max ra's and dec's
'''
input: top left corner and bottom right corner of two different rectangles
output: boolean true if overlap, false otherwise
'''
def doOverlap(minRa1, maxRa1, minDec1, maxDec1,
                minRa2, maxRa2, minDec2, maxDec2):
    if(minRa1 > maxRa2 or minRa2 > maxRa1):
        return False
    if(minDec1 > maxDec2 or minDec2 > maxDec1):
        return False
    return True

#link up every detection in the cone
'''
input: a detection and a list of regions in a season
output: none, but every detection satisfying certain requirements
        is appended to the list of linked objects
'''
def linkDetection(detection, regions, detDict):
    minRa, maxRa, minDec, maxDec = detection.bounds()
    expnum = detection.expnum
    mjd = detection.mjd
    linkList = []
    #checkTime = 0
    #addTime = 0
    #coneTime = 0
    #overlapTime = 0
    for temp in regions:
        for reg in temp:
        # still checks if they overlap
            #time0 = time.time()
            '''
            overlap = doOverlap(minRa, maxRa, minDec, maxDec,
                        reg.raLo, reg.raHi, reg.decLo, reg.decHi)
            '''
            #overlap = True
            #time1 = time.time()
            if (True):
                for d in reg.detections:
                    d = detDict.get(d)
                    #time12 = time.time()
                    withinCone = detection.withinCone(d)
                    #time15 = time.time()
                    toAdd = d.expnum > expnum and d.mjd-mjd > 0.1 and withinCone
                    toAdd = toAdd and detection.magFit(d)
                    #time2 = time.time()
                    if (toAdd):
                        linkList.append(d.objid)
                    #time3 = time.time()
                    #coneTime += time15-time12
                    #checkTime += time2-time15
                    #addTime += time3-time2
            #overlapTime += time1-time0
    #print('overlapTIme ' + str(overlapTime))
    #print('coneTime: ' + str(coneTime))
    #print('checkTime: ' + str(checkTime))
    #print('addTime: ' + str(addTime))
    return linkList

## determins the subregions that overlap with the bounds of the detection
def getSubregion(det, regions, boundary, step):
    #subregion = []
    numcol = len(regions)
    numrow = len(regions[0])
    minRa, maxRa, minDec, maxDec = det.bounds()
    
    #time1 = time.time()
    ((j1,j2),(i1,i2)) = getRegionIndex(minRa, minDec, maxRa, maxDec, boundary, step, regions)
    #time2 = time.time()
    #if((j2-j1)*(i2-i1)>100):
    #    print('j1:' + str(j1) + ', j2:' + str(j2) + ', i1:' + str(i1) + ', i2:' + str(i2))
    subregion = regions[i1:i2, j1:j2]
    #print(subregion)
    #subregion = [item for sublist in subregion for item in sublist]
    '''
    for i in range(int(i1), int(i2)):
        for j in range(int(j1), int(j2)):
                subregion.append(regions[i][j])
    '''
    #print(subregion)
    #time3 = time.time()
    #print('getIndexTime: ' + str(time2-time1))
    #print('appendTime: ' + str(time3-time2))
    return subregion

"""
Function returns indices of regions (2d array) that correspond to corners
of detection bounds
"""
def getRegionIndex(ra1, dec1, ra2, dec2, boundary, step, regions):
    stepRa = step[0]
    stepDec = step[1]
    maxRa = boundary['maxRa']
    minRa = boundary['minRa']
    maxDec = boundary['maxDec']
    minDec = boundary['minDec']

    """
    ra1 and dec1 is smaller than ra2 and dec2 respectively.
    Items from boundary are from the corners of entire region(entire 2d array).
    ra1, dec1, ra2, dec2 points are the two corners for detection bounds

    If statement for corners outside the season bounds. In this case, 
    all regions are checked. Else statement for corners within season bounds.
    Returns index of regions encompassing the corners of detection bounds
    """

    if ra1 < minRa:
        ra1 = minRa
    if ra2 >= maxRa:
        ra2 = maxRa
    if dec1 < minDec: 
        dec1 = minDec
    if dec2 >= maxDec:
        dec2 = maxDec
        
    
    delRa1 = ra1 - minRa #Find distance between ra1 and dec1
    ra1Steps = math.floor(delRa1 / stepRa) #Floor to find index - index must be integer
    delRa2 = ra2 - minRa
    ra2Steps = math.floor(delRa2 / stepRa)
        
    delDec1 = dec1 - minDec
    dec1Steps = math.floor(delDec1 / stepDec)
    delDec2 = dec2 - minDec
    dec2Steps = math.floor(delDec2 / stepDec)
        #print(ra1Steps, ra2Steps, dec1Steps, dec2Steps)
    return ((int(ra1Steps), int(ra2Steps + 1)), (int(dec1Steps), int(dec2Steps + 1)))


#links detections to every other detection
'''
input: list of regions with detections in each region
output: a list of detections with their linked detections (objid)
'''        
def linkDetections(regions, detDict):
    '''
    regions[0][-1] has largest RA and lowest Dec
    regions[-1][0] has lowest RA and highest Dec
    To summarize:
    for regions[i][j], increasing i would increase Dec
    and increasing j would increase RA
    '''
    startT = time.time()
    detectionLinks = []

    print('counting')
    counter = 0
    for r in regions.flatten():
        counter += len(r.detections)
    print('done counting:' + str(counter))
    totalsize = counter

    Lo = regions[0][0]
    Hi = regions[-1][-1]

    RaLo = Lo.raLo
    RaHi = Hi.raHi
    DecLo = Lo.decLo
    DecHi = Hi.decHi

    boundary = {'minRa':RaLo, 'maxRa':RaHi, 'minDec':DecLo, 'maxDec':DecHi}

    deltaRa =  RaHi - RaLo
    deltaDec = DecHi - DecLo


    # row corresponds to dec
    # column corresponds to ra
    numrow = len(regions)
    numcol = len(regions[0])

    step = (deltaRa/numcol, deltaDec/numrow)

    getTime = 0
    getRegTime = 0
    linkTime = 0
    for regionlst in regions:
        for region in regionlst:
            startT2 = time.time()
            for detID in region.detections:
                det = detDict.get(detID)
                subregion = getSubregion(det, regions, boundary, step) 
                time3 = time.time()
                numRegions = len(subregion)*len(subregion[0])
                linkedList = linkDetection(det, subregion, detDict)
                time4 = time.time()
                detectionLinks.append((det.objid, linkedList))

                counter -= 1
                linkTime += time4-time3
                if(counter%10000 == 0):
                    printPercentage(totalsize - counter, totalsize, time.time() - startT)
                    
    return detectionLinks

#writes the links to a file
'''
input: a detection, the file to be written
output: none, file is written
'''
def printDet(detection, txtfile):
    linkedList = detection[1]
    txtfile.write('\n' + str(detection[0]) + '\n')
    txtfile.write('***********links:************** \n')
    for d in linkedList:
        txtfile.write(str(d)+'\n')

def main():
    print('Running')
    args = argparse.ArgumentParser()
    args.add_argument('regionsFile', 
                        help='path to pickle file for regions in the season; ' + 
                            'filename has format regions+SNOBS_SEASON###_ML0#.pickle')
    args.add_argument('detectionList', help='path to csvFile that contains list of detections')
    args.add_argument('-l', '--lookAhead', help='number of nights to search')
    args.add_argument('-p', '--printP', help='whether to print out percentages')
    args = args.parse_args()
    print('\nloading region file')
    time1 = time.time()
    regions = pickle.load(open(args.regionsFile))  
    print('done loading after ' + str(time.time()-time1) + ' sec')
    
    lookAhead = 50
    if(args.lookAhead):
        lookAhead=int(args.lookAhead)
    # load the regions that split the season
    print('look ahead nites: ' + str(lookAhead))
    print('making object dictionary')
    pr = True
    if(args.printP):
        pr = True
    
    detDict = LL.objidDictionary(args.detectionList, lookAhead, pr)

    orgFile = args.regionsFile.split('+')[-1].split('.')[0]
    # link up each detection with potential pairs
    print('\nlinking detections')
    startTime = time.time()
    detectionLinks = linkDetections(np.array(regions), detDict)
    print('\n')
    print('Time taken: ' + str(time.time() - startTime))
    txtfile = 'detectionLinks+' + orgFile + '.txt'
    
    #write to a text file
    with open(txtfile, 'w+') as f:
        counter = len(detectionLinks)
        startT = time.time()
        for d in detectionLinks:
            printPercentage(len(detectionLinks) - counter, len(detectionLinks),
                time.time() - startT)
            printDet(d,f)
            counter -= 1
    #save list as a picle file
    saveName = 'detectionLinks+' + orgFile + '.pickle'
    with open(saveName, 'wb') as f:
        pickle.dump(detectionLinks, f)   

    print '\n'
if __name__ == '__main__':
    main()


