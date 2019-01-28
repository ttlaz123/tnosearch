import os
import sys
import gc
try:
    tnopath = os.environ['TNO_PATH']
except KeyError:
    print('*******\nERROR: need to set environment variables TNO_PATH')
    print('TNO_PATH: location of linker libraries')
    exit()
sys.path.insert(0, tnopath)

import numpy as np
import pandas as pd
import scipy.spatial as sp
import argparse
import pickle
import time

from astropy.table import Table
from astropy.table import unique

import random
import heapq
from operator import itemgetter

import LinkerLib as LL
from LinkerLib import Triplet, Detection

def sizeof_fmt(num, suffix='B'):
    ''' By Fred Cirera, after https://stackoverflow.com/a/1094933/1870254'''
    for unit in ['','Ki','Mi','Gi','Ti','Pi','Ei','Zi']:
        if abs(num) < 1024.0:
            return "%3.1f%s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f%s%s" % (num, 'Yi', suffix)

'''
1. Make a dictionary that maps MJD range to list of corresponding detections
2. Make a dictionary that maps MJD range to KD tree for the detections from 1
3. Make prediction for mid point in each MJD range for each triplet and search for nearby points using KD tree
4. Validate each iteration

Note: I should check how much the object moves in a given time span
'''


def mjd_det_dict(dets, interval=20):
    mjd_dict = dict()
    for det in dets['mjd']:
        mjd = int(dets['mjd'][det] / interval) * interval
        mjd_dict.setdefault(mjd, []).append({'objid':det, 
                        'ra':dets['ra'][det], 
                        'dec':dets['dec'][det],
                        'mjd':dets['mjd'][det],
                        'err':dets['err'][det],
                        'expnum':dets['expnum'][det]})
    return mjd_dict


def mjd_kd_tree_dict(mjd_det):
    kd_dict = dict()
    detlist_dict = dict()
    for key, value in mjd_det.iteritems():
        kd_dict[key] = sp.cKDTree([(x['ra'], x['dec']) for x in value])
        detlist_dict[kd_dict[key]] = value
    return kd_dict, detlist_dict
    


def mjd_generator(mjd_det):
    arr = []
    for mjd in mjd_det:
        arr.append(mjd)
    return np.array(arr)


'''
input: --trackid to mjd to position and error dictionary
       --trackid of triplet
       --mjd of radius that needs to be found
       --interval of the mjd range
output: --maximum radius of where a detection should in a certain mjd
'''
def search_radius(trackMJDtoPos, trackid, mjd, interval=2, errSize=3):
    posC = trackMJDtoPos[trackid][mjd]
    pos2 = trackMJDtoPos[trackid][mjd+interval]
    pos3 = trackMJDtoPos[trackid][mjd-interval]
    #distance from the two 
    dist2 = np.sqrt((pos2['RA']-posC['RA'])**2 + (pos2['DEC']-posC['DEC'])**2)
    dist3 = np.sqrt((pos3['RA']-posC['RA'])**2 + (pos3['DEC']-posC['DEC'])**2) 
    dist2 += pos2['ERR']*errSize/3600 
    dist3 += pos3['ERR']*errSize/3600

    return max(dist2, dist3)


def withinEllipse(erra, errb, pa, delta_ra, delta_dec, errSize=2):
    erra /= 3600
    errb /= 3600
    pa = pa + 90
    pa = np.radians(pa)
    x = np.cos(pa) * delta_ra - np.sin(pa) * delta_dec
    y = np.sin(pa) * delta_ra + np.sin(pa) * delta_dec
    return x ** 2 / (errSize * errb) + y ** 2 / (errSize * erra) <= 1


'''
input: --a list of triplets
       --an array of mjd's to find the position of the prediciton
       --the interval of the range of the mjds
       --name of the outfile

output: --name of outfile
        --writes the file necessary for orbit position prediction
'''
def writeNites(trips, mjd_arr, interval, outfile):
   
    time0 = time.time()
    trackCol = []
    mjdCol = []
    for trip in trips:
        mjdDict = {}
        for mjd in mjd_arr:
            if(mjd not in mjdDict):
                mjdCol.append(mjd)
                trackCol.append(int(trip.trackid))
                mjdDict[mjd] = 1
            if((mjd + interval) not in mjdDict):
                mjdCol.append(mjd+interval)
                trackCol.append(int(trip.trackid))
                mjdDict[mjd+interval] = 1
            if((mjd - interval) not in mjdDict):
                mjdCol.append(mjd-interval)
                trackCol.append(int(trip.trackid))
                mjdDict[mjd-interval] = 1
    outTable = Table([trackCol, mjdCol], 
                names=('ORBITID', 'MJD'), dtype=['int64', 'f8'])
    print('writing to: ' + outfile) 
    outTable.write(outfile, format="fits", overwrite=True)
    print('total time: ' + str(time.time()-time0))
    return outfile 
    
'''
input: --input file to C code for orbit position predictor
       --orbital parameter fits file
output: --dictionary from trackid and mjd to ra, dec, and error
'''
def callMjdPrediction(inputFile, outputname, orbitFile, overwrite=True):
    if(not os.path.isfile(outputname) or overwrite):
        cmd = 'BulkPredict -observationFile=' + (inputFile)
        cmd = cmd + ' -orbitFile=' + orbitFile
        cmd = cmd + ' -predictFile=' + outputname
        #TODO call the code in a way that isn't a disgrace to progammers everywhere
        print('running BulkPredict...')    
        print(cmd)
        time0 = time.time()
        os.system(cmd)
        print('done after ' + str(time.time()-time0) + ' seconds')
    else:
        print('file already exists: ' + outputname)
    print('reading in ' + outputname)
    data = Table.read(outputname, format="fits")
    
    trackToDf = {}
    time0 = time.time()
    trackids = data['ORBITID'].tolist()
    mjds = data['MJD'].tolist()
    ras = data['RA'].tolist()
    decs = data['DEC'].tolist()
    errs = data['ERROR_A'].tolist()
    nextUp = 60
    for x in range(len(trackids)):
        if(time.time() - time0 > nextUp):
            LL.printPercentage(x, len(trackids), time.time()-time0)
            nextUp+=60
        ra = ras[x]
        dec = decs[x]
        err = errs[x]
        mjd = mjds[x]
        mjdDict = {'RA': ra, 'DEC': dec, 'ERR': err}
        if(trackids[x] in trackToDf):
            trackToDf[trackids[x]][mjd] = mjdDict
        else:
            trackToDf[trackids[x]] = {mjd: mjdDict}
    print('Dictionary time: ' + str(time.time() - time0))
##############################################3
    for name, size in sorted(((name, sys.getsizeof(value)) for name,value in locals().items()),
                         key= lambda x: -x[1])[:10]:
        print("{:>30}: {:>8}".format(name,sizeof_fmt(size)))
##############################################TODO
 
    return trackToDf

'''
input: --a list of triplets to grow
       --a dictionary from trackid and mjd to ra and dec
       --a kdtree for detections at each mjd interval
       --a dictionary from kdtree to a list of detections

output: --a dictionary from trackid to a list of detections
'''
def determineCandsInRadius(trips, trackMJDtoPos, 
            mjd_det, mjd_arr, interval=2, errSize=3):
    print('making kd_trees')
    mjd_kd_tree, kd_tree_detlist = mjd_kd_tree_dict(mjd_det)
    trackDict = {}
    counter = 0
    time0 = time.time()
    nextUp = 60
    maxCands = 30
    for trip in trips:
        counter += 1
        if(time.time()-time0 > nextUp):
            LL.printPercentage(counter, len(trips), time.time()-time0)
            nextUp += 60
        trackDict[trip.trackid] = []
        for mjd in mjd_arr:
            radius = search_radius(trackMJDtoPos, trip.trackid, 
                                        mjd, interval, errSize)
            kdtree = mjd_kd_tree[mjd]
            dets = kd_tree_detlist[kdtree]
            
            pos_err = trackMJDtoPos[trip.trackid][mjd]
            pos = [pos_err['RA'], pos_err['DEC']]
            dists, candKeys = kdtree.query(pos, k=maxCands, distance_upper_bound=radius)
            candidates = [] 
            for i in candKeys:
                try:
                    candidates.append(dets[i])
                except IndexError:
                    pass
            #print(len(candidates))
            if(len(candidates) > maxCands):
                print('Overflow: num=' + str(len(candidates)))
                print('radius: ' + str(radius))
            #elif(len(candidates) > 1):
            #    print('not overflow: num=' + str(len(candidates)))
            trackDict[trip.trackid].extend(candidates)
        #print(len(trackDict[trip.trackid]))
##############################################3
    for name, size in sorted(((name, sys.getsizeof(value)) for name,value in locals().items()),
                         key= lambda x: -x[1])[:10]:
        print("{:>30}: {:>8}".format(name,sizeof_fmt(size)))
##############################################TODO
 
    return trackDict

'''
input: --trackid to cands list dictionary
       --name of input file
output:
       --writes the file
'''
def writeEllipses(trackToCandsDict, outfile):
    count = sum(len(v) for v in trackToCandsDict.itervalues())
    print(count)
    trackList = np.empty(count)
    objidList = np.empty(count)
    expList = np.empty(count)
    raList = np.empty(count)
    decList = np.empty(count)
    errList = np.empty(count)
    time0 = time.time()
    nextUp = 60
    counter = 0
        
    for trackid in sorted(trackToCandsDict.iterkeys()):
        candsList = trackToCandsDict[trackid]

        if(time.time()-time0 > nextUp):
            LL.printPercentage(counter, count, time.time()-time0)
            nextUp += 60
        for cand in candsList:
            trackList[counter] = (int(trackid))
            objidList[counter] = (cand['objid'])
            expList[counter] = (cand['expnum'])
            raList[counter] = (cand['ra'])
            decList[counter] = (cand['dec'])
            errList[counter] = (cand['err'])*3600
            counter+=1
    print('writing to fits table')    
##############################################3
    for name, size in sorted(((name, sys.getsizeof(value)) for name,value in locals().items()),
                         key= lambda x: -x[1])[:10]:
        print("{:>30}: {:>8}".format(name,sizeof_fmt(size)))
##############################################TODO
    outTable = Table([trackList, objidList, expList, raList, decList, errList], 
                    names=('ORBITID', 'OBJ_ID', 'EXPNUM', 'RA', 'DEC', 'SIGMA'), 
                    dtype = ('int64', 'i8', 'i4', 'f8', 'f8', 'f8'))
    print('writing to: ' + str(outfile))
    print('total time: ' + str(time.time()-time0))
    outTable.write(outfile, format = "fits", overwrite=True)
    return outfile

'''
input: --input file to orbit prediction fitter code
output: --a dictionary from track_id to objid to how many
            sigmas away a detection is from the prediction ellipse

'''
def callSigmaDet(inputFile, outputname, orbitFile, overwrite=True):
    if(not os.path.isfile(outputname) or overwrite):
        #TODO call the c function in a way that's not a disgrace to programmers everywhere
        cmd = 'BulkProximity -observationFile=' + (inputFile)
        cmd = cmd + ' -orbitFile=' + orbitFile
        cmd = cmd + ' -chisqFile=' + outputname
        print('running BulkProximity...')
        time0 = time.time()
        os.system(cmd)     
        print('done after ' + str(time.time()-time0) + ' seconds')
    else:
        print('file already exists: ' + outputname)
    data = Table.read(outputname, format='fits')
    
    print('making dictionary...')
    time0 = time.time()
    trackDict = {}
    trackids = data['ORBITID']
    objids = data['OBJ_ID']
    chisq = data['CHISQ']
    nextUp = 60
    for x in range(len(trackids)):
        if(time.time()-time0 > nextUp):
            LL.printPercentage(x, len(objids), time.time()-time0)
            nextUp+=60
        if(trackids[x] in trackDict):
            trackDict[trackids[x]][objids[x]] = chisq[x]
        else:
            trackDict[trackids[x]] = {objids[x]: chisq[x]}
    print('done after ' + str(time.time()-time0) + ' seconds')
    return trackDict

'''
input: --list of triplets
       --list of detections
       --name of fits file with orbital parameters
       --interval parameter
       --errSize
       --name of chunk
       --name of season
output: --the same triplets but with a list of candidates added to their cands field
            these candidates are every detection in the list of detections that 
            falls into their prediciton ellipses 

'''
def find_candidates(trips, dets, orbitFile, interval=2, errSize=2, 
                    chunkname="", savename="", overwrite = True):
    maxCands = 100
    print('creating dictionaries...')
    # dict from mjd range to detections
    mjd_det = mjd_det_dict(dets, interval)
    # a list of mjdi
    mjd_arr = mjd_generator(mjd_det)
    
    mjdPred = 'mjdPredRequest+' + chunkname + '+' + savename + '.fits'
    predictfile = mjdPred.split('+',1)[-1].split('.')[0] + '.predict'
    ellRequest = 'ellSigmaRequest+' + chunkname + '+' + savename + '.fits'
    proxFile = ellRequest.split('+',1)[-1].split('.')[0] + '.prox'
    #prepare to call C function to predict mjds
    if(not os.path.isfile(mjdPred) or overwrite):
        print('\nwriting to file for prediction C function...')
        writeNites(trips, mjd_arr, interval, mjdPred)
#############################TODO
    gc.collect()    
    if(not os.path.isfile(ellRequest) or overwrite):
        # call C function, get dictionary from trackid to mjd to positions and errors
        trackMjdPos = callMjdPrediction(mjdPred, predictfile, orbitFile, overwrite)
        print('\ndetermining candidates in the maximum radius...')
        #determine the good candidates, get dictionary from trackid to list of cands
        trackToCandsDict = determineCandsInRadius(trips, 
                trackMjdPos, mjd_det, mjd_arr, interval, errSize)
        # prepare to call C function to predict ellipses
        '''
        for cands in trackToCandsDict:
            print('trackid: ' + str(cands))
            print([x.objid for x in trackToCandsDict[cands]])
        '''
        print('\nwriting to file for ellipse C function...')    
        gc.collect()
        inputFile = writeEllipses(trackToCandsDict, ellRequest)
    gc.collect()
    #call C function, get dictionary from trackid to detections to their sigmas
    trackCandsSigma = callSigmaDet(ellRequest, proxFile, orbitFile, overwrite)
    grownTrips = []

    # get all cands that are within errSize sigma of their respective predicitons
    # if greater than maxCands, then just get the smallest sigmas
    print('\nlooping through triplets...')
    time0 = time.time()
    for trip in trips:
        #print('*********')
        try:
            candsToSigma = trackCandsSigma[trip.trackid]
            #print(candsToSigma)
        except KeyError:
            candsToSigma = {}
            #print('error: ' + str(trip))
        withinThresh = dict((key, value) for key, value 
                in candsToSigma.iteritems() if value < errSize*10)
        if(len(withinThresh) > maxCands):
            withinThresh = heapq.nsmallest(maxCands, withinThresh.items()) 
        else:
            withinThresh = withinThresh.items()
        #print(withinThresh)
        trip.cands = [i[0] for i in withinThresh]
        objids = [x.objid for x in trip.dets]
        trip.cands = [x for x in trip.cands if x not in objids]
        #print(trip)
        if(len(trip.cands) >1):
            grownTrips.append(trip)
    print('done after ' + str(time.time()-time0) + ' seconds')
    return grownTrips

def efficientWrap(csvFile):
    print('reading csvFile')
    df = pd.read_csv(csvFile)
    print('converting')
    df.rename(str.lower, axis='columns', inplace=True)
    df.rename(columns={'snobjid':'objid', 'snfake_id':'fakeid',
                    'ccdnum':'ccd', 'errawin_world': 'err'}, inplace=True)
    df = df[['mjd', 'err', 'ra', 'dec', 'objid', 'expnum']]   
    detDict = df.set_index('objid').to_dict()
    return detDict

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('triplets', help='path to triplet pickle file')
    parser.add_argument('detections', help='path to detection list csv file')
    parser.add_argument('orbitFile', help='path to fits file with orbital parameters')
    parser.add_argument('-w', '--overwrite', action='store_true', 
                    help='whether to overwrite existing C files')
    args = parser.parse_args()
    
    splitName = args.triplets.split('+')
    if('/' in splitName[0]):
        chunkName = splitName[0].split('/')[-1]
    else:
        chunkName = splitName[0]
    
    saveName = splitName[-1].split('.')[0]
    savename = chunkName +'+crossCampaignTriplets+' + saveName
    
    print("Loading triplets and detections")
    triplets = []
    time0 = time.time()
    with open(args.triplets,'rb') as f:
        triplets = pickle.load(f)
    #TODO remove this
    '''
    #######
    newList = []
    for trip in triplets:
        if(trip.sameFake() == 185112789 or trip.sameFake() == 180107837):
            newList.append(trip)
    triplets = newList
    #######
    '''

   # detections = LL.wrapDets(args.detections)
    detections = efficientWrap(args.detections)
    print(detections.keys())
    print('loading and wrapping done after ' + str(time.time()-time0) + ' seconds')
    print('Finding candidates')
    t0 = time.time()
    errSize = 2
    interval = 2
    grownTriplets = find_candidates(triplets, detections, args.orbitFile,
                                interval, errSize, chunkName, saveName, args.overwrite)
    t = time.time()-t0
    print('Completed after ' + str(t) + ' seconds for ' + str(len(triplets)) + ' triplets')

    print('\nsaving predictions to '+savename)

    LL.writeTriplets(grownTriplets, savename +'.txt', False)
    LL.pickleTriplets(grownTriplets, savename+'.pickle')


if __name__ == '__main__':
    main()

