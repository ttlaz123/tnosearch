import sys
import os
tnopath = os.environ['TNO_PATH']
sys.path.insert(0, tnopath)

import LinkerLib as LL
from LinkerLib import Triplet
from LinkerLib import Detection
import numpy as np
from astropy.table import Table
import time
import pickle
import argparse


'''
Converts a list of Triplet objects into two tables: 
a table with orbital element values
a table with detection values

'''
def writeFitsTable(tripList, orbit=None):
    trackList = []
    aList = []
    eList = []
    iList = []
    aopList = []
    topList = []
    lanList = []
    chiList = []
    numList = []
    uniList = []
    
    detTrackList = []
    objidList = []
    raList = []
    decList = []
    magList = []
    mjdList = []
    expList = []
    bandList = []
    errList = []
    abgList = []
    abgCov = []
    aeiCov = []
    trackid = 0
    for trip in tripList:
        trackid += 1
        trackList.append(trackid)
        try:
            aList.append(trip.elements['a'])
            eList.append(trip.elements['e'])
            iList.append(trip.elements['i'])
            lanList.append(trip.elements['lan'])
            topList.append(trip.elements['top'])
            aopList.append(trip.elements['aop'])
        except TypeError:
            aList.append(-1)
            eList.append(-1)
            iList.append(-1)
            lanList.append(-1)
            topList.append(-1)
            aopList.append(-1)
        chiList.append(trip.chiSq)
        numList.append(len(trip.dets))
        uniList.append(trip.realLength())
        abgList.append(trip.abg)
        abgCov.append(trip.cov)
        aeiCov.append(trip.aeiCov)
        for det in trip.dets:
            detTrackList.append(trackid)
            objidList.append(det.objid)
            raList.append(det.ra)
            decList.append(det.dec)
            magList.append(det.mag)
            mjdList.append(det.mjd)
            expList.append(det.expnum)
            bandList.append(det.band)
            errList.append(det.posErr)
        # add to tables
    abgList = np.array(abgList)
    abgCov = np.array(abgCov)
    aeiCov = np.array(aeiCov)
    table1 = Table((trackList, aList, eList, iList, 
                        lanList,topList,aopList, 
                        chiList, numList, uniList, abgList, abgCov, aeiCov),
                names=('orbitid', 'a', 'e', 'i', 'lan', 'top', 'aop',
                            'chisq', 'num', 'nunique', 'abg', 'abgcov', 'aeicov'))
    table2 = Table([detTrackList, objidList, raList, decList, 
                        magList,mjdList,expList, 
                        bandList, errList],
                names=('orbitid', 'detectionid', 'ra', 'dec', 
                        'mag', 'mjd', 'exposure',
                            'band', 'poserr'),
                dtype=('int64', 'int64', 'f8', 'f8', 
                        'f8', 'f8', 'i8',
                            'a1', 'f8'))
    orbit = Table.read(orbit)
    table1.meta = orbit.meta
    return table1, table2


def main():
    args = argparse.ArgumentParser()
    args.add_argument('-t', '--triplets', help='list of tno orbits')
    args.add_argument('-r', '--orbit', help='fits file with orbital parameters')
    args.add_argument('-o', '--outname1', help='name for orbital elements outfile')
    args.add_argument('-d', '--outname2', help='name for detections info outfile')
    args = args.parse_args()

    triplets = []
    base = args.triplets.split('.')[0]
    if('/' in base):
        base = base.split('/')[-1]
    outname1 = 'orbitParams+' + base.split('+')[-1] + '.fits'
    outname2 = 'detParams+' + base.split('+')[-1] + '.fits'
    if(args.outname1):
        outname1 = args.outname1
    if(args.outname2):
        outname2 = args.outname2

    print('saving to ' + outname1 + ' and ' + outname2)
    print('loading from ' + base)
    with open(args.triplets, 'rb') as f:
        triplets = pickle.load(f)

    print(args.orbit) 
    orbs, dets = writeFitsTable(triplets, args.orbit)
    print('saving tables')
    orbs.write(outname1, format='fits', overwrite=True)
    dets.write(outname2, format='fits', overwrite=True)

    print('done')
if __name__=='__main__':
    main()
