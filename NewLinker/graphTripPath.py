import sys
import os
tnopath = os.environ['TNO_PATH']
sys.path.insert(0, tnopath)
import matplotlib
matplotlib.use('Agg')
import LinkerLib as ll
from LinkerLib import Triplet
from LinkerLib import Detection
from LinkerLib import printPercentage
import growTriplets as gt

import time
import pickle
import numpy as np
import matplotlib.pyplot as plt
import re
import argparse


def flatten(coords_list):
    xs = []
    ys = []
    for coord in coords_list:
        x = coord[0]
        if(x > 180):
            x -= 360
        xs.append(x)  
        ys.append(coord[1])
    return xs, ys

def graph_points(new, orig, path, savename, 
                a=0, e=0, i=0, chi=0, miss=[], fakePreds=[], realLength=0):
    bluex, bluey = flatten(new)
    greenx, greeny = flatten(orig)
    redx, redy = flatten(path)
    missx, missy = flatten(miss)
    fakex, fakey = flatten(fakePreds)
    plt.figure()
    if(len(miss) > 0):
        plt.scatter(missx, missy, color='purple', label='missed triplets')
    if(len(fakePreds) > 0):
        plt.plot(fakex, fakey, color='red', label='predicted fake')
    plt.scatter(bluex, bluey, color='red', label='new cands')
    plt.plot(redx, redy, color='gray', label='predicted orbit')
    plt.scatter(greenx, greeny, color='red', label='original cands')
    plt.xlabel('Right Ascension', fontsize = 22, color = 'black')
    plt.ylabel('Declination', fontsize = 22, color = 'black')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.title('a=' + str(a) + ' e=' + str(e) + ' i=' + str(i) +  
                '\nn=' + str(len(bluex)+len(greenx)) + ' chisq=' + str(chi) + 
                ' realLength=' + str(realLength)) 
   # plt.legend(loc=1)
    plt.savefig(savename+'.png')
    plt.close('all')

def findMiss(trip, fakes):
    missFakes = []
    for fake in fakes:
        if(not fake in trip.dets):
            missFakes.append((fake.ra, fake.dec))
    return missFakes

def graph_triplets(triplets, savename1, orbs, fakeDict=0):
    print('\nplotting graphs')
    counter = 0
    time0 = time.time()
    for trip in triplets:
        printPercentage(counter, len(triplets), time.time()-time0)
        savename = 'triplet' + str(counter) + '+' + savename1
        counter += 1
        orgPoints = []
        newPoints = []
        for det in trip.dets:
            if(det.lookAhead==-1):
                newPoints.append((det.ra, det.dec))
            else:
                orgPoints.append((det.ra, det.dec))
        if(not fakeDict == 0):
            fakeid = [x.fakeid for x in trip.dets if x.lookAhead != -1]
            missPoints = findMiss(trip, fakeDict[fakeid[0]])
        else:
            missPoints = []
        trip.sortByMjd()
        mjd1 = trip.dets[0].mjd
        mjd2 = trip.dets[-1].mjd
        predPoints = generate_predictions(trip, mjd1-365, mjd2+365, orbs)
        if(not fakeDict == 0):
            fakePoints = generate_fake_preds(trip, mjd1-365, mjd2+365)
        else:
            fakePoints = []
        try:
            a = trip.elements['a']
            e = trip.elements['e']
            i = trip.elements['i']
        except TypeError:
            a = -1
            e = -1
            i = -1
        graph_points(newPoints, orgPoints, predPoints, savename, 
            a, e, i, trip.getChiSq(), missPoints, fakePoints, trip.realLength())

def generate_predictions(trip, mjd1, mjd2, orbfile):
    outname = gt.writeNites([trip], 
                range(int(mjd1), int(mjd2), 2), 2, 'graph_path.fits')
    trackToDf = gt.callMjdPrediction(outname, outname.split('.')[0] + '.pred', orbfile)
    coords =[ ]
    for trackid in trackToDf.iterkeys():
        for mjd in sorted(trackToDf[trackid].iterkeys()):
            ra = trackToDf[trackid][mjd].RA
            dec = trackToDf[trackid][mjd].DEC
            coords.append((ra,dec))
    return coords

def generate_fake_preds(trip, mjd1, mjd2):
    fakeid = trip.majFake()
    fakeDets = []
    for det in trip.dets:
        if(det.lookAhead != -1):
            fakeDets.append(det)
    #fakeDets = [x for x in trip.dets if x.fakeid == fakeid]
    fakeTrip = Triplet(fakeDets)
    coords = []
    for days in range(int(mjd1), int(mjd2), 1):
        coordinate = fakeTrip.predictPos(days)
        coor = list(coordinate[0])
        if coor[0] > 180:
            coor[0] = coor[0]-360
        coords.append(coor)
        #print(str(coor) + ' ' + str(days))
    return coords

def main():
    args = argparse.ArgumentParser()
    args.add_argument('triplets', help='pickle file of list of triplets')
    args.add_argument('-c', '--csvDetectionFile', 
                help='path to list of detection')
    args = args.parse_args()
    triplets = pickle.load(open(args.triplets, 'rb'))
#    triplets = triplets[:50]
    fakeDict = []
    if(args.csvDetectionFile):
        fakeDict = ll.fakeDict(args.csvDetectionFile)
    else:
        fakeDict = 0
    orbs = args.triplets
    ext = orbs.split('+')[-1].split('.')[0]
    name = 'merged+'+ext+'.orbit'
    if('/' in orbs):
        orbs = orbs.rsplit('/',1)[0] + '/' + name 
    else:
        orbs = name
    
#    triplets = [triplets[0]]     

    graph_triplets(triplets, 
        args.triplets.split('+')[-1].split('.')[0], orbs, fakeDict)

if __name__ == '__main__':
    main()
