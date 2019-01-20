import os
import sys
tnopath = os.environ['TNO_PATH']
sys.path.insert(0, tnopath)
import numpy as np
import pandas as pd
import scipy.spatial as sp
import argparse
import pickle
import time
import LinkerLib as LL
from LinkerLib import Triplet, Detection

'''
1. Make a dictionary that maps MJD range to list of corresponding detections
2. Make a dictionary that maps MJD range to KD tree for the detections from 1
3. Make prediction for mid point in each MJD range for each triplet and search for nearby points using KD tree
4. Validate each iteration

Note: I should check how much the object moves in a given time span
'''


def mjd_det_dict(dets, interval=20):
    mjd_dict = dict()
    for det in dets:
        mjd = (int(det.mjd) / interval) * interval
        mjd_dict.setdefault(mjd, []).append(det)
    return mjd_dict


def mjd_kd_tree_dict(mjd_det):
    kd_dict = dict()
    detlist_dict = dict()
    for key, value in mjd_det.iteritems():
        kd_dict[key] = sp.cKDTree([(x.ra, x.dec) for x in value])
        detlist_dict[kd_dict[key]] = value
    return kd_dict, detlist_dict
    


def mjd_generator(mjd_det):
    arr = []
    for mjd in mjd_det:
        arr.append(mjd)
    return np.array(arr)


def angle_diff(ang1, ang2):
    larger = max(ang1, ang2)
    smaller = min(ang1, ang2)
    return larger-smaller


def search_radius(trip, mjd, ra, dec, interval=20):
    (ra_b, dec_b), erra_b, errb_b, pa_b = trip.predictPos(mjd - interval)
    (ra_a, dec_a), erra_a, errb_a, pa_a = trip.predictPos(mjd + interval)
    erra_b /= 3600
    erra_a /= 3600
    
    max_distance = np.sqrt(max((2*erra_b+angle_diff(ra, ra_b)) ** 2 + 
                                (2*erra_b+angle_diff(dec, dec_b)) ** 2,
                                (2*erra_a+angle_diff(ra, ra_a)) ** 2 + 
                                (2*erra_a+angle_diff(dec, dec_a)) ** 2))
    '''
    if max_distance > 10:
        print('ra: ' + str(ra))
        print('ra_b: ' + str(ra_b))
        print('ra_a: ' + str(ra_a))
        print('dec: ' + str(dec))
        print('dec_b: ' + str(dec_b))
        print('dec_a: ' + str(dec_a))
        print('erra_b: ' + str(erra_b))
        print('erra_a: ' + str(erra_a))
    '''
    return max_distance


def withinEllipse(erra, errb, pa, delta_ra, delta_dec, errSize=2):
    erra /= 3600
    errb /= 3600
    pa = pa + 90
    pa = np.radians(pa)
    x = np.cos(pa) * delta_ra - np.sin(pa) * delta_dec
    y = np.sin(pa) * delta_ra + np.sin(pa) * delta_dec
    return x ** 2 / (errSize * errb) + y ** 2 / (errSize * erra) <= 1


def filter_candidates(trip, candidates, errSize=2):
    filtered = []
    for candidate in candidates:
        '''
        # create dummy detection
        coord, erra, errb, pa = trip.predictPos(candidate.mjd)
        predictionDet = Detection(coord[0], coord[1], candidate.mjd)
        predictionDet.erra = erra
        predictionDet.errb = errb
        predictionDet.pa = pa
        if(withinEllipse(candidate, predictionDet, errSize) and (not candidate in trip.dets)):
            filtered.append(candidate)
        '''
        (ra, dec), erra, errb, pa = trip.predictPos(candidate.mjd)
        if withinEllipse(erra, errb, pa, angle_diff(candidate.ra, ra), 
                angle_diff(candidate.dec, dec), errSize):
            filtered.append(candidate)
    return filtered


def find_candidates(trips, dets, interval=2, chithresh=50):
    mjd_det = mjd_det_dict(dets, interval)
    mjd_kd_tree, kd_tree_detlist = mjd_kd_tree_dict(mjd_det)

    mjd_arr = mjd_generator(mjd_det)
    half = float(interval / 2)
    cand_count = 0

    grownTrips = []
    counter = 0
    time1 = time.time()
    size = len(trips)
    
    #TODO delete
    f = open('tempFile.predict', 'w+')
    f.write('RA,DEC,ERR_A,MJD,TRACK_ID\n')
    for trip in trips:
        counter += 1
        LL.printPercentage(counter, size, time.time()-time1)
        if trip.getChiSq() > chithresh:
            #print(trip)
            continue
        t0 = time.time()
        count_first_stage = 0
        count_second_stage = 0
        totalCands = []
        #print(mjd_arr)
        for mjd in mjd_arr:
            midpoint = mjd + half
            (ra, dec), erra, errb, pa = trip.predictPos(midpoint)
            
            #TODO delete
            f.write(str(ra) + ',')
            f.write(str(dec) +',')
            f.write(str(erra) + ',')
            f.write(str(mjd) + ',')
            f.write(str(trip.trackid) + '\n')

            if(erra>1000):
                #print('mjd:'+str(mjd))
                #print(trip)
                continue
            ell = ((ra, dec), erra, errb, pa)
                        #radius = min(5, search_radius(trip, midpoint, ra, dec, interval))
            radius = search_radius(trip, midpoint, ra, dec, interval)

            kdtree = mjd_kd_tree[mjd]
            dets = kd_tree_detlist[kdtree]
            candidates = [dets[i] for i in kdtree.query_ball_point((ra, dec), radius)]
            #print('ra:' + str(ra) + ' dec:' + str(dec) + ' radius:'+str(radius))
            '''
            if((mjd <= 57277 and mjd >=57275) or 
                (mjd <= 56540) and (mjd >= 56538)):
                #print('search radius' + str(radius))
                for det in dets:
                    print(det.toStr())
                print('************mjd:'+str(mjd)+'*************')
                print(trip)
                print(ell)
                print('number of cands:' + str(len(candidates)))
                for cand in candidates:
                    print(cand.toStr())
            '''

            if(len(candidates) >= 100):
                print('num cands:' + str(len(candidates)) + 
                        ' search radius:' + str(radius) + 
                        ' ellipse:' + str(ell))
                #print(trip)
            filtered = filter_candidates(trip, candidates)
            count_first_stage += len(candidates)
            count_second_stage += len(filtered)
            totalCands.extend(filtered)
            #if len(filtered) > 0:

            #    print("Found " + str(len(filtered)) + " candidate detections")
        #print(str(count_first_stage) + ' through the first stage')
        #print(str(count_second_stage) + ' through the second stage')
        #print('Time taken: ' + str(time.time() - t0) + '\n')
        cand_count += count_second_stage
        for cand in totalCands:
            #print(cand.toStr())
            cand.lookAhead = -1
            trip.addDetection(cand)
        if(trip.realLength()>3):
            grownTrips.append(trip)
            #print(trip)
        #print('******************************')
    f.close()
    return grownTrips

def adjust_trip_ra(trips):
    for trip in trips:
        for det in trip.dets:
            if det.ra > 180:
                det.ra = det.ra - 360


def grow_directory(directory, detections):
    print("Finding candidates")
    for filename in os.listdir(directory):
        if filename.endswith(".pickle"): 
            print(os.path.join(directory, filename))
            print("********************************************************")
            triplets = pickle.load(open(directory + '/' + filename, 'rb'))
            #adjust_trip_ra(triplets)
            t0 = time.time()
            grownTrips = find_candidates(triplets, detections, 2)
            t = time.time() - t0
            print('On average ' + str(t / len(triplets)) + ' seconds per triplet')
        else:
            continue


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('triplets', help='path to triplet pickle file')
    parser.add_argument('detections', help='path to detection list csv file')
    args = parser.parse_args()

    #dets = LL.wrapDets(args.detections[0])

    
    #print("Loading triplets and detections")
    triplets = []
    with open(args.triplets,'rb') as f:
        triplets = pickle.load(f)
    #triplets = [triplets[303], triplets[304]]
    detections = LL.wrapDets(args.detections)
    #adjust_trip_ra(triplets)

    splitName = args.triplets.split('+')
    chunkName = splitName[0].split('/')[-1]
    savename = splitName[-1].split('.')[0]
    savename = chunkName +'+crossCampaignTriplets+' + savename
    
    print('Finding candidates')
    t0 = time.time()
    grownTriplets = find_candidates(triplets, detections, 2)
    t = time.time()-t0
    print('On average ' + str(t / len(triplets)) + ' seconds per triplet')

    print('\nsaving predictions to '+savename)

    LL.writeTriplets(grownTriplets, savename +'.txt', False)
    LL.pickleTriplets(grownTriplets, savename+'.pickle')

    #grow_directory(args.triplets[0], dets)




if __name__ == '__main__':
    main()

