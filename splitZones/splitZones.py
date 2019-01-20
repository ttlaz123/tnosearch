import os
import time
import glob

from collections import Counter
import numpy as np
from astropy.table import Table
from astropy.io import fits
from astropy.io.fits.hdu.hdulist import HDUList
import argparse

def getPatchExp(patch, ra0, dec0, tdb0, gamma0, dgamma, radius=3.5):
    
    print('getting patch ' + str(patch))
    params = (' ' + str(ra0) + ' ' + str(dec0) + ' ' + str(radius) + ' ' + 
                str(tdb0) + ' ' + str(gamma0) + ' ' + str(dgamma)) 
    print('GetExposurePool' + params)
    toFile = 'patch' + str(patch) + '_yr' + str(tdb0) + '.txt'
    os.system('GetExposurePool' + params + ' > ' + toFile)

def iterPatches(table, p, radius=3):
    if(not p is None):
        patches = [p]
    else:
        patches = np.unique(table['PATCH'].tolist())
    queue = {}
    for patch in patches:
        mask = table[table['PATCH']==patch]
        ra = mask['RA0'].tolist()[0]
        dec = mask['DEC0'].tolist()[0]
        tdbs = np.unique(mask['TDB0'].tolist())
        gamma = mask['GAMMA'].tolist()[0]
        dgamma = mask['DGAMMA'].tolist()[0]
        queue[patch] = [(ra, dec, tdb, 1./30, 0.05/30) for tdb in tdbs]
    time0 = time.time()
    for patch in patches:
        queues = queue[patch]
        for q in queues:
            print('time passed: ' + str(time.time()-time0))
            getPatchExp(patch, q[0], q[1], q[2], q[3], q[4], radius)

def getExposures(patch):
    patchFiles = glob.glob('patch' + str(patch) + '_yr*.txt')
    if(len(patchFiles) < 1):
        raise IOError
    exposures = set()
    for patchfile in patchFiles:
        with open(patchfile, 'r') as f:
            for line in f.readlines():
                if line[0] == '#':
                    continue
                exp = line.split()[0]
                use = line.split()[2]
                try:
                    if(int(use) == 1):
                        exposures.add(int(exp))
                except ValueError:
                    print(line.split())
    
    return exposures

def inRange(coord1, coord2, radius):
    coord1 = list(coord1)
    coord2 = list(coord2)
    if(coord1[0] > 180):
        coord1[0] -= 360
    if(coord2[0] > 180):
        coord2[0] -= 360
    '''
    term1 = np.sin(coord1[1]*np.pi/180)*np.sin(coord2[1]*np.pi/180)
    term2 = np.cos(coord1[1]*np.pi/180)*np.cos(coord2[1]*np.pi/180)
    term3 = np.cos((coord1[0]-coord1[1])*np.pi/180)

    cosa = term1+term2*term3
    dist = np.arccos(cosa)*180/np.pi
    #print('cosa: ' + str(cosa))
    #print('thres: ' + str(thresh))
    '''
    dist1 = (coord1[1]-coord2[1])
    dist2 = ((coord1[0]-coord2[0])*np.cos(dist1*np.pi/180))
    dist = dist1**2 + dist2**2

    return dist < radius**2

def expInRadius(exptable, coords, radius):
    expDict = {}
    print('finding exposures in radius:' + str(radius))
    explist = exptable['EXPNUM'].tolist()
    ralist = exptable['RADEG'].tolist()
    declist = exptable['DECDEG'].tolist()
    totalcount = 0
    track = []
    for coord in coords:
        patchnum = coord[0]
        coord = coord[1]
        exps = [explist[x] for x in range(len(explist)) 
                if inRange((ralist[x], declist[x]), coord, radius)]
        track.extend(exps)
        expDict[patchnum] = exps 
        print('patch' + str(patchnum) + ' exps: ' + str(len(exps)))
        totalcount += len(exps)
    print('unique exps: ' + str(len(np.unique(track))))
    print('total exps: ' + str(totalcount))
    return expDict

def getDetections(table, patches, patchExps=None):
    ralist = table['RA'].tolist()
    declist = table['DEC'].tolist()
    explist = table['EXPNUM'].tolist()
    ccdlist = table['CCDNUM'].tolist()
    bandlist = table['BAND'].tolist()
    mjdlist = table['MJD_OBS'].tolist()
    errlist = table['ERRAWIN_WORLD'].tolist()
    fluxlist = table['FLUX_AUTO'].tolist()
    ferrlist = table['FLUXERR_AUTO'].tolist()
    magzlist = table['MAG_ZERO'].tolist()
    maglist = table['MAG'].tolist()
    catlist = table['CATALOG_ID'].tolist()
    expDict ={}
    for x in range(len(ralist)):
        if(int(errlist[x]) > 1./3600):
            continue
        if(explist[x] in expDict):
            expDict[int(explist[x])].append(x)
        else:
            expDict[int(explist[x])] = [x]
    uniqueDict = {}
    count = 0
    for patch in patches:
        print('getting detections in patch ' + str(patch))
        outname = 'SNFIT_SEASON{0:03d}'.format(patch) + '_ML06.fits'
        if(not patchExps is None):
            exposures = patchExps[patch]
        else:
            try:
                outname = 'grow+' + outname
                exposures = getExposures(patch)
            except IOError:
                print('patch not found: ' + str(patch))
                return
        ral = []    
        decl =[]
        expl = []
        ccdl = []
        bandl = []
        mjdl = []
        errl = []
        fluxl= []
        ferrl = []
        magzl = []
        magl = []
        catl = []
        print('original size: ' + str(len(ralist)))
#    print('exposures**************' + str(sorted(exposures)))
#    print(sorted(set(explist)))
        emptyCount = 0
        fullCount = 0
        time0 = time.time()
        for exp in exposures:
            try:
                indList = expDict[int(exp)]
            except KeyError:
                #print('exp has no detections: ' + str(exp))
                emptyCount += 1
                continue
            fullCount += 1
            for x in indList:
                if(int(errlist[x]) > 1./3600):
                    if(catlist[x] in uniqueDict):
                        uniqueDict[catlist[x]] += 1
                    else:
                        uniqueDict[catlist[x]] = 1
                    continue
                if(ralist[x] is None or declist[x] is None or
                    explist[x] is None or ccdlist[x] is None or
                    bandlist[x] is None or mjdlist[x] is None or
                    errlist[x] is None or fluxlist[x] is None or
                    ferrlist[x] is None or magzlist[x] is None or
                    maglist[x] is None or catlist[x] is None):
                    print(x)
                    continue
                if(catlist[x] in uniqueDict):
                    uniqueDict[catlist[x]] += 1
                else:
                    uniqueDict[catlist[x]] = 1
                    
                ral.append(ralist[x])
                decl.append(declist[x])
                expl.append(explist[x])
                ccdl.append(ccdlist[x])
                bandl.append(bandlist[x])
                mjdl.append(mjdlist[x])
                errl.append(errlist[x])
                fluxl.append(fluxlist[x])
                ferrl.append(ferrlist[x])
                magzl.append(magzlist[x])
                magl.append(maglist[x])
                catl.append(catlist[x])
        #   else:
        #       print(str(ralist[x]) + ',' + str(declist[x]))
        time1 = time.time()
        print('time to insert detections: ' + str(time1-time0))
        print('empty exposures: ' + str(emptyCount))
        print('nonempty exposures: ' + str(fullCount))
        print('final list: ' + str(len(ral))) 
        totalLength = len(ral)
        assert(len(decl) == totalLength)
        assert(len(expl) == totalLength)
        assert(len(ccdl) == totalLength)
        assert(len(bandl) == totalLength)
        assert(len(mjdl) == totalLength)
        assert(len(errl) == totalLength)
        assert(len(fluxl) == totalLength)
        assert(len(ferrl) == totalLength)
        assert(len(magzl) == totalLength)
        assert(len(magl) == totalLength)
        assert(len(catl) == totalLength)

        outTable = Table([ral, decl, expl, ccdl, bandl, mjdl, errl, 
                         fluxl, ferrl, magzl, magl, catl],
                    names=('RA', 'DEC', 'EXPNUM', 'CCDNUM', 'BAND', 'MJD_OBS',
                            'ERRAWIN_WORLD', 'FLUX_AUTO', 'FLUXERR_AUTO', 
                            'MAG_ZERO', 'MAG', 'CATALOG_ID'),
                    dtype=('f8', 'f8', 'i4', 'i2', 'S5', 'f8', 
                            'f8', 'f8', 'f8', 'f8', 'f8', 'i8'))
        binTable = fits.BinTableHDU(outTable)
        hdu = fits.PrimaryHDU()
        hdul = HDUList([hdu, binTable])
        print('writing to ' + outname)
        hdul.writeto(outname, overwrite=True)
        print('total time for patch ' + str(patch) +
                    ': ' + str(time.time()-time0))
        print('next patch')
        count += totalLength
    print('total number of detections: '+ str(count))
    print('total number of detections: ' + str(len(uniqueDict)))
    return


    def countTNOs(tablenames):
        count = {}
        for table in tablenames:
            tab = Table.read(table)
            catlist = tab['CATALOG_ID'].tolist()
            for cat in catlist:
                if(cat in count):
                    count[cat] += 1
                else:
                    count[cat] = 1
        return count

def getPatchExps(exptable, patchtable, patches, radius, pool=False):
    coords = [] 
    print('finding coords of patches')
    gamma = 1./30.
    dgamma = 0.05/30
    for patch in patches:
        mask = patchtable[patchtable['PATCH']==patch]
        ra = mask['RA0'].tolist()[0]
        dec = mask['DEC0'].tolist()[0]
        coords.append((patch,(ra,dec)))
        if(pool):
            years = range(4)
            for y in years:
                tdb = 13.9+y
                getPatchExp(patch, ra, dec, tdb, gamma, dgamma, radius) 
    
    if(pool):
        return
    exptable = Table.read(exptable)
    expDict = expInRadius(exptable, coords, radius)
    return expDict


def main():
    args = argparse.ArgumentParser()
    args.add_argument('-e', '--exp', help='fits table for exposure center')
    args.add_argument('-d', '--dets', help='detection file to split')
    args.add_argument('-r', '--radius')
    args.add_argument('-p', '--patch', help='table with patch centers')
    args.add_argument('-s', '--start')
    args.add_argument('-t', '--end')
    args.add_argument('-g', '--pool', action='store_true')
    args.add_argument('-x', '--skip', action='store_true')
    args = args.parse_args()
    print('************start*************')
    print('reading in ' + args.patch)
    patchCenters = Table.read(args.patch)

    radius = 3.5
    if(args.radius):
        radius = float(args.radius)
    start = 0
    end = 240
    if(args.start and args.end):
        start = int(args.start)
        end = int(args.end)
    patches = range(start, end)
    print('patches:' + str(patches))
    if(args.skip):
        patchExps = None
    else:
        patchExps = getPatchExps(args.exp, patchCenters, patches, radius, args.pool)
    print('reading in ' + args.dets)
    dets = Table.read(args.dets)
    getDetections(dets, patches, patchExps)
        

if __name__ == '__main__':
    main()
