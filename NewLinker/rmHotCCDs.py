import sys
import os
tnopath = os.environ['TNO_PATH']
sys.path.insert(0, tnopath)
import argparse
import pandas as pd
import time
from LinkerLib import printPercentage
from astropy.table import Table
# removes elements with ccd
def rmCCDs(df, thresh):
    print('number of detections: ' + str(len(df)))
    #a dictionary of dictionaries
    ccdict = {}
    
    df.rename(str.lower, axis='columns', inplace=True)
    df.rename(columns={'ccdnum': 'ccd'}, inplace=True)
    df.rename(columns={'mjd_obs': 'mjd'}, inplace=True)
    df.rename(columns={'catalog_id': 'objid'}, inplace=True)
    print('Dataframe columns: ' + df.columns)

    expnums = df['expnum'].tolist()
    ccdnums = df['ccd'].tolist()
    time0 = time.time()
    #count the number of things with certain ccd
    for x in range(len(df)):
        if(x%10000 == 0):
            printPercentage(x, len(df), time.time()-time0)
        exp = expnums[x]
        ccd = ccdnums[x]
        if exp in ccdict:
            if ccd in ccdict[exp]:
                ccdict[exp][ccd] += 1
            else:
                ccdict[exp][ccd] = 1 
        else:
            ccdict[exp] = {}
    #list of (expnum, ccdnum)
    overThresh = []
    #check if any are over thresh
    for key, value in ccdict.iteritems():
        for key2, value2 in value.iteritems(): 
            if value2 > thresh:
                print('expnum:' + str(key) + ' ccdnum:' + str(key2) + ' hits:' + str(value2))
                overThresh.append((key, key2))
    #to do remove detections with corresponding overthresh
    rmList = []
    for x in range(len(df)):
        if ((expnums[x], ccdnums[x]) in overThresh):
            rmList.append(x)
    print('\nrows removed:' + str(len(rmList)))
    df = df.drop(df.index[rmList])
    print('size after remove: ' + str(len(df)))
    return df

def insertFakeid(df, truthName):
    fakeTable = Table.read(truthName, format='fits')
    objidList = fakeTable['CATALOG_ID'].tolist()
    fakeList = fakeTable['ORBITID'].tolist()
    fakeDict ={}
    for x in range(len(fakeList)):
        fakeDict[objidList[x]] = fakeList[x]
        if(fakeList[x] == -1):
            fakeDict[objidList[x]] = 0
    oldList = df['objid'].tolist()

    newFakeList = []
    for objid in oldList:
        newFakeList.append(fakeDict[int(objid)])
    print('length: ' + str(len(newFakeList)))
    df = df.assign(fakeid=(newFakeList))
    print('old length:' + str(len(df.index)))
    return df

def main():
    args = argparse.ArgumentParser()
    #pass in varRM+ file
    args.add_argument('detections', help='path to csv detection file; ' + 
                        'filename has format SNOBS_SEASON###_ML0#.csv')
    args.add_argument('-t', '--ccdThresh',
            help='threshold for disregarding ccd')
    args.add_argument('-f', '--fakeid', help='extra fake ids')
    args = args.parse_args()
    savename = args.detections.split('+')[-1].split('/')[-1].split('.')[0]
    dets = args.detections
    if(dets.split('.')[-1] == 'csv'):
        try:
            df = pd.read_csv(args.detections)
        except IOError:
            print('try fits file')
            dets = dets.split('.')[0] + '.fits'

    if(dets.split('.')[-1] == 'fits'):
        try:
            print('reading fits file: ' + dets)
            fitsTable = Table.read(dets, format='fits') 
            df = fitsTable.to_pandas() 
        except IOError:
            print('file does not exist: ' + dets)
            return
    if(not dets.split('.')[-1]=='fits' and not dets.split('.')[-1]=='csv'):
        print('not a proper extension: ' + args.detections)
        return 

    ccdThresh = 200
    if(args.ccdThresh):
        ccdThresh = int(args.ccdThresh)
    print('ccdthresh' + str(ccdThresh))
    df = rmCCDs(df, ccdThresh)
    if(args.fakeid):
        df = insertFakeid(df, args.fakeid)
    ext = 'clean'
    if(args.detections.split('/')[-1].split('+')[0] == 'grow'):
        ext = 'grow'
    df.to_csv(ext + '+' + savename + '.csv', index=False)

if __name__=='__main__':
    main()
