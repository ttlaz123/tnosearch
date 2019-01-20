import sys
import os
tnopath = os.environ['TNO_PATH']
sys.path.insert(0, tnopath)
import argparse
import pandas as pd
import time
import numpy as np
from LinkerLib import printPercentage

'''
converts ra and dec from units of radians to degrees
'''
def convertToDegrees(ra, dec):
    ra = float(ra)
    dec = float(dec)
    newra = ra/np.pi*180.0
    newdec = dec/np.pi*180.0
    assert(newdec < 90 and newdec > -90)
    if(newra > 180):
        newra -=360
    return newra, newdec
'''
converts all ra's and dec's in dataframe from radians to degrees
'''
def radiansToDegrees(df):
    df.rename(str.lower, axis='columns', inplace=True)
    raList = df['ra'].tolist()
    decList = df['dec'].tolist()
    raDegrees = []
    decDegrees = []
    for x in range(len(df)):
        rad, decd = convertToDegrees(raList[x], decList[x])
        raDegrees.append(rad)
        decDegrees.append(decd)
    df['ra'] = raDegrees
    df['dec'] = decDegrees
    return df


# removes elements under ml thresh
def rmML(df, thresh):
    mlthresh = int(thresh)/10.0
    print('number of detections: ' + str(len(df)))
    
    df.rename(str.lower, axis='columns', inplace=True)
    print('Dataframe columns: ' + df.columns)

    mls = df['ml_score'].tolist()
    time0 = time.time()
    rmList = []
    #remove everything with ml less than certain range
    for x in range(len(df)):
        printPercentage(x, len(df), time.time()-time0)
        if (float(mls[x]) < mlthresh):
            rmList.append(x)
    print('\nrows removed:' + str(len(rmList)))
    df = df.drop(df.index[rmList])
    print('size after remove: ' + str(len(df)))
    return df

def main():
    args = argparse.ArgumentParser()
    #pass in varRM+ file
    args.add_argument('detections', help='path to csv detection file; ' + 
                        'filename has format SNOBS_SEASON###_ML0#.csv')
    args.add_argument('-m', '--mlThresh',
            help='threshold for disregarding detection')
    args.add_argument('-r', '--rad', action='store_true', help='whether to convert to degrees from rad')
    args = args.parse_args()
    savename = args.detections.split('+')[-1].split('/')[-1].split('.')[0]
    parts = savename.split('_')
    objtype = parts[0]
    season = parts[1]
    mlThresh = '06'
    if(args.mlThresh):
        mlThresh = args.mlThresh
        savename = objtype + '_' + season + '_ML' + mlThresh
    else:
        savename = savename
    df = pd.read_csv(args.detections)
    if(args.mlThresh):
        df = rmML(df, mlThresh)
    if(args.rad):
        df = radiansToDegrees(df)
    df.to_csv(savename + '.csv', index=False)

if __name__=='__main__':
    main()
