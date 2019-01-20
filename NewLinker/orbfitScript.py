'''
This program write a script for orbit fitter
'''
import sys
import os
tnopath = os.environ['TNO_PATH']
sys.path.insert(0, tnopath)
import argparse
import pickle
import math
from GammaTPlotwStatTNOExFaster import mjdToDate

def helpFormatTime(t):
    if len(t) == 1:
        t = '0' + t
    return t

# formats time in hh:mm:ss.ss
def formatTime(h,m,s):
    h = helpFormatTime(h)
    m = helpFormatTime(m)
    if len(s.split('.')[0]) == 1:
        s = '0' + s
    return h + ':' + m + ':' + s

def degToHour(deg):
    abval = abs(deg)
    h = int(math.floor(abval / 15))
    m = int(math.floor((abval % 15) * 4))
    s = round((abval - 15 * h - m / 4.0) * 240.0, 3)
    result = formatTime(str(h), str(m), str(s))
    if deg < 0:
        result = '-' + result
    return result

def degToDMS(degree):
    abval = abs(degree)
    deg = int(abval)
    arcmin = int((abval - deg) * 60)
    arcsec = round(3600 * (abval - deg) - 60 * arcmin, 2)
    result = formatTime(str(deg), str(arcmin), str(arcsec))
    if degree < 0:
        result = '-' + result
    return result

def decimalDay(day, time):
    t = time.split(':')
    h = float(t[0])
    m = float(t[1])
    s = float(t[2])
    return str(int(day) + h / 24 + m / 1440 + s / 86400)

def scriptWriter(triplets, saveName):
    for trip in triplets:
        saveas = saveName + '_' + str(trip.dets[0].fakeid) + '.dat'
        f = open(saveas, 'w+')
        for det in trip.dets:
            date = mjdToDate(det.mjd).split('/')
            daytime = date[2].split('T')
            f.write(date[0] + ' ' + helpFormatTime(date[1]) + ' ' + decimalDay(daytime[0],daytime[1]) +
                    ' ' + degToHour(det.ra) + ' ' + degToHour(det.dec) + ' 0.1 807\n')
        f.close()




def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('pickle', nargs=1, help='path to pickle file')
    args = parser.parse_args()
    saveName = args.pickle[0].split('+')[-1].split('.')[0]
    triplets = pickle.load(open(args.pickle[0],'rb'))
    print(degToHour(-37.13241))
    scriptWriter(triplets, saveName)

if __name__ == '__main__':
    main()
