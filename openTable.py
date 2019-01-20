from astropy.table import Table
from astropy.io import fits
import argparse

def main():
    args = argparse.ArgumentParser()
    args.add_argument('-t', '--start', help='start trackid')
    args.add_argument('-e', '--end', help='end trackid')
    args.add_argument('table', help='path to fits table')
    args = args.parse_args()
    tab = Table.read(args.table)
    if(args.start and args.end):
        tabmask1 = (tab['ORBITID'] >= int(args.start)) 
        tabmask2 = (tab['ORBITID'] <= int(args.end))
        tabmask = (tabmask1 & tabmask2)
        rel = tab[tabmask]
    else:
        rel = tab
    data = fits.open(args.table)
    '''
    print(data[1].header['RA0'])
    print(data[1].header['DEC0'])
    print(data[1].header['MJD0'])
    '''
    print(rel)
    print(tab.info)
    print(tab)
    rel.write('sampleOrbits.fits', format='fits', overwrite=True)
if __name__ == '__main__':
    main()
