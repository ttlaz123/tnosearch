
from astropy.table import Table
import numpy as np
import argparse

def main():
    args = argparse.ArgumentParser()
    args.add_argument('-f', '--fitsTable', help='path to fits table')
    args = args.parse_args()

    tab = Table.read(args.fitsTable)
    '''
    for row in tab:
        ele = row['ELEMENTS']
        row['ELEMENTS'] = tuple([-999.9 if np.isnan(x) else x for x in ele])
    '''
    '''
    for row in tab:
        ele = row['ELEMENTS']
        for e in ele:
            if(np.isnan(e)):
                e = -999
    '''
    for row in tab:
        ele = row['ELEMENTS']
        for e in ele:
            if(np.isnan(e)):
                row = 0
                print(ele)
    tab.write(args.fitsTable, format='fits', overwrite=True)
if __name__ == '__main__':
    main()
