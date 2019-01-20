import os
import sys
import pickle
import argparse
import time

TNO_PATH='/scratch2/scratchdirs/liuto/NewLinkerTester/'
Other_Path='/global/project/projectdirs/dessn/diffim/TNOsearch/NewLinker/'
sys.path.insert(0, Other_Path)

import LinkerLib as LL
from LinkerLib import printPercentage

def submitJobsChunks(chunks):
    for c in range(chunks):
        os.system('sbatch ' + 'RUN_growTrips.sl ' + '{0:04d}'.format(c))

def submitJobsLink(seasons, objtype, run=False, overwrite=True):
    for c in seasons:
        os.system('sbatch ' + ' -J link' + str(c) + 
            ' RUN_linkPairs.sl ' + str(c) + ' ' + objtype)

def submitJobsSift(season, objType, end, start=0, trun=True, overwrite=False, ml='06', submitSize = 24):
    c = season
    end = int(end)
    start = int(start)
    lists = []
    print('start: ' + str(start))
    print('end: '+ str(end))
    for x in range(start, end):
        
        fileG = 'chunk{0:06d}+'.format(x) + objType + '_SEASON' + str(season)
        fileG = fileG + '_ML06.pickle'
        pathG = TNO_PATH+objType + '_SEASON' + str(c) + '_ML06_chunks/'
        filename = 'chunk{0:06d}+'.format(x)+'goodtriplets+' + objType + '_SEASON' + str(season)
        filename = filename + '_ML06.pickle'
        path = TNO_PATH+objType + '_SEASON' + str(c) + '_ML06_siftedChunks/'
        if((not os.path.isfile(path+filename) or overwrite) and os.path.isfile(pathG+fileG)):
            #if(run):
            lists.append(x)
        if(len(lists)%submitSize == 0 and not len(lists)==0):
            chunks = ''
            for y in lists:
                chunks = chunks + ' {0:06d}'.format(y) 
            run = TNO_PATH + 'RUN_siftTrips.sl ' + objType + ' ' + str(c) + ' ' + ml +  chunks
            cmd = ('sbatch -J extr' + str(c) + ' ' + run) 
            #cmd = 'sh ' + run
            print(run)
            if(trun):
                os.system(cmd)
            del lists[:]
    chunks = ''
    print('*******************')
    for y in lists:
        chunks = chunks + ' {0:06d}'.format(y) 
    run = TNO_PATH + 'RUN_siftTrips.sl ' + objType + ' ' + str(c) + ' ' + ml +  chunks

    cmd = ('sbatch -J sift' + str(c) + ' ' + run)
    #cmd = 'sh ' + run
    if(not chunks == ''):
        print(cmd)
        if(trun):
            os.system(cmd)


def submitJobsSift1(season, objtype, chunkNum, start=0, run=True, overwrite=False,ml='06'):
    print('running sift stage:')
    for x in range(start, chunkNum):
        filename = ('chunk{0:06d}'.format(x)+
                    '+goodtriplets+'+
                    objtype+'_SEASON' + str(season))
        filename = filename + '_ML'+ml+'.pickle'
        path = TNO_PATH+objtype+'_SEASON' + str(season) + '_ML'+ml+'_siftedChunks/'

        fileG = ('chunk{0:06d}'.format(x)+ '+'+ objtype+
                '_SEASON' + str(season))
        fileG = fileG + '_ML'+ml+'.pickle'
        pathG = TNO_PATH+objtype+'_SEASON' + str(season) + '_ML'+ml+'_chunks/'
        if((overwrite or not os.path.isfile(path + filename) )
            and os.path.isfile(pathG+fileG)):
            chunk = '{0:06d}'.format(x)
           
            cmd = ('sbatch -J sift' + str(x) + ' ' + TNO_PATH + 
                    'RUN_siftTrips.sl ' +  objtype+ ' ' + str(season) + 
                    ' ' + ml + ' ' + chunk)
            cmd2 = ('sh ' + TNO_PATH + 'RUN_siftTrips.sl ' + objtype + 
                    ' ' + str(season) + ' ' + ml + ' ' + chunk )
            print(cmd)
            if(run):
                os.system(cmd)

def submitJobsGrow1(season, objtype, end, start=0, run=False, overwrite=False, ml='06', opt=True):
    print('running grow stage:')
    c = season
    missingCounter = 0
    end = int(end)
    start = int(start)
    for x in range(start, end):
        
        filename = ('chunk{0:06d}'.format(x)+
                    '+crossCampaignTriplets+'+
                    objtype+'_SEASON' + str(season))
        filename = filename + '_ML'+ml+'.pickle'
        path = TNO_PATH+objtype+'_SEASON' + str(c) + '_ML'+ml+'_growChunks/'

        fileG = ('chunk{0:06d}'.format(x)+
                '+goodtriplets+'+objtype+
                '_SEASON' + str(season))
        fileG = fileG + '_ML'+ml+'.pickle'
        pathG = TNO_PATH+objtype+'_SEASON' + str(c) + '_ML'+ml+'_siftedChunks/'
        if((not os.path.isfile(path+filename) or overwrite)
            and os.path.isfile(pathG+fileG)):
            run_cmd =  (TNO_PATH + 'RUN_growTrips.sl ' + 
                         objtype + ' ' + str(c) + ' ' + ml)
            if(opt):             
                run_cmd = run_cmd + ' {0:06d}'.format(x)
            else:
                run_cmd = run_cmd + ' x {0:06d}'.format(x)
            cmd = ('sbatch -J grow' + str(x) + ' ' + run_cmd)
            cmd2 = 'sh ' + run_cmd
            print(cmd)
            if(run):
                os.system(cmd)

            
    return missingCounter

def submitJobsGrow(season, objType, end, start=0, run1=True, overwrite=False, ml='06', submitSize = 24):
    c = season
    end = int(end)
    start = int(start)
    lists = []

    print('start: ' + str(start))
    print('end: '+ str(end))
    for x in range(start, end):
        
        filename = 'chunk{0:06d}'.format(x)+'+goodtriplets+' + objType + '_SEASON' + str(season)
        filename = filename + '_ML06.pickle'
        path = TNO_PATH+objType + '_SEASON' + str(c) + '_ML06_siftedChunks/'
        fileG = 'chunk{0:06d}'.format(x)+'+crossCampaignTriplets+' + objType + '_SEASON' + str(season)
        fileG = fileG + '_ML06.pickle'
        pathG = TNO_PATH+objType + '_SEASON' + str(c) + '_ML06_growChunks/'
        if((not os.path.isfile(pathG+fileG) or overwrite) and os.path.isfile(path+filename)):
            #if(run):
            lists.append(x)
        if(len(lists)%submitSize == 0 and not len(lists)==0):
            chunks = ''
            for y in lists:
                chunks = chunks + ' {0:06d}'.format(y) 
            run = TNO_PATH + 'RUN_growTrips.sl ' + objType + ' ' + str(c) + ' ' + ml +  chunks
            cmd = ('sbatch -J grow' + str(c) + ' ' + run) 
            #cmd = 'sh ' + run
            print(cmd)
            if(run1):
                os.system(cmd)
            del lists[:]
    chunks = ''
    for y in lists:
        chunks = chunks + ' {0:06d}'.format(y) 
    run = TNO_PATH + 'RUN_growTrips.sl ' + objType + ' ' + str(c) + ' ' + ml +  chunks

    cmd = ('sbatch -J grow' + str(c) + ' ' + run)
    #cmd = 'sh ' + run
    if(not chunks == ''):
        print(cmd)
        if(run1):
            os.system(cmd)

def submitJobsExtract(season, objType, end, start=0, run1=True, overwrite=False, ml='06', submitSize=24):
    c = season
    end = int(end)
    start = int(start)
    lists = []
    print('start: ' + str(start))
    print('end: '+ str(end))
    for x in range(start, end):
        
        filename = 'splitchunk{0:06d}'.format(x)+'+extractedTriplets+' + objType + '_SEASON' + str(season)
        filename = filename + '_ML06.pickle'
        path = TNO_PATH+objType + '_SEASON' + str(c) + '_ML06_extractChunks/'
        fileG = 'splitchunk{0:06d}'.format(x)+'+crossCampaignTriplets+' + objType + '_SEASON' + str(season)
        fileG = fileG + '_ML06.pickle'
        pathG = TNO_PATH+objType + '_SEASON' + str(c) + '_ML06_splitChunks/'
        if((not os.path.isfile(path+filename) or overwrite) and os.path.isfile(pathG+fileG)):
            #if(run):
            lists.append(x)
        if(len(lists)%submitSize == 0 and not len(lists)==0):
            chunks = ''
            for y in lists:
                chunks = chunks + ' {0:06d}'.format(y) 
            run = TNO_PATH + 'RUN_extractTrips.sl ' + objType + ' ' + str(c) + ' ' + ml +  chunks
            cmd = ('sbatch -J extr' + str(c) + ' ' + run) 
            #cmd = 'sh ' + run
            print(cmd)
            if(run1):
                os.system(cmd)
            del lists[:]
    chunks = ''
    for y in lists:
        chunks = chunks + ' {0:06d}'.format(y) 
    run = TNO_PATH + 'RUN_extractTrips.sl ' + objType + ' ' + str(c) + ' ' + ml +  chunks

    cmd = ('sbatch -J extr' + str(c) + ' ' + run)
    #cmd = 'sh ' + run
    if(not chunks == ''):
        print(cmd)
        if(run1):
            os.system(cmd)

def submitJobsExtract1(season, objtype, end, start=0, run=False, overwrite=False, ml='06'):
    print('running extract stage')
    end = int(end)
    start = int(start)
    counter = 0
    for x in range(start, end, 1):
        filename = ('splitchunk'+'{0:06d}'.format(x)+
                    '+extractedTriplets+' + objtype + '_SEASON' + str(season))
        filename = filename + '_ML'+ml+'.pickle'
        path = TNO_PATH+objtype + '_SEASON' + str(season) + '_ML'+ml+'_extractChunks/'
        fileG = ('splitchunk'+'{0:06d}'.format(x)+
                '+crossCampaignTriplets+' + objtype + '_SEASON' + str(season))
        fileG = fileG + '_ML'+ml+'.pickle'
        pathG = TNO_PATH+objtype + '_SEASON' + str(season) + '_ML'+ml+'_splitChunks/'
        if((overwrite or not os.path.isfile(path+filename))
                and os.path.isfile(pathG+fileG)):
            counter+=1
            print(path+filename)
            runj = (TNO_PATH + 'RUN_extractTrips.sl ' + objtype +
                ' ' + str(season) + ' ' + ml + ' {0:06d}'.format(x)) 
            cmd = 'sbatch -J s' + str(season) + 'c' + str(x) + ' ' + runj
            cmd2 = 'sh ' + runj
            print(cmd)
            if(run):
                os.system(cmd)
    print('number: ' + str(counter))

def main():
    args = argparse.ArgumentParser()
    args.add_argument('-s', '--season', help='season to submit')
    args.add_argument('-e', '--end', help='end chunk')
    args.add_argument('-t', '--start', help='start chunk')
    args.add_argument('-j', '--job', help='extract, grow, sift, link')
    args.add_argument('-d', '--folder', help='folder to chunks')
    args.add_argument('-b', '--objtype', help='object type')
    args.add_argument('-r', '--run', action='store_true', help='whether to run')
    args.add_argument('-w', '--overwrite', action='store_true', 
                help='whether to overwrite existing files')
    args.add_argument('-m', '--ml', help='machine learning score')
    args.add_argument('-o', '--blah', action='store_false', help='random extra parameter')
    args.add_argument('-c', '--submitSize', help='number of things to run per job')
    args = args.parse_args()
    if(args.folder):
        submitJobsExtract1(args.folder)
        return
    #if(args.folder):
        #moveChunks(args.folder) 
        #countTriplets(args.folder)
        #rmChunks(args.folder, int(args.start), int(args.end), args.season)
        #rmBadGrows(args.folder)
    #else:
    submitSize = 24
    if(args.submitSize):
        submitSize = int(args.submitSize)
    print('run:')
    print(args.run)
    ml = '06'
    if(args.ml):
        ml = args.ml
    job = args.job 
    if(job == 'extract'):
        print('extracting:')
        submitJobsExtract(args.season, args.objtype, 
                    args.end, args.start, args.run, args.overwrite, ml, submitSize)
    elif(job == 'grow'):
        submitJobsGrow(args.season, args.objtype, args.end, 
                    args.start, args.run, args.overwrite, ml, submitSize)
    elif(job == 'sift'):
        submitJobsSift(args.season, args.objtype, int(args.end), 
                        int(args.start),args.run, args.overwrite, ml,submitSize)
    elif(job == 'link'):
        submitJobsLink(args.season, args.objtype, args.overwrite)
    else:
        print('not an accepted job: ' + args.job)

if __name__=='__main__':
    main()

