#!/bin/bash

starttime=`date +%s`
source /global/common/software/dessn/edison/tno/setup.sh
source /global/common/software/dessn/edison/destnos/setup_dg.sh
###############################################################
#FILEPATH EDITS
###############################################################
#module load PrgEnv-gnu
#module swap PrgEnv-intel PrgEnv-gnu
#module load gcc/7.3.0
#module swap gcc/4.9.3 gcc/7.3.0
#TNO_DIR=/global/project/projectdirs/dessn/diffim/TNOsearch/
FILES=/scratch2/scratchdirs/liuto/tnosearch/
TNO_DIR=${FILES}
#FILES=/scratch3/scratchdirs/masao/tno/
FILES_DIR=${FILES}Data/
SCRIPTS_DIR=${TNO_DIR}NewLinker/
#detDir=${TNO_DIR}TTcsvDetectionFiles/
detDir=${FILES_DIR}TTcsvDetectionFiles/

################################################
# DO NOT EDIT BELOW THIS LINE
#################################################


script=$1
objtype=$2
season=$3
ml=$4
chunk=$5
optionalParam=$6

if [ ! ${#ml} -eq 2 ]
then
    echo "not proper ml error: ml is ${ml}"
    exit
elif [ ! ${#season} -eq 3 ] 
then 
    echo "not proper season error: season is ${season}"
    exit
elif [ ${objtype} != 'SNALL' ] && [ ${objtype} != 'SNFAKE' ] && [ ${objtype} != 'SNOBS' ] && [ ${objtype} != 'red_SNALL' ] && [ ${objtype} != 'SNFIT' ]
then
    echo "not proper objtype error: objtype is ${objtype}"
    exit
fi

export TNO_PATH=$SCRIPTS_DIR

startDir=${FILES_DIR}${objtype}_SEASON${season}_ML${ml}_start/
chunkDir=${FILES_DIR}${objtype}_SEASON${season}_ML${ml}_chunks/
siftChunkDir=${FILES_DIR}${objtype}_SEASON${season}_ML${ml}_siftedChunks/
growChunkDir=${FILES_DIR}${objtype}_SEASON${season}_ML${ml}_growChunks/
splitChunkDir=${FILES_DIR}${objtype}_SEASON${season}_ML${ml}_splitChunks/
extractChunkDir=${FILES_DIR}${objtype}_SEASON${season}_ML${ml}_extractChunks/
finalDir=${FILES_DIR}${objtype}_SEASON${season}_ML${ml}_finalStretch/

regionsFile=${startDir}regions+${objtype}_SEASON${season}_ML${ml}.pickle
linkFile=${startDir}detectionLinks+${objtype}_SEASON${season}_ML${ml}.pickle
chunkFile=${chunkDir}chunk${chunk}+${objtype}_SEASON${season}_ML${ml}.pickle
goodchunkFile=${siftChunkDir}chunk${chunk}+goodtriplets+${objtype}_SEASON${season}_ML${ml}.pickle
orbitFile=${siftChunkDir}chunk${chunk}+${objtype}_SEASON${season}_ML${ml}.orbit
growFile=${growChunkDir}chunk${chunk}+crossCampaignTriplets+${objtype}_SEASON${season}_ML${ml}.pickle
splitFile=${splitChunkDir}splitchunk${chunk}+crossCampaignTriplets+${objtype}_SEASON${season}_ML${ml}.pickle
extractFile=${extractChunkDir}splitchunk${chunk}+extractedTriplets+${objtype}_SEASON${season}_ML${ml}.pickle
combFile=${finalDir}chunkCombined+${objtype}_SEASON${season}_ML${ml}.pickle
mergeFile=${finalDir}merged+chunkCombined+${objtype}_SEASON${season}_ML${ml}.pickle
goodFile=${finalDir}merged+goodtriplets+${objtype}_SEASON${season}_ML${ml}.pickle
orbFile=${finalDir}merged+${objtype}_SEASON${season}_ML${ml}.orbit

orgDetFile=${detDir}${objtype}_SEASON${season}_ML${ml}.csv
orgGrowFile=${detDir}grow+${objtype}_SEASON${season}_ML${ml}.csv
detFile=${startDir}clean+${objtype}_SEASON${season}_ML${ml}.csv
growFile=${startDir}grow+${objtype}_SEASON${season}_ML${ml}.csv
expFile=${TNO_DIR}exposures_true_griz.dat

if [ $script == 'rmccd' ]
then
    pyScript=rmHotCCDs.py
    pythonPath=${SCRIPTS_DIR}${pyScript}
    if [ ${#optionalParam} -gt 0 ] 
    then
        cmd="${pythonPath} ${orgDetFile} -t ${optionalParam}"
    else
        cmd="${pythonPath} ${orgDetFile}"
    fi
    toDir=$startDir

elif [ $script == 'conv' ]
then 
    pyScript=rmHotCCDs.py
    pythonPath=${SCRIPTS_DIR}${pyScript}
    cmd="${pythonPath} ${orgGrowFile}"
    toDir=$startDir

elif [ $script == 'part' ] 
then
    pyScript=partitionSeason.py
    pythonPath=${SCRIPTS_DIR}${pyScript}

    cmd="${pythonPath} ${detFile}"
    toDir=$startDir

elif [ $script == 'linkdet' ] 
then
    pyScript=linkDetections.py
    pythonPath=${SCRIPTS_DIR}${pyScript}
    if [ ! ${#optionalParam} -gt 0 ]
    then
        cmd="${pythonPath} ${regionsFile} ${detFile}"
    else
        cmd="${pythonPath} ${regionsFile} ${detFile} -l ${optionalParam}"
    fi
    toDir=$startDir

elif [ $script == 'linkpair' ]
then
    pyScript=linkPairs.py
    pythonPath=${SCRIPTS_DIR}${pyScript}
    if [ ${#optionalParam} -gt 0 ] 
    then
        cmd="${pythonPath} ${linkFile} -n ${optionalParam}"
    else
        cmd="${pythonPath} ${linkFile} -n 100000"
    fi
    if [[ $chunk =~ ^-?[0-9]+$ ]]; then
        cmd="$cmd -c $chunkFile"
    fi
    toDir=$chunkDir
   
elif [ $script == 'sift' ] 
then
    if [ ${#chunk} -lt 6 ]
    then
        echo "not proper chunk: ${chunk}"
        exit
    fi
    pyScript=siftTriplets.py
    pythonPath=${SCRIPTS_DIR}${pyScript}

    cmd="${pythonPath} ${chunkFile} ${detFile}"
    toDir=$siftChunkDir

elif [ $script == 'grow' ] 
then
    if [ ${#chunk} -lt 6 ]
    then
        echo "not proper chunk: ${chunk}"
        exit
    fi
    pyScript=growTriplets.py
    pythonPath=${SCRIPTS_DIR}${pyScript}
    
    if [ ${#optionalParam} -gt 0 ] 
    then 
        cmd="${pythonPath} ${goodchunkFile} ${growFile} ${orbitFile}"  
    else
        cmd="${pythonPath} ${goodchunkFile} ${growFile} ${orbitFile} -w"
    fi
    toDir=$growChunkDir

elif [ $script == 'split' ] 
then
    pyScript=splitTriplets.py
    pythonPath=${SCRIPTS_DIR}${pyScript}
    if [ ${#optionalParam} -gt 0 ] 
    then
        cmd="${pythonPath} -f ${growChunkDir} -n ${optionalParam} -p crossCampaignTriplets"
    else
        cmd="${pythonPath} -f ${growChunkDir} -p crossCampaignTriplets"
    fi
    toDir=$splitChunkDir
    
elif [ $script == 'extract' ] 
then
    if [ ${#chunk} -lt 6 ]
    then
        echo "not proper chunk: ${chunk}"
        exit
    fi
    pyScript=extractTriplets.py
    pythonPath=${SCRIPTS_DIR}${pyScript}

    cmd="${pythonPath} ${splitFile} ${growFile}"
    toDir=$extractChunkDir

elif [ $script == 'combine' ] 
then
    pyScript=combineChunks.py
    pythonPath=${SCRIPTS_DIR}${pyScript}

    cmd="${pythonPath} -f ${extractChunkDir} -d ${growFile}"
    toDir=$finalDir

elif [ $script == 'merge' ] 
then
    pyScript=mergeTrips.py
    pythonPath=${SCRIPTS_DIR}${pyScript}

    cmd="${pythonPath} ${combFile}"
    toDir=$finalDir

elif [ $script == 'element' ]
then
    pyScript=siftTriplets.py
    pythonPath=${SCRIPTS_DIR}${pyScript}

    cmd="${pythonPath} ${mergeFile} ${growFile} -o"
    toDir=$finalDir    

elif [ $script == 'final' ]
then
    pyScript=finalConvert.py
    pythonPath=${SCRIPTS_DIR}${pyScript}
  
    cmd="${pythonPath} -t ${goodFile} -r ${orbFile}"
    toDir=$finalDir

elif [ $script == 'graph' ]
then
    pyScript=graphTripPath.py                                                                       
    pythonPath=${SCRIPTS_DIR}${pyScript}

    cmd="${pythonPath} ${goodFile}"
    toDir=$finalDir

else
    echo "not a proper script in pipeline"
    exit 1
fi


echo $cmd

mkdir -p $toDir
cd $toDir
python $cmd 

endtime=`date +%s`
echo $endtime $starttime | awk '{printf "Runtime = %d seconds\n", $1-$2}'
echo "JOB FINISHED"
