#!/bin/sh

# Starts a set of BLAT servers required for 3PD primer design tomcat webapplication

GFSERVER=/home/tjakobi/scratch/threepd/3PD/3CPrimerDesign/bin/gfServer
DATADIR=../WebContent/WEB-INF/data/twobit
MINMATCH=1
STEPSIZE=1
STEPSIZEBIGSEQ=2
TILESIZE=11
TILESIZEBIGSEQ=8

#$GFSERVER -tileSize=${TILESIZE} -stepSize=${STEPSIZE} -minMatch=${MINMATCH} start localhost 10100 $DATADIR/Celegans.2bit > /dev/null & 
$GFSERVER -tileSize=${TILESIZE} -stepSize=${STEPSIZE} -minMatch=${MINMATCH} start localhost 10110 $DATADIR/Ppacificus.2bit > /dev/null & 
#$GFSERVER -tileSize=${TILESIZE} -stepSize=${STEPSIZE} -minMatch=${MINMATCH} start localhost 10140 $DATADIR/Dmelanogaster.2bit > /dev/null & 
#$GFSERVER -tileSize=${TILESIZE} -stepSize=${STEPSIZE} -minMatch=${MINMATCH} start localhost 10150 $DATADIR/Scerevisiae.2bit > /dev/null &
#$GFSERVER -tileSize=${TILESIZEBIGSEQ} -stepSize=${STEPSIZEBIGSEQ} -minMatch=${MINMATCH} start localhost 10120 $DATADIR/Mmusculus.2bit > /dev/null & 
#$GFSERVER -tileSize=${TILESIZEBIGSEQ} -stepSize=${STEPSIZEBIGSEQ} -minMatch=${MINMATCH} start localhost 10130 $DATADIR/Hsapiens.2bit > /dev/null & 
