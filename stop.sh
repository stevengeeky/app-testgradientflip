#!/bin/bash

if [ -f finished ]; then
    echo "job has already finished"
    exit 1
fi

if [ -f jobid ]; then
    jobid=`cat jobid`
    echo "running qdel $jobid"
    qdel $jobid
fi

if [ -f pid ]; then
    pid=`cat pid`
    echo "running kill" 
    #kill $pid
    kill -- -$(ps -o pgid= $pid | grep -o [0-9]*)
fi
