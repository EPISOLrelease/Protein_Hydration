#!/usr/bin/env bash

#set -e

echo "=== HTCondor Job Information ==="
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "System: $(uname -spo)"
echo "_CONDOR_JOB_IWD: $_CONDOR_JOB_IWD"
echo "Cluster: $CLUSTER"
echo "Process: $PROCESS"
echo "RunningOn: $RUNNINGON"
echo "User: $(whoami)" # Fails 
#echo "Working Directory: $(pwd)"
echo "============================="
echo ""

echo "============================="
echo "PROGRAM RUNNING WITH " ${PDB} ${RHO}
echo "============================="
python3 run_rmsd_rho_centroid.py ${PDB} ${RHO}
echo "============================= FILES"

#ls -la
echo "============================= END FILES"
echo "removing files you dont want"

echo "KYRIE ELESON"
