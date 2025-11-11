#!/bin/sh
export OMP_NUM_THREADS=80

piks=/home/xq/Softwares/invKS/src
nohup ${piks}/iks.x > log 2>&1 &
